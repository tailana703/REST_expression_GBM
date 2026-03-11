############################################
### gene names conversion functions

###function to convert gene IDs to ensembl Ids###
SYMBOLtoENS <- function(genes = character()) {
  require("EnsDb.Hsapiens.v86")
  mapIds <- mapIds(EnsDb.Hsapiens.v86, keys = genes, keytype = "SYMBOL")
  mapIds_df <- as.data.frame(mapIds)
  mapIds_df$mapIds
}

###opposite conversion###
ENStoSYMBOL <- function(genes = character()) {
  require("org.Hs.eg.db")
  mapIds <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="SYMBOL")
  mapIds_df <- as.data.frame(mapIds)
  mapIds_df$mapIds
}


### original query on TCGA-GBM RNA-Seq dataset
library(TCGAbiolinks)
library(SummarizedExperiment)
query_GBM <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM-UQ",
                      legacy = FALSE,
                      sample.type = c("Primary Tumor", "Recurrent Tumor"),
                      access = "open")

GDCdownload(query_GBM)
###.rda file was created on 2/13/2021 and can be found in repo###
e_gbm <- load(file = "TCGA_exp.gbm.rda")
colData(data)
assays(data)

### cleaning and filtering expression data
exp_Prep <- TCGAanalyze_Preprocessing(object = data, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - FPKM-UQ")                      

exp_Filt <- TCGAanalyze_Filtering(tabDF = exp_Prep,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

### preparing expression matrix and locating REST gene
t_exp_matrix <- t(exp_Filt)
expr <- as.data.frame(t_exp_matrix)
dim(expr)

# Ensure numeric
expr[] <- lapply(expr, as.numeric)
REST_position <- which(colnames(expr)  == "ENSG00000084093") #locating REST gene

### uploading canonical REST target genes (obtained via integration of ChIP-Seq and RNA-Seq data [PMID:38609948])
## The full list is located in this file: REST_target_genes.csv; REST_genes_qPCR_validated.csv contains experimentally-validated genes that respond to pharmacological REST degradation
canonical_REST_targets <- as.data.frame(read.csv("REST_genes_qPCR_validated.csv"))
canonical_REST_targets <- canonical_REST_targets[,2]

# Ensure data is clean
canonical_REST_targets <- canonical_REST_targets[!is.na(canonical_REST_targets) & canonical_REST_targets != ""]

# converting to ensembl simbols
RESTTargetGenes <- SYMBOLtoENS(canonical_REST_targets)

### feature engineering = REST activity z-scores
# Z-score each gene across samples
REST_targets_z <- scale(expr[, RESTTargetGenes])

# REST represses targets → invert sign
REST_activity <- -rowMeans(REST_targets_z, na.rm = TRUE)

# Add to metadata
df <- colData(data)
df$REST_activity <- REST_activity[rownames(df)]

###sanity checks

# REST activity vs REST mRNA (supposedly, is not too strong)
REST_mRNA <- expr[,REST_position]
names(REST_mRNA) <- rownames(expr)
REST_protein <- df$REST_activity

REST_mRNA_aligned <- REST_mRNA[names(REST_protein)]
cor(REST_protein, REST_mRNA_aligned, use = "complete.obs")

# REST activity by subtype
library(ggplot2)
ggplot(df, aes(x = subtype_Original.Subtype, y = REST_activity, fill = subtype_Original.Subtype)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(title = "REST activity across GBM subtypes",
       y = "REST activity (mean z-score)") +
  scale_y_continuous(limits = c(-3, 1.5))

## Fitting GLM models for REST activity to find shared clinical characteristics significantly associated with REST 
require(tidyverse)
require(caret)
require(MASS)
require(nestfs)
require(data.table)
require(dplyr)

# Checking normality of the outcome variable
shapiro.test(REST_activity)
hist(REST_activity, breaks = 30)
qqnorm(REST_activity); qqline(REST_activity, col = "red") ## distribution is not normal; hence, GLM is more appropriate

## Cleaning missing values in GBM subtype
clean_df_no_na_subtype <- as.data.frame(df) %>% 
  dplyr::filter(subtype_Original.Subtype == "Classical" |
                subtype_Original.Subtype == "Mesenchymal" |
                subtype_Original.Subtype == "Neural" |
                 subtype_Original.Subtype == "Proneural") %>% 
  dplyr::select(REST_activity, subtype_Original.Subtype,
         age_at_diagnosis, alcohol_history, gender,
         subtype_Mutation.Count, subtype_IDH.status, subtype_X1p.19q.codeletion,
         subtype_MGMT.promoter.status, subtype_ATRX.status, subtype_BRAF.V600E.status)

# drop unused levels just in case
clean_df_no_na_subtype$subtype_Original.Subtype <-
  droplevels(clean_df_no_na_subtype$subtype_Original.Subtype)

## modeling with Lasso regression (GLM)
  
library(caret)
library(glmnet)

set.seed(1)

# predictors + outcome in a single data frame
df <- clean_df_no_na_subtype
df <- df[complete.cases(df), ]

# find factors with <2 levels
one_level <- sapply(df, function(x) is.factor(x) && nlevels(x) < 2)
df <- df[, !one_level]

# also drop near-zero-variance predictors (optional but helpful)
nzv <- nearZeroVar(df, saveMetrics = TRUE)
df <- df[, !nzv$nzv]

## training
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3
)

grid <- expand.grid(
  alpha = 1,                 # LASSO
  lambda = 10^seq(-4, 1, length.out = 50)
)

lasso_fit <- train(
  REST_activity ~ .,
  data = df,
  method = "glmnet",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneGrid = grid
)

lasso_fit
coef(lasso_fit$finalModel, s = lasso_fit$bestTune$lambda)

## model performance on a test set (20% of data)
set.seed(1)
idx <- createDataPartition(df$REST_activity, p = 0.8, list = FALSE)
train_df <- df[idx, ]
test_df  <- df[-idx, ]

lasso_fit <- train(
  REST_activity ~ .,
  data = train_df,
  method = "glmnet",
  trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3),
  preProcess = c("center", "scale"),
  tuneGrid = grid
)

pred <- predict(lasso_fit, newdata = test_df)
postResample(pred, test_df$REST_activity)   # RMSE, Rsquared, MAE

sd(REST_activity) ## RMSE << sd, model is better than null (predicting mean)

## subseting patients (n = 145) into low- and high-REST groups (cutoff=median)
clean_df_no_na_subtype$REST_group <- ifelse(
  clean_df_no_na_subtype$REST_activity > median(clean_df_no_na_subtype$REST_activity, na.rm = TRUE),
  "High-REST", "Low-REST"
)


## comparison of cell signaling (proteomic data) between REST-low and REST-high samples using RPPA data to identify pathways for synergetic inhibition

## uploading open RPPA data, level 3 (normalized+scaled)
RPPA <- as.data.frame(read_csv("TCGA-GBM-L3-S42.csv"))

## merge RPPA data with REST activity

###wrangling RPPA sample-IDs and mRNA-samples-IDs so they are compatible for joining
RPPA$Sample_ID <- substr(RPPA$Sample_ID,1,nchar(RPPA$Sample_ID)-15)
clean_df_no_na_subtype$Sample_ID <- substr(rownames(clean_df_no_na_subtype),1,12)

merged_RPPA <- merge(RPPA, clean_df_no_na_subtype[, c("Sample_ID", "REST_activity", "REST_group")],
                     by = "Sample_ID")

## assessing proliferation scores in every sample based on proliferative biomarkers 
prolif_markers <- c("MKI67", "PCNA", "CCNB1", "CDK1", "FOXM1")
prolif_markers <- prolif_markers[prolif_markers %in% colnames(merged_RPPA)]

merged_RPPA$Prolif_score <- rowMeans(
  (merged_RPPA[, prolif_markers]),
  na.rm = TRUE
)

## correlation between REST activity and proliferation
cor.test(merged_RPPA$REST_activity,
         merged_RPPA$Prolif_score,
         method = "spearman") #NS correlation

## pathway-level RPPA analysis - what pathways can be beneficial to target simultaneously with REST?
library(limma)

design <- model.matrix(~ REST_group, data = merged_RPPA)
fit <- lmFit(t(merged_RPPA[, 4:(ncol(RPPA))]), design)
fit <- eBayes(fit)

topTable(fit, coef = "REST_groupLow-REST")

## vis of top differentially expressed protein markers
ggplot(merged_RPPA, aes(x = REST_group, y = HER3_pY1289, col = REST_group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "HER3_pY1289",
       y = "Protein amount") 

ggplot(merged_RPPA, aes(x = REST_group, y = Bad_pS112, col = REST_group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Bad_pS112",
       y = "Protein amount") 

## Conclusion: apoptosis activators may be promising to combine with REST targeting
