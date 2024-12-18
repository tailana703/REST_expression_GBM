############################################
###Goal1: to find top REST-target genes that are negatively correlated with REST on mRNA level in GBM samples###

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


library(TCGAbiolinks)
library(SummarizedExperiment)
###querying TCGA-GBM RNA-Seq dataset###
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

###cleaning and filtering expression data###
exp_Prep <- TCGAanalyze_Preprocessing(object = data, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - FPKM-UQ")                      

exp_Filt <- TCGAanalyze_Filtering(tabDF = exp_Prep,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   

###transposing the matrix of counts###
t_exp_matrix <- t(exp_Filt)
dim(t_exp_matrix)

###obtaining correlation coefficients###
REST_position <- which(colnames(t_exp_matrix)  == "ENSG00000084093")
out_cor <- t(cor(t_exp_matrix[, REST_position], t_exp_matrix[, -REST_position]))
out_cor <- as.data.frame(out_cor)

###ranking genes###
library(data.table)
setDT(out_cor, keep.rownames = "Genes")
out_cor_ranked <- out_cor[order(out_cor$V1),]
head(out_cor_ranked, 100)
colnames(out_cor_ranked)[2] <- "Pearson r"

###obtaining p-values and q-values###
library(dplyr)
out_cor_ranked %>%
  mutate(p_value = NA) %>%
  mutate(q_value = NA)

for (i in 1:nrow(out_cor_ranked)) {
  gene_correl <- which(colnames(t_exp_matrix)  == out_cor_ranked$Genes[i])
  out_cor_ranked$p_value[i] = cor.test(t_exp_matrix[, REST_position], t_exp_matrix[, gene_correl])$p.value
}

library(WGCNA)
qqs <- qvalue(out_cor_ranked$p_value)
out_cor_ranked$q_value <- qqs$qvalues

###writing all ranked correlations to a csv file###
write.csv(as.data.frame(out_cor_ranked), "Output_TCGA_correlations_REST.csv")

###Goal2: to build regression-based predictive model to assess functional (silencing) activity of REST transcription factor based on its target genes###
###Method: forward regression with nested CV###

###feature selection: use top-100 negatively correlated genes and screen for those that have canonical/conserved RE-1 motif###
###list of genes with RE-1 motifs is from https://pubmed.ncbi.nlm.nih.gov/25990720/###

###convert ensembl top-100 REST-target gene IDs to hgnc-symbol###
feature_genes_ensembl <- out_cor_ranked$Genes[1:100]
feature_genes <- ENStoSYMBOL(feature_genes_ensembl)

###accessing published binding sites of REST in ESC cells and filtering genes with conserved RE-1 motifs###
RE1 <- as.data.frame(read_excel("Binding sites_REST_ESC_PMID25990720.xlsx", sheet = "Table S2 human (hg19)"))
RE1 <- RE1 %>%
  filter(`non-alignable, alignable conserved or alignable non-conserved` == "conserved" & is.na(`gene name`) == F)
dim(RE1) ##1252 genes left

###finding overlap between our feature genes and RE-1 element-containing genes###
feature_genes_RE1 <- intersect(RE1$`gene name`, feature_genes)

###converting back to ensembl ids###
RESTTargetGenes <- SYMBOLtoENS(feature_genes_RE1)

###EDA:exploring GBM patients metrics that correlate with REST mRNA###
REST_mRNA <- c(t_exp_matrix[, "ENSG00000084093"])
IDH_status <- data$subtype_IDH.status
GBM_subtype <- data$subtype_Transcriptome.Subtype
age <- data$age_at_diagnosis
gender <- data$gender

df <- data.frame(REST= REST_mRNA, IDH_status = IDH_status, 
                 GBM_subtype = GBM_subtype, age = age,
                 gender = gender)

mod = lm(REST~GBM_subtype+gender+age+IDH_status, data=df)
anova(mod) ## REST expression does not correlate with age, gender, and IDH status, but REST level is dependent on GBM molecular subtype###

REST_CL <- df %>%
  dplyr::filter(GBM_subtype == "CL") #classical GBM
REST_ME <- df  %>%
  dplyr::filter(GBM_subtype == "ME")  #mesenchymal GBM
REST_NE <- df %>%
  dplyr::filter(GBM_subtype == "NE")  #neural GBM
REST_PN <- df %>%
  dplyr::filter(GBM_subtype == "PN") #proneural GBM

types <- c(rep("CL", length(REST_CL$REST)), rep("ME", length(REST_ME$REST)),
           rep("NE", length(REST_NE$REST)), rep("PN", length(REST_PN$REST)))
REST_mRNA <- c(REST_CL$REST, REST_ME$REST, REST_NE$REST, REST_PN$REST)

library(ggplot2)
types_df <- data.frame(Type = types, REST_exp = REST_mRNA)
p <- ggplot(types_df, aes(x = Type, y=REST_exp, color = Type)) +
  geom_boxplot() + labs(title="REST expression (TCGA-GBM)", x=NULL, y="FPKM")
p <- p+ theme(plot.title = element_text(size = 20))
p<- p+ theme(text = element_text(size=20),
             axis.text.x = element_text(size=15),
             axis.text.y = element_text(size=15)) 
p
###based on this plot, REST expression is the highest in most aggressive GBMs: classical and mesenchymal types.

###adding GBM subtype to expression data###
t_exp_matrix <- as.data.frame(t_exp_matrix)
t_exp_data <- cbind(t_exp_matrix, df$GBM_subtype)
names(t_exp_data)[names(t_exp_data) == 'df$GBM_subtype'] <- 'GBM_subtype'

t_exp_data <- t_exp_data %>%
  dplyr::filter(is.na(GBM_subtype) == FALSE)
dim(t_exp_data) ##29 patients without GBM subtype info were removed###

###Fitting models for every REST target gene to find shared predictors###
FitFSwithNestedCV <- function(gene, data) {
  require(tidyverse)
  require(caret)
  require(MASS)
  require(nestfs)
  require(data.table)
  require(dplyr)
  
  ###Looking for top-10 correlated genes (pos/neg)###
  gene_position <- which(colnames(data)  == gene)
  cor <- t(cor(data[, gene_position], data[, -gene_position]))
  cor <- as.data.frame(cor)
  cor[,"Genes"] <- rownames(cor)
  cor_ranked <- cor[order(cor$V1),]
  cor_ranked <- na.omit(cor_ranked)
  predictors <- c(cor_ranked$Genes[1:10], cor_ranked$Genes[(nrow(cor_ranked)-9):nrow(cor_ranked)])
  data_reduced <- data[,c(gene, predictors, "GBM_subtype")]

  ###Modeling with nested CV###
  folds <- create.folds(5, nrow(data))
  nest.res <- nested.fs(as.formula(paste(gene, "~ GBM_subtype")), data= data_reduced, family=gaussian(), folds=folds,
                        min.llk.diff = 10)
  summary(nest.res)
}

###this loop may take a while.
results_mod <- list()
for (i in 1:length(RESTTargetGenes)) {
  results_mod[i] <- FitFSwithNestedCV(RESTTargetGenes[i], data = t_exp_data)
}

total_res <- unlist(results_mod)
tail(sort(table(total_res)), 10)

###saving top genes predicting REST protein level (occur more than 4 times): MAPK8IP2, ATP1A3, RUNDC3A, CDK5R1, REST, TLCD3B
max_overlap <- c("ENSG00000008735", "ENSG00000105409", "ENSG00000108309",
                 "ENSG00000176749", "ENSG00000084093", "ENSG00000149926")

###Goal3: using top-overlapping predictors, to discriminate all 169 tumor samples into REST-low and REST-high###
cor(t_exp_data[,c(max_overlap)]) ###negative correlation between REST mRNA and all these genes
exp_data_final <- t_exp_matrix[,c(max_overlap)]

###k-means clustering of samples into two categories###
cl <- kmeans(exp_data_final, centers = 2, nstart = 10)
cl$centers

library(ggplot2)
p <- ggplot(data = exp_data_final,
            aes(x= ENSG00000084093, y=ENSG00000108309, color = factor(cl$cluster))) +
  geom_point() +
  labs(title = "Clustering of TCGA-GBM samples", y = "RUNDC3A mRNA", x = "REST mRNA")
p

###Goal 4: to compare signaling pathways (proteomic data) between REST-low and REST-high samples using RPPA data to identify pathways for synergetic inhibition###
###uploading TCGA samples (mRNA expression, 6 top-predictor genes) from cBioPortal and RPPA data, separately###

RPPA <- as.data.frame(read_csv("TCGA-GBM-L3-S42.csv"))
mRNA_tested_samples <- as.data.frame(read_delim("TCGA_mRNA_tested_samples.txt"))

###filtering out missing values, e.g. samples with Not Performed (NP) expression testing
mRNA_tested_samples_filt <- mRNA_tested_samples %>%
  dplyr::filter(REST != "NP")
mRNA_tested_samples_filt <- mRNA_tested_samples_filt[,-1]
matrix <- mRNA_tested_samples_filt[,2:7]
rownames(matrix) <- mRNA_tested_samples_filt[,1]
  
###clustering these samples into high- and low-REST
cl <- kmeans(matrix, centers = 2, nstart = 10)
cl$centers
res_cl <- as.data.frame(cl$cluster)
res_cl[,"Sample_ID"] <- rownames(res_cl)
head(res_cl)
rownames(res_cl) <- NULL
colnames(res_cl) <- c("Cluster", "Sample_ID")

###wrangling RPPA sample-IDs and mRNA-samples-IDs so they are compatible for joining
RPPA$Sample_ID <- substr(RPPA$Sample_ID,1,nchar(RPPA$Sample_ID)-15)
res_cl$Sample_ID <- substr(res_cl$Sample_ID,1,nchar(res_cl$Sample_ID)-3)

###merging RPPA and gene expression data###
merged_data <- merge(RPPA, res_cl, by = "Sample_ID")
dim(merged_data)
table(merged_data$Cluster)

###Analysis of both clusters separately + calculating average protein expression for every protein
merged_split <- split(merged_data, merged_data$Cluster)
proteins <- sapply(merged_split, function(x) colMeans(x[, 4:ncol(merged_data)]))

###sorting by cluster1 (high-REST) in the order of increasing the concentration###
proteins <- as.data.frame(proteins)
colnames(proteins) <- c("High-REST", "Low-REST")
proteins_sorted <- proteins[order(proteins$`High-REST`, decreasing = F),]
proteins_sorted

###printing out ratios
ratios <- proteins_sorted$`High-REST`/proteins_sorted$`Low-REST`
names(ratios) <- rownames(proteins_sorted)

###sorting in order of decreasing ratio High/Low-REST
ratios <- ratios[order(ratios, decreasing = T)]
par(mar=c(15,4,4,2))
barplot(ratios[1:20], col = "aquamarine",
        main="High-REST associated proteins",ylab="Ratio High/Low REST",las=2)

###Conclusion: several pathways like EGFR or PAI inhibition may be promising to combine with REST targeting###

