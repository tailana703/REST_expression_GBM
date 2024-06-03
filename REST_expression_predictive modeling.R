############################################
###Goal1: to find top REST-target genes that are negatively correlated with REST on mRNA level in GBM samples###

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
###.rda file was created on 2/13/2021###
e_gbm <- load(file = "exp.gbm.rda")
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
write.csv(as.data.frame(out_cor_ranked), "Output_TCGA_ALL correlations with REST.csv")

###Goal2: to build regression-based predictive model to assess functional (silencing) activity of REST transcription factor based on its target genes###
###Method: forward regression with nested CV###

###feature selection: use top-100 negatively correlated genes and screen for those that have canonical/conserved RE-1 motif###
###list of genes with RE-1 motifs is from https://pubmed.ncbi.nlm.nih.gov/25990720/###

###convert ensembl REST-target gene IDs to hgnc-symbol###
feature_genes_ensembl <- out_cor_ranked$Genes[1:100]
require(biomaRt)
mart <- useEnsembl("ensembl",mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol"),
  filter = "ensembl_gene_id",
  values = feature_genes_ensembl,
  uniqueRows=FALSE)      
feature_genes <- annotLookup$hgnc_symbol

###overlap with RE-1 motif genes###
RE1 <- as.data.frame(read_excel("Binding sites_REST in ESC_ suppl data.xlsx", sheet = "Table S2 human (hg19)"))
RE1 <- RE1 %>%
  filter(`non-alignable, alignable conserved or alignable non-conserved` == "conserved" & is.na(`gene name`) == F)
dim(RE1) ##1252 genes left

feature_genes_RE1 <- intersect(RE1$`gene name`, feature_genes)

###converting back to ensembl ids###
feature_genes_ensembl <- character()
for (i in 1:nrow(annotLookup)) {
  if (annotLookup$hgnc_symbol[i] %in% feature_genes_RE1) {
    feature_genes_ensembl <- append(feature_genes_ensembl, annotLookup$ensembl_gene_id[i])
  }
}

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
  filter(GBM_subtype == "CL") #classical GBM
REST_ME <- df  %>%
  filter(GBM_subtype == "ME")  #mesenchymal GBM
REST_NE <- df %>%
  filter(GBM_subtype == "NE")  #neural GBM
REST_PN <- df %>%
  filter(GBM_subtype == "PN") #proneural GBM

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

###reduce expression data to genes of interest and metrics of interest###
t_exp_matrix <- as.data.frame(t_exp_matrix)
reduced_data <- t_exp_matrix %>%
  dplyr::select(all_of(feature_genes_ensembl))

reduced_data_final <- cbind(reduced_data, df[,c(1,3)])
head(reduced_data_final)

reduced_data_final <- reduced_data_final %>%
  filter(is.na(GBM_subtype) == FALSE)
dim(reduced_data_final) ##29 patients without GBM subtype info###

###Fitting models###
library(tidyverse)
library(caret)
library(MASS)
library(nestfs)

###Baseline model without CV - not possible to estimate performance as all the data is used for training###
###Init.model is REST ~ GBM_subtype
fs.res <- fs(REST ~ GBM_subtype, data= reduced_data_final, family=gaussian(),
             min.llk.diff = 2, seed =4)
summary(fs.res)

###Modeling with nested CV###
folds <- create.folds(10, nrow(reduced_data_final), seed=1)
nest.res <- nested.fs(REST ~ GBM_subtype, data= reduced_data_final, family=gaussian(), folds=folds,
                      min.llk.diff = 2, sel.crit = "both", seed = 40)
summary(nest.res)
nested.performance(nest.res)
nest.res[[1]]$final.model
nest.res[[1]]$model





