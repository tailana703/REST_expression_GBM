# REST_expression_GBM
## Rationale and algorithm

Repressor element-1 silencing transcription factor (REST) is a transcriptional repressor that has been identified as an oncogenic protein driving cell proliferation and migration in various brain tumor types, including neuroblastoma, medulloblastoma, and glioblastoma [<https://pubmed.ncbi.nlm.nih.gov/7697725>; <https://pubmed.ncbi.nlm.nih.gov/23414932>; <https://pubmed.ncbi.nlm.nih.gov/38609948>]. 
However, not all the glioblastoma tumors are dependent on REST, hence predicting functional activity of this protein to identify high-REST tumors may represent a beneficial personalized treatment and/or prognosis strategy.

### Goal 
In this project, we focused on the search of transcriptomic biomarkers predicting functional activity of REST among its target genes, followed by clustering of samples into high- and low-REST to identify potential molecular pathways to be inhibited simultaneously with REST to achieve synergetic eradication of glioblastoma cells.

### Algorithm
The problem is that REST is controlled on post-translational level by upstream phosphatases and kinases, hence it is hard to predict its protein level simply via its mRNA expression. Here, we developed an algorithm to predict REST activity via expression of its target genes. The algorithm includes the following steps:
* Analysis of TCGA-GBM open RNA-Seq data and deriving top-100 genes that negatively correlate with REST mRNA level;
* Overlapping these genes with a published gene subset that have conserved RE-1 elements - DNA motifs that are canonical binding sites for REST <https://pubmed.ncbi.nlm.nih.gov/25990720>;
* The core function, *FitFSwithNestedCV()*, builds a regression model with forward selection of genes predicting expression level of every REST-target gene from step 2. Overfitting is controlled by nested cross-validation. This step may take a while.
* Next, overlapping genes (that were selected by several predictive models) are saved and used as a final panel in order to cluster patients' samples into predicted high-REST and low-REST subsets.
* Transcriptomic/proteomic data from the clusters can be further compared in order to identify signaling pathways co-occurring with REST upregulation.


## Usage

This code was written and tested using R version 4.0.5 (2021-03-31) and publicly available RNA-Seq and RPPA data from TCGA-GBM project:
<https://www.cancer.gov/ccg/research/genome-sequencing/tcga> 

## Contributing

Pull requests (especially employing different types of predictive models with superior performance) are welcome.
