# REST_expression_GBM
## Rationale and algorithm

Repressor element-1 silencing transcription factor (REST) is a transcriptional repressor that has been identified as an oncogenic protein driving cell proliferation and migration in various brain tumor types, including neuroblastoma, medulloblastoma, and glioblastoma [<https://pubmed.ncbi.nlm.nih.gov/7697725/>; <https://pubmed.ncbi.nlm.nih.gov/23414932/>; <https://pubmed.ncbi.nlm.nih.gov/38609948/>]]. 
However, not all the glioblastoma tumors are dependent on REST, hence predicting functional activity of this protein to identify high-REST tumors may represent a beneficial personalized treatment and/or prognosis strategy.

### Goal 
In this project, we focused on search of transcriptomic biomarkers predicting functional activity of REST among its target genes, followed by 1) identification of clinical parameters predicting REST activity and 2) identification of potential signaling pathways that can be inhibited simultaneously with REST to achieve synergetic eradication of glioblastoma cells.

### Algorithm
The problem is that REST is controlled on post-translational level by upstream phosphatases and kinases, hence it is hard to predict its protein level simply via its mRNA expression. Here, we developed an algorithm to predict REST activity via expression of its target genes. The algorithm includes the following steps:
* Integration of whole-genome REST binding sites (ChIP-Seq) with transcriptomic changes in REST-knockout glioblastoma cells <https://pubmed.ncbi.nlm.nih.gov/38609948/>.
* Expression of top REST-target genes was used as a proxy for REST activity estimation using formula REST_activity = -mean(z-score, REST-target genes).
* Lasso regression model was fitted to identify clinical and genetic variables that influence REST activity.
* Next, patients samples were stratified into high-REST (REST_activity > median) and low-REST (REST_activity < median).
* Limma model on proteomic data was used to identify signaling pathways co-occurring with REST upregulation.


## Usage

This code was written and tested using R version 4.5.2 (2025-10-31) and publicly available tarnscriptomic (RNA-Seq) and proteomic (RPPA) data from TCGA-GBM project:
<https://www.cancer.gov/ccg/research/genome-sequencing/tcga> 

## Contributing

Pull requests (especially employing different types of predictive models with superior performance) are welcome.
