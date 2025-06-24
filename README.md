# lncNSCLCDiag
### Exosome.R
#### In exosome.R, the aim is to find lncRNAs from RNAseq data with significant differential expression, defined by statistical thresholds of p < 0.05 and |logFC| > 1. lncRNAs are considered that are clearly annotated as such in both the gene_type and transcript_type categories. 
### Environmental installation for Exosome.R
#### install.packages("plyr") == 1.8.9
#### install.packages("readr") == 2.1.5
#### BiocManager::install(“rtracklayer”,force = T) == 1.60.1
#### install.packages("tidyr") == 1.3.1
#### install.packages("data.table") == 1.15.2
#### BiocManager::install(“org.Hs.eg.db”,force = T) == 3.17.0
#### install.packages("data.table") == 1.15.2
#### install.packages("ggplot2") == 3.5.0
#### install.packages("gplots") == 3.1.3.1
#### install.packages("pheatmap") == 1.0.12
#### BiocManager::install(“tidyverse”) == 2.0.0
#### install.packages("limma") == 3.56.2  
#### BiocManager::install(“clusterProfiler”) == 4.8.3
#### install.packages("BiocManager") == 1.30.22
#### install.packages("stats") == 4.3.1
#### install.packages("pacman") == 0.5.1
#### install.packages("pROC") == 1.18.5
#### BiocManager::install(“EnhancedVolcano”) == 1.18.0
#### BiocManager::install(“patchwork”) == 1.2.0
#### BiocManager::install(“DT ”) == 0.32
#### install.packages("pROC") == 1.18.5
#### install.packages("dplyr") == 1.1.4
#### install.packages("tibble") == 3.2.1

### LASSO.R
### Environmental installation for LASSO.R
#### install.packages("ggvenn") == 0.1.10
#### install.packages("glmnet") == 4.1-8
#### install.packages("survival") == 3.5-8
#### install.packages("randomForestSRC") == 3.2.3

#### In LASSO.R, a LASSO regression model is constructed using differentially expressed lncRNAs (with both gene_type and transcript_type annotated as lncRNA) derived from exosome.R. The optimal regularization parameter is obtained via the lambda.min attribute of the model, and three key lncRNAs are selected.
