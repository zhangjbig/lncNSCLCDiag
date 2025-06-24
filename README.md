# lncNSCLCDiag
## R = 4.3.1 is used.
## Exosome.R

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

#### In exosome.R, the aim is to find lncRNAs from RNAseq data with significant differential expression, defined by statistical thresholds of p < 0.05 and |logFC| > 1. lncRNAs are considered that are clearly annotated as such in both the gene_type and transcript_type categories. 

## LASSO.R
### Environmental installation for LASSO.R
#### install.packages("ggvenn") == 0.1.10
#### install.packages("glmnet") == 4.1-8
#### install.packages("survival") == 3.5-8
#### install.packages("randomForestSRC") == 3.2.3

#### In LASSO.R, a LASSO regression model is constructed using differentially expressed lncRNAs (with both gene_type and transcript_type annotated as lncRNA) derived from exosome.R. The optimal regularization parameter is obtained via the lambda.min attribute of the model, and three key lncRNAs are selected.

## plasma.R
#### In plasma.R, differential expression analysis are performed on plasma data to identify lncRNAs with statistical significance (p < 0.05 and |logFC| > 1).Subsequently,  the three key lncRNAs selected in LASSO.R are examined to see whether they also differentially express in the plasma dataset.

## public validation.R
### Environmental installation
#### install.packages("future") == 1.33.2  
#### install.packages("future.apply") == 1.11.2 
#### install.packages("pbapply") == 1.7-2
#### install.packages("Parallel") == 4.3.1
#### install.packages("doParallel") == 1.0.17
#### install.packages("RColorBrewer") == 1.1-3

#### In public_validation.R, RNA-seq data of 12 different cancer types are collected from the TCGA dataset, along with tissue RNA-seq data from 8 independent cohorts of lung cancer in the GEO database.

## overlap.R
#### In overlap.R, differential expression analysis are conducted on RNA-seq data of lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) from the TCGA database (p < 0.05 and |logFC| > 1). LncRNAs exhibiting the same expression pattern (upregulation or downregulation) as those in exosomal data are filtered, leading to the identification of 9 eligible lncRNAs.

#### After completing the analysis in overlap.R and screening out 9 target lncRNA biomarkers for deep learning, we further extracted and saved the expression profiles of these 9 lncRNAs from the sequencing data processed in exosome.R, plasma.R, and public_validation.R respectively.

## python=3.9.13 is used.
## AERMLP.py
#### pip install pandas == 1.4.4
#### pip install numpy == 1.26.4
#### pip install tensorflow == 2.18.0
#### pip install scikit-learn==1.0.2
#### pip install matplotlib == 3.5.2
#### pip install tqdm == 4.64.1
#### pip install imbalanced-learn == 0.12.4

#### Input the RNA-seq data of 9 target lncRNAs from LUAD and LUSC cohorts in the TCGA database. First, address the class imbalance issue using SMOTE oversampling independently for each fold in 10-fold cross-validation. Specifically, apply SMOTE separately to the training set of each fold to generate synthetic samples, while keeping the validation set untouched to maintain real-world distribution. During model training, dynamically save the parameters of the best-performing fold and the corresponding autoencoder weights. After training, visualize the loss curve, accuracy curve, and AUC values for both the internal training and validation sets.

## MLP.py
#### pip install joblib == 1.4.2
#### pip install scipy == 1.13.1

#### Input the RNA sequencing data of 9 target lncRNAs of LUAD and LUSC from TCGA database, and train the model using ten-fold cross-validation method. During the training process, the model parameters with the best fold performance are saved. After the completion of training, the loss function change curve, the accuracy fluctuation curve, and the AUC values of the internal training and validation sets are visualized.

## lncNSCLCDiag.zip
## main.py
### Environmental installation
#### pip install pyqt5 == 5.15.11

#### main.py is used as the main program and four custom modules, PlaACC.py (plasma analysis), TissACC.py (tissue analysis), Batchpla.py (batch plasma processing), and BathTiss.py (batch tissue processing), are loaded via import statements to implement the complete data processing flow.

