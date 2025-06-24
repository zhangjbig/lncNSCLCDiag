library(plyr)
library(readr)
setwd("E:/github/data/exosome_gene_/")
path<-"E:/github/data/exosome_gene_/"
filenames <- dir(path)
n = length(filenames)

merge.data = read.table(file = filenames[1],header = T,sep = "\t")
merge.data1 <- merge.data[,c("gene_id","FPKM")]
colnames(merge.data1)[2] <- filenames[1]
merge.data1 <- merge.data1[!duplicated(merge.data1$gene_id),]

merge.data2 = read.table(file = filenames[2],header = T,sep = "\t")
merge.data2 <- merge.data2[,c("gene_id","FPKM")]
colnames(merge.data2)[2] <- filenames[2]
merge.data2<- merge.data2[!duplicated(merge.data2$gene_id),]
data <- merge(merge.data1,merge.data2,by='gene_id')

for(i in 3:n){
  shuju = read.table(file = filenames[i],header = T,sep = "\t")
  shuju = shuju[,c("gene_id","FPKM")]
  colnames(shuju)[2] <- filenames[i]
  shuju <- shuju[!duplicated(shuju$gene_id),]
  data <- merge(shuju,data,by="gene_id")
}

write.table(data,file = "E:/github/R/merge_all.csv")

BiocManager::install("rtracklayer",force = TRUE)
gtf <- rtracklayer::import("E:/github/R/gencode.v31lift37.annotation.gtf")
gtf_df <- as.data.frame(gtf)
geneid_df1 <- dplyr::select(gtf_df,c(gene_id,transcript_type,type,transcript_name,gene_type))

merge_data <- merge(data,geneid_df1,by = "gene_id")
data_lncRNA <- subset(merge_data,type=="transcript"&transcript_type=="lncRNA")
data_mRNA <- subset(merge_data,type!="transcript"&transcript_type=="protein_coding")

library(tidyr)
library(data.table)
BiocManager::install("org.Hs.eg.db",force = TRUE)
#install.packages("ggsignif")
#install.packages("RColorBrewer")
install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("beepr")
install.packages("gplots")
install.packages("pheatmap")
BiocManager::install("tidyverse")
library(org.Hs.eg.db)
library(limma)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler",force = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

exp <- data
head(exp)[1:6,1:6]
exp <- exp%>%separate(gene_id,c("ENSEMBL","meiyong"),"[.]")
exp <- exp[,-2]

#exp_lncRNA <- data_lncRNA
#exp_lncRNA  <- exp_lncRNA %>%separate(gene_id,c("ENSEMBL","meiyong"),"[.]")
#exp_lncRNA  <- exp_lncRNA [,-2]

exp_mRNA <- data_mRNA
exp_mRNA <- exp_mRNA%>%separate(gene_id,c("ENSEMBL","meiyong"),"[.]")
exp_mRNA <- exp_mRNA[,-2]

ID <- bitr(exp[,1],
           fromType = 'ENSEMBL',
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb = org.Hs.eg.db)
table(is.na(ID))

#ID_lncRNA <- bitr(exp_lncRNA[,1],
           #fromType = 'ENSEMBL',
           #toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           #OrgDb = org.Hs.eg.db)
#table(is.na(ID_lncRNA))

ID_mRNA <- bitr(exp_mRNA[,1],
                  fromType = 'ENSEMBL',
                  toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                  OrgDb = org.Hs.eg.db)
table(is.na(ID_mRNA))

exp_ID <- merge(exp,
                ID,
                by.x = "ENSEMBL",
                by.y = "ENSEMBL",all = F)

exp_ID <- exp_ID[,-1]
colnames(exp_ID)
exp_ID <- exp_ID[which(!duplicated(exp_ID$SYMBOL)),]
rownames(exp_ID) <- exp_ID$SYMBOL
head(sort(table(exp_ID$SYMBOL),decreasing = T),100)

#exp_ID1 <- merge(exp_lncRNA,
                #ID_lncRNA,
                #by.x = "ENSEMBL",
                #by.y = "ENSEMBL",all = F)

exp_ID1 <- data_lncRNA[,-1]
colnames(exp_ID1)
exp_ID1 <- exp_ID1[which(!duplicated(exp_ID1$transcript_name)),]
rownames(exp_ID1) <- exp_ID1$transcript_name

exp_ID2 <- merge(exp_mRNA,
                 ID_mRNA,
                 by.x = "ENSEMBL",
                 by.y = "ENSEMBL",all = F)

exp_ID2 <- exp_ID2[,-1]
colnames(exp_ID2)
exp_ID2 <- exp_ID2[which(!duplicated(exp_ID2$SYMBOL)),]
rownames(exp_ID2) <- exp_ID2$SYMBOL

exp_fpkm <- exp_ID[,c(1:19)]
exp_fpkm1 <- exp_ID1[,c(1:19)]
exp_fpkm2 <- exp_ID2[,c(1:19)]

library(stats)
# exp_fpkm_gene <- aggregate(exp_fpkm[1:20],
                          # by=list(exp_fpkm$SYMBOL),
                          # sum)
# exp_fpkm_gene <- column_to_rownames(exp_fpkm_gene,"Group.1")
exp_fpkm_gene <- exp_fpkm
exp_fpkm_gene1 <- exp_fpkm1
exp_fpkm_gene2 <- exp_fpkm2

# apply(exp_fpkm_gene,1,function(x){sum(floor(x)==0)>3})
exp_fpkm_gene <- exp_fpkm_gene[!apply(exp_fpkm_gene,1,function(x){sum(floor(x)==0)>3}),]
# 很多表达量为0的样本，直接选择某个基因如果在3个样本中的表达量为0，则直接舍去
exp_fpkm_gene1 <- exp_fpkm_gene1[!apply(exp_fpkm_gene1,1,function(x){sum(floor(x)==0)>19}),]
write.csv(exp_fpkm_gene1,"E:/github/R/exosome_lncRNA.csv")
exp_fpkm_gene2 <- exp_fpkm_gene2[!apply(exp_fpkm_gene2,1,function(x){sum(floor(x)==0)>19}),]
dim(exp_fpkm_gene)
dim(exp_fpkm_gene1)
dim(exp_fpkm_gene2)
boxplot(exp_fpkm_gene)
boxplot(exp_fpkm_gene1)
boxplot(exp_fpkm_gene2)

#limma归一化
exp_fpkm_gene <- normalizeBetweenArrays(exp_fpkm_gene)
boxplot(exp_fpkm_gene,las=2)
dat4 <- log2(exp_fpkm_gene + 1) 
boxplot(dat4,outline = F)

exp_fpkm_gene1 <- normalizeBetweenArrays(exp_fpkm_gene1)
boxplot(exp_fpkm_gene1,las=2)
dat4_1 <- log2(exp_fpkm_gene1 + 1) 
boxplot(dat4_1,outline = F)

exp_fpkm_gene2 <- normalizeBetweenArrays(exp_fpkm_gene2)
boxplot(exp_fpkm_gene2,las=2)
dat4_2 <- log2(exp_fpkm_gene2 + 1) 
boxplot(dat4_2,outline = F)

#表达矩阵数据校正
install.packages("pacman")
library(pacman)

# # 数据的FPKM 向TPM 的转换
expMatrix <- dat4
# # 构建一个将FPKM 转化为TPM 的函数，这是个固定公式，TPM是一种改进算法，适用于不同样本之间的比较
fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# 利用apply()函数对expMatrix的每一列进行FPKM 到TPM 的转化。
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,] #观察前三行
colSums(tpms) #加和每一列

expMatrix1 <- dat4_1
# 利用apply()函数对expMatrix的每一列进行FPKM 到TPM 的转化。
tpms1 <- apply(expMatrix1,2,fpkmToTpm)
tpms1[1:3,] #观察前三行
colSums(tpms1) #加和每一列

expMatrix2 <- dat4_2
# 利用apply()函数对expMatrix的每一列进行FPKM 到TPM 的转化。
tpms2 <- apply(expMatrix2,2,fpkmToTpm)
tpms2[1:3,] #观察前三行
colSums(tpms2) #加和每一列

#表达矩阵数据校正
exprSet <- tpms
boxplot(exprSet,outline=F,las=2)
exprSet1 <- tpms1
boxplot(exprSet1,outline=F,las=2)
exprSet2 <- tpms2
boxplot(exprSet2,outline=F,las=2)

# exprSet <- exprSet[!apply(exprSet,1,function(x){sum(floor(x)==0)>3}),]
dim(exprSet)
# boxplot(exprSet,las=2)

exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet, outline=F,las=2)
exprSet1=normalizeBetweenArrays(exprSet1)
boxplot(exprSet1, outline=F,las=2)
exprSet2=normalizeBetweenArrays(exprSet2)
boxplot(exprSet2, outline=F,las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
boxplot(exprSet,las=2,outline = F) 
exprSet1 <- log2(exprSet1+1)
boxplot(exprSet1,las=2,outline = F) 
exprSet2 <- log2(exprSet2+1)
boxplot(exprSet2,las=2,outline = F) 


# limma 的差异分析的过程
# 构建group_list
dat <- exprSet
dat_1 <- exprSet1
dat_2 <- exprSet2
group_list <- read.table("E:/github/R/gene_list",header = F,sep = "\t")

# colnames(group_list) <- group_list[order(group_list[,1]),]
new_group <- group_list[order(group_list[,1]),]
group <- group_list[,2]

#加载limma包，构造如下矩阵
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group)) #对样本实验设置进行计算。设计矩阵
colnames(design) <- levels(factor(group))
rownames(design) <- colnames(dat)

design1 <- model.matrix(~0+factor(group))
colnames(design1) <- levels(factor(group))
rownames(design1) <- colnames(dat_1)

design2 <- model.matrix(~0+factor(group))
colnames(design2) <- levels(factor(group))
rownames(design2) <- colnames(dat_2)

exp_fpkm_gene1.1 <- t(exp_fpkm_gene1)
os <- c(1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
os <- as.data.frame(os)
exp_fpkm_gene1.1  <- cbind(os,exp_fpkm_gene1.1)
boxplot(exp_fpkm_gene1.1$`ZNF667-AS1-205` ~ exp_fpkm_gene1.1$os,
        data = exp_fpkm_gene1.1,
        col = c("#8ECFC9","#FFBE7A"),
        xlab = "label",ylab = "Gene Expression",
        main = "STRA6LP")
axis(side = 1, at=c(1:2),labels = c("healthy","LUAD"))

exosome_BMS1P4 <- exp_fpkm_gene1.1$"AC022400.3-201"
exosome_BMS1P4 <- as.data.frame(exosome_BMS1P4)
rownames(exosome_BMS1P4) <- rownames(exp_fpkm_gene1.1)
exosome_AC108010.1 <- exp_fpkm_gene1.1$"AC108010.1-201"
exosome_AC108010.1 <- as.data.frame(exosome_AC108010.1)
rownames(exosome_AC108010.1) <- rownames(exp_fpkm_gene1.1)
exosome_STRA6LP <- exp_fpkm_gene1.1$"STRA6LP-201"
exosome_STRA6LP <- as.data.frame(exosome_STRA6LP)
rownames(exosome_STRA6LP) <- rownames(exp_fpkm_gene1.1)

install.packages("pROC")
library(pROC)
roc_data1 <- roc(exp_fpkm_gene1.1$os,as.numeric((unlist(exp_fpkm_gene1.1$"STRA6LP-201"))))
plot(roc_data1, main = "STRA6LP", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data1), 2)), col = "blue", cex = 1.2)
Youden <- roc_data1$sensitivities + roc_data1$specificities - 1
best_Youden <- which.max(Youden)
best_threshold <- roc_data1$thresholds[best_Youden]
print(best_threshold)

roc_data10 <- roc(exp_fpkm_gene1.1$os,as.numeric((unlist(exp_fpkm_gene1.1$"AC108010.1-201"))),levels=c("0", "1"))
roc_data11 <- roc(exp_fpkm_gene1.1$os,as.numeric((unlist(exp_fpkm_gene1.1$"AC022400.3-201"))),levels=c("0", "1"))
plot(roc_data10, add = T, col = "green", lwd = 2)
plot(roc_data11, add = T, col = "#80B1D3", lwd = 2)
legend("bottomright",legend = c("SUGT1P4-STRA6LP-CCDC180","AC108010.1","BMS1P4-AGAP5"),col = c("#80B1D3","#FDB462","#FB8072"),lwd = 2)

exp_fpkm_gene2.1 <- t(exp_fpkm_gene2)
exp_fpkm_gene2.1  <- cbind(os,exp_fpkm_gene2.1)
boxplot(exp_fpkm_gene2.1$"CD274"~ exp_fpkm_gene2.1$os,
        data = exp_fpkm_gene2.1,
        col = c("#8ECFC9","#FFBE7A"),
        xlab = "label",ylab = "Gene Expression",
        main = "SOX18")
roc_data2 <- roc(exp_fpkm_gene2.1$os,as.numeric((unlist(exp_fpkm_gene2.1$"CD274"))),levels=c("0", "1"))
plot(roc_data2, main = "KRT19", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data2), 2)), col = "blue", cex = 1.2)

exo.2 <- cbind(exp_fpkm_gene2.1$KRT19,exp_fpkm_gene2.1$CD274)
write.csv(exo.2,"E:/github/R/exo.2.csv")

model_1 <- glm(exp_fpkm_gene2.1[,1]~exp_fpkm_gene2.1$CD274+exp_fpkm_gene2.1$KRT19+exp_fpkm_gene1.1$`STRA6LP-201`, family = binomial(link ="logit"))
pred <- model_1$fitted.values
roc_multivar_1<-roc(exp_fpkm_gene2.1[,1],pred)
plot.roc(roc_multivar_1,col="red")   
text(0.8, 0.2, paste("AUC =", round(auc(roc_multivar_1), 2)), col = "blue", cex = 1.2)

model_2 <- glm(exp_fpkm_gene2.1[,1]~exp_fpkm_gene2.1$CD274+exp_fpkm_gene1.1$`STRA6LP-201`, family = binomial(link ="logit"))
pred <- model_2$fitted.values
roc_multivar_2<-roc(exp_fpkm_gene2.1[,1],pred)
plot.roc(roc_multivar_2,col="red")   
text(0.8, 0.2, paste("AUC =", round(auc(roc_multivar_2), 2)), col = "blue", cex = 1.2)

model_3 <- glm(exp_fpkm_gene2.1[,1]~exp_fpkm_gene2.1$KRT19+exp_fpkm_gene1.1$`STRA6LP-201`, family = binomial(link ="logit"))
pred <- model_3$fitted.values
roc_multivar_3<-roc(exp_fpkm_gene2.1[,1],pred)

roc_data1 <- roc(exp_fpkm_gene1.1$os,as.numeric((unlist(exp_fpkm_gene1.1$"STRA6LP-201"))))
plot(roc_data1, main = "STRA6LP", col = "#FB8072", lwd = 2, legacy.axes = TRUE)
# 重新绘制边框，只绘制下边框和左边框
box(which = "plot", bty = "L")
roc_data3 <- roc(exp_fpkm_gene2.1$os,as.numeric((unlist(exp_fpkm_gene2.1$CD274))),levels=c("0", "1"))
roc_data4 <- roc(exp_fpkm_gene2.1$os,as.numeric((unlist(exp_fpkm_gene2.1$KRT19))),levels=c("0", "1"))
plot(roc_data3, add = T, col = "#B3DE69", lwd = 2)
plot(roc_data4, add = T, col = "#BC80BD", lwd = 2)
plot(roc_multivar_2, add = T, col = "#FFA500", lwd = 2)
plot(roc_multivar_3, add = T, col = "#80B1D3", lwd = 2)
plot(roc_multivar_1, add = T, col = "red", lwd = 2)
legend("bottomright",legend = c("AUC=0.74 LNC-03","AUC=0.64 PD-L1","AUC=0.58 Cyfra21-1","AUC=0.77 LNC-03+PD-L1","AUC=0.79 LNC-03+Cyfra21-1","AUC=0.82 LNC-03+PD-L1+Cyfra21-1"),
       col = c("#FB8072","#B3DE69","#BC80BD","#FFA500","#80B1D3","red"),lwd = 2)

#比较design矩阵中两者的level
contrasts.matrix <- makeContrasts(NSCLC-healthy,
                                  levels = design)


# logFCcutoff <- log2(0)
# adjPvalueCutoff <- 0.05

#利用voom归一化
#library("edgeR")
#dge <- DGEList(counts = dat)
#dge <- calcNormFactors(dge)
#v <- voom(dge,design,plot=TRUE)

fit <- lmFit(dat,design) ##构建比对模型，比较两个条件下的表达数据
fit2 <- contrasts.fit(fit,contrasts.matrix) #contrasts.matrix对比矩阵，希望RNA样本之间进行那些比较
fit2 <- eBayes(fit2)
plotSA(fit2)
tempOutput <- topTable(fit2,coef = 1,n=Inf)
nrDEG = na.omit(tempOutput)

fit_1 <- lmFit(dat_1,design1) 
fit2_1 <- contrasts.fit(fit_1,contrasts.matrix) 
fit2_1 <- eBayes(fit2_1)
plotSA(fit2_1)
tempOutput1 <- topTable(fit2_1,coef = 1,n=Inf)
nrDEG1 = na.omit(tempOutput1)

fit_2 <- lmFit(dat_2,design2) 
fit2_2 <- contrasts.fit(fit_2,contrasts.matrix) 
fit2_2 <- eBayes(fit2_2)
plotSA(fit2_2)
tempOutput2 <- topTable(fit2_2,coef = 1,n=Inf)
nrDEG2 = na.omit(tempOutput2)


#查看差异基因数目
summary(decideTests(fit2,adjust.method = "none",p.value = 0.05,lfc = 0.58))#adjust.method = "none" p.value,而不是adj.p.value
summary(decideTests(fit2_1,adjust.method = "none",p.value = 0.05,lfc = 1))
summary(decideTests(fit2_2,adjust.method = "none",p.value = 0.05,lfc = 1))

Sdt <- decideTests(fit2,adjust.method = "none",p.value = 0.05,lfc = 0.58) #lfc=logfc
dt <- as.data.frame(dt)
de.common <- which(dt[,1]!=0)
length(de.common)

vennDiagram(dt[,1])
# vennDiagram(dt[,3:5])

colnames(fit2)

res_healthy_vs_NSCLC<- topTable(fit2,coef = 1,n=Inf) #求差异基因使用toptable函数 coef：列号或列名，指定线性模型的哪个系数或对比是感兴趣的 Inf表示无穷大
res_healthy_vs_NSCLC1<- topTable(fit2_1,coef = 1,n=Inf)
res_healthy_vs_NSCLC2<- topTable(fit2_2,coef = 1,n=Inf)

diffsig_healthy_vs_NSCLC <- res_healthy_vs_NSCLC[with(res_healthy_vs_NSCLC,(abs(logFC)>0.58 & P.Value<0.05)),]
diffsig_healthy_vs_NSCLC1 <- res_healthy_vs_NSCLC1[with(res_healthy_vs_NSCLC1,(abs(logFC)>1& P.Value<0.05)),]
diffsig_healthy_vs_NSCLC2 <- res_healthy_vs_NSCLC2[with(res_healthy_vs_NSCLC2,(abs(logFC)>1& P.Value<0.05)),]

diffsig_healthy_vs_NSCLC1.1 <- res_healthy_vs_NSCLC1[with(res_healthy_vs_NSCLC1,(abs(logFC)>0.58& P.Value<0.05)),]
write.csv(diffsig_healthy_vs_NSCLC1.1,"E:/github/R/0.58_exosome.lncRNA.csv")

BiocManager::install('EnhancedVolcano',force = TRUE)
BiocManager::install('patchwork',force = TRUE)
library(patchwork)
library(EnhancedVolcano)


p1 <- EnhancedVolcano(res_healthy_vs_NSCLC,
                      lab = rownames(res_healthy_vs_NSCLC),
                      x = 'logFC',
                      y = 'P.Value',
                      title = 'res_healthy_vs_NSCLC',
                      pointSize = 1.0,
                      legendPosition = 'right',
                      pCutoff = 0.05,
                      FCcutoff = 0.58
)

p2 <- EnhancedVolcano(res_healthy_vs_NSCLC1,
                      lab = NA,
                      x = 'logFC',
                      y = 'P.Value',
                      title = 'exosome_lncRNA',
                      pointSize = 1.0,
                      labSize = 3.0,
                      legendPosition = 'right',
                      pCutoff = 0.05,
                      FCcutoff = 1)

p3 <- EnhancedVolcano(res_healthy_vs_NSCLC2,
                      lab = rownames(res_healthy_vs_NSCLC2),
                      x = 'logFC',
                      y = 'P.Value',
                      title = 'exosome_mRNA',
                      pointSize = 1.0,
                      labSize = 3.0,
                      legendPosition = 'right',
                      pCutoff = 0.05,
                      FCcutoff = 1)

BiocManager::install('DT',force = TRUE)
library('DT')
library(pheatmap)
library(gplots)

diff1 <- cbind(row.names(diffsig_healthy_vs_NSCLC),diffsig_healthy_vs_NSCLC)
diff1_expr <- exp_fpkm_gene[rownames(diff1),]
diff <- diff1_expr
diff1_1 <- cbind(row.names(diffsig_healthy_vs_NSCLC1),diffsig_healthy_vs_NSCLC1)
diff1_expr1 <- exp_fpkm_gene1[rownames(diff1_1),]
diff_1 <- diff1_expr1
diff1_2 <- cbind(row.names(diffsig_healthy_vs_NSCLC2),diffsig_healthy_vs_NSCLC2)
diff1_expr2 <- exp_fpkm_gene2[rownames(diff1_2),]
diff_2 <- diff1_expr2

annotation_col <- group_list
label <- group_list[,-1]
label <- as.data.frame(label)
rownames(label) <- colnames(diff_2)
rownames(annotation_col) <- colnames(diff)

p_healthy <- pheatmap(diff,
                      annotation_col = annotation_col,
                      color = colorRampPalette(c("navy","white","red"))(50),
                      scale = "row",
                      treeheight_row = 200,
                      treeheight_col = 30,
                      border_color = NA,
                      fontsize = 10,
                      show_rownames = T
)
#write.csv(diff,"E:/github/R/diff_healthy.csv")

p_healthy1 <- pheatmap(diff_1,
                      annotation_col = label,
                      color = colorRampPalette(c("navy","white","red"))(50),
                      scale = "row",
                      treeheight_row = 100,
                      treeheight_col = 30,
                      border_color = NA,
                      fontsize = 9,
                      show_rownames = T
)
write.csv(diff_1,"E:/github/R/exosome_lncRNA1.csv")

p_healthy2 <- pheatmap(diff_2,
                       annotation_col = label,
                       color = colorRampPalette(c("navy","white","red"))(50),
                       scale = "row",
                       treeheight_row = 100,
                       treeheight_col = 30,
                       border_color = NA,
                       fontsize = 9,
                       show_rownames = T
)

library(dplyr)
diff_1 <- as.data.frame(diff_1)
lnc1 <- cbind(row.names(diff_1),diff_1)
colnames(lnc1)[1] <- "transcript_name"
data_lncRNA1 <- data_lncRNA[,21:24]
lnc <- merge(lnc1,data_lncRNA1,by = "transcript_name")
lnc <- subset(lnc,gene_type=="lncRNA")
lnc_1 <- lnc
rownames(lnc_1) <- lnc$transcript_name
lnc_1 <- lnc_1[,c(-1,-21:-23)]

write.csv(lnc_1,"E:/github/R/exosome.lncRNA.csv")

exosome <- read_csv("E:/github/R/exosome.csv")
exosome <- exosome[,-1]
rownames(exosome) <- rownames(label)
exosome <- t(exosome)
p_exosome <- pheatmap(exosome,
                      annotation_col = label,
                      color = colorRampPalette(c("navy","white","red"))(50),
                      scale = "row",
                      treeheight_row = 200,
                      treeheight_col = 30,
                      border_color = NA,
                      fontsize = 10,
                      show_rownames = T
)
# 假设label是一个数据框，包含样本名和对应的标签列
# 这里假设标签列名为"label_column_name"
ordered_indices <- order(label[,1])
exosome_ordered <- exosome[, ordered_indices]
label_ordered <- label[ordered_indices, ]
label_ordered <- as.data.frame(label_ordered)
exosome_ordered <- t(exosome_ordered)
rownames(label_ordered) <- rownames(exosome_ordered)
exosome_ordered <- t(exosome_ordered)
p_exosome <- pheatmap(exosome_ordered,
                      annotation_col = label_ordered,
                      color = colorRampPalette(c("navy","white","red"))(50),
                      scale = "row",
                      treeheight_row = 200,
                      treeheight_col = 30,
                      border_color = NA,
                      fontsize = 10,
                      show_rownames = T,
                      cluster_cols = FALSE
)

tcga <- read_csv("E:/github/R/tcga.csv")
# 生成行名
row_names <- paste0("sample", 1:nrow(tcga))
# 为数据框添加行名
rownames(tcga) <- row_names
tcga.label <- tcga[,1]
rownames(tcga.label) <- rownames(tcga)
tcga <- t(tcga)
tcga <- tcga[-1,]
tcga.label.1 <- tcga.label
tcga.label.1[,1] <- ifelse(tcga.label[,1] == "0","Healtyh","NSCLC")
rownames(tcga.label.1) <- rownames(tcga.label)
tcga.label.1 <- as.data.frame(tcga.label.1)
rownames(tcga.label.1) <- rownames(tcga.label)
p_tcga <- pheatmap(tcga,
                   annotation_col = tcga.label.1,
                   color = colorRampPalette(c("navy","white","red"))(6),
                   scale = "row",
                   treeheight_row = 200,
                   treeheight_col = 30,
                   border_color = NA,
                   fontsize = 10,
                   show_rownames = T)
# 假设label是一个数据框，包含样本名和对应的标签列
# 这里假设标签列名为"label_column_name"
# 根据标签排序
ordered_indices <- order(tcga.label.1[,1])
tcga_ordered <- tcga[, ordered_indices]
tcga.label.1_ordered <- tcga.label.1[ordered_indices, , drop = FALSE]
# 绘制热图，同时禁止列聚类，保留排序顺序
p_tcga <- pheatmap(tcga_ordered,
                   annotation_col = tcga.label.1_ordered,
                   color = colorRampPalette(c("navy","white","red"))(6),
                   scale = "row",
                   treeheight_row = 200,
                   treeheight_col = 30,
                   border_color = NA,
                   fontsize = 10,
                   show_rownames = TRUE,
                   cluster_cols = FALSE)
# 自定义基因顺序
custom_order <- c("OVCH1-AS1", "RTCA-AS1", "STRA6LP", "MMP25-AS1", "TRG-AS1", "AL445524.1", "OR2A1-AS1", "SNHG1", "TTC28-AS1")
tcga_ordered.1 <- tcga_ordered[match(custom_order, rownames(tcga_ordered)), ]
p_tcga.1 <- pheatmap(tcga_ordered.1,
                     annotation_col = tcga.label.1_ordered,
                     color = colorRampPalette(c("navy","white","red"))(6),
                     scale = "row",
                     treeheight_row = 20,
                     treeheight_col = 30,
                     border_color = NA,
                     fontsize = 10,
                     show_rownames = TRUE,
                     cluster_cols = FALSE)
# 自定义基因顺序
custom_order <- c("TRG-AS1","OVCH1-AS1", "RTCA-AS1", "STRA6LP", "MMP25-AS1", "AL445524.1", "OR2A1-AS1", "SNHG1", "TTC28-AS1")
# 根据自定义顺序重新排列tcga数据框的行
tcga <- tcga[match(custom_order, rownames(tcga)), ]

# 根据标签排序
ordered_indices <- order(tcga.label.1[,1])
tcga_ordered <- tcga[, ordered_indices]
tcga.label.1_ordered <- tcga.label.1[ordered_indices, , drop = FALSE]
breaks_seq <- seq(-3, 3, length.out = 11) 
# 绘制热图，同时禁止列聚类，保留排序顺序
p_tcga <- pheatmap(tcga_ordered,
                   annotation_col = tcga.label.1_ordered,
                   color = colorRampPalette(c("navy","white","red"))(10),
                   scale = "row",
                   treeheight_row = 200,
                   treeheight_col = 30,
                   border_color = NA,
                   fontsize = 10,
                   show_rownames = TRUE,
                   cluster_cols = FALSE,
                   cluster_rows = F,
                   breaks = breaks_seq)

# 按行将数据归一化到 -1 到 1 的范围
norm_tcga_ordered <- t(apply(tcga_ordered, 1, function(x) {
  min_val <- min(x)
  max_val <- max(x)
  if (min_val == max_val) {
    return(rep(0, length(x)))  # 处理最小值和最大值相等的情况
  }
  return(2 * ((x - min_val) / (max_val - min_val)) - 1)
}))
p_tcga <- pheatmap(norm_tcga_ordered,
                   annotation_col = tcga.label.1_ordered,
                   color = colorRampPalette(c("navy","white","red"))(20),
                   scale = "none",  # 这里因为已经归一化，可不进行额外缩放
                   treeheight_row = 200,
                   treeheight_col = 30,
                   border_color = NA,
                   fontsize = 10,
                   show_rownames = TRUE,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE)
