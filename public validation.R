setwd("E:/github/")
path<-"E:/github/"
install.packages("pROC")
library(pROC)

data4 <- read.table("E:/github/data/GSE40275/GSE40275_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T)
data4 <- data4[c(-1:-36),]
data4.1 <- read.table("E:/github/data/GSE40275/GPL15974_family.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
BiocManager::install("org.Hs.eg.db",force = TRUE)
library(org.Hs.eg.db)
BiocManager::install("clusterProfiler",force = TRUE)
library(clusterProfiler)
data4.2 <- bitr(data4.1[,2],
           fromType = 'ENTREZID',
           toType='SYMBOL',
           OrgDb = org.Hs.eg.db)
data4.1 <- merge(data4.1,data4.2,by.x="ENTREZ_GENE_ID",by.y="ENTREZID")
data.GSE40275 <- merge(data4.1,data4,by.x="ID",by.y="X.Sample_title")
data.GSE40275 <- data.GSE40275[which(!duplicated(data.GSE40275$SYMBOL)),]
rownames(data.GSE40275) <- data.GSE40275[,3]
data.GSE40275 <- data.GSE40275[,c(-1:-3)]
data.GSE40275 <- t(data.GSE40275)
os4 <- c(rep(0, 37),rep(1,16))
data.GSE40275 <- cbind(os4,data.GSE40275)
data.GSE40275 <- as.data.frame(data.GSE40275)

roc_data4 <- roc(data.GSE40275$os4,as.numeric((unlist(data.GSE40275$"SUGT1P4-STRA6LP-CCDC180"))),levels=c("0", "1"))
plot(roc_data4, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data4), 2)), col = "blue", cex = 1.2)
data.GSE40275.1 <- data.GSE40275
data.GSE40275.1[,1] <- ifelse(data.GSE40275[,1] == "0","healthy","NSCLC")
bp4 <- boxplot(as.numeric((unlist(data.GSE40275.1$"SUGT1P4-STRA6LP-CCDC180"))) ~ data.GSE40275.1[,1],
               data = data.GSE40275.1,
               col = c("#8ECFC9","#FFBE7A"),
               xlab = "label",ylab = "Gene Expression",
               main = "SUGT1P4-STRA6LP-CCDC180")
bp4$stats[3,]
means4 <- tapply(as.numeric((unlist(data.GSE40275$"SUGT1P4-STRA6LP-CCDC180"))),data.GSE40275[,1],mean)
means4
Data8.GSE40275 <- cbind(data.GSE40275[,1],data.GSE40275$SNHG1,data.GSE40275$"SUGT1P4-STRA6LP-CCDC180")
write.csv(Data8.GSE40275,"E:/github/R/Data8.GSE40275.csv")

data7 <- read.table("E:/github/data/GSE19188/GSE19188_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T)
data7 <- t(data7)
data7 <- as.data.frame(data7)
colnames(data7) <- data7[1,]
#data7 <- subset(data7,data7[,2]=="cell type: ADC" | data7[,2]=="cell type: healthy")
data7.1 <- read.table("E:/github/data/GSE19188/GPL570-55999.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data7 <- data7[,-1]
data7 <- t(data7)
data7.2 <- rownames(data7)
data7.2 <- as.data.frame(data7.2)
data7.2 <- cbind(data7.2,data7)
data.GSE19188 <- merge(data7.1,data7.2,by.x="ID",by.y="data7.2")
data.GSE19188 <- data.GSE19188[which(!duplicated(data.GSE19188$Gene.Symbol)),]
rownames(data.GSE19188) <- data.GSE19188[,2]
data.GSE19188 <- data.GSE19188[,c(-1:-3)]
data7.3 <- t(data7)
data7.3 <- data7.3[-1,]
data7.3 <- data7.3[,1]
data7.3 <- as.data.frame(data7.3)
data.GSE19188 <- t(data.GSE19188)
data.GSE19188 <- cbind(data7.3,data.GSE19188)
data.GSE19188 <- as.data.frame(data.GSE19188)
means7 <- tapply(as.numeric((unlist(data.GSE19188$"CCDC180 /// LOC100499484-C9ORF174"))),data.GSE19188[,1],mean)
means7
bp7 <- boxplot(as.numeric((unlist(data.GSE19188$"CCDC180 /// LOC100499484-C9ORF174"))) ~ data.GSE19188[,1],
               data = data.GSE19188,
               col = c("#FFBE7A","#8ECFC9","#FA7F6F","#82B0D2"),
               xlab = "label",ylab = "Gene Expression",
               main = "SUGT1P4-STRA6LP-CCDC180")
bp7$stats[3,]
data.GSE19188.1 <- data.GSE19188
data.GSE19188.1[,1] <- ifelse(data.GSE19188.1[,1] == "cell type: healthy",0,1)
roc_data7 <- roc(data.GSE19188.1[,1],as.numeric((unlist(data.GSE19188.1$"CCDC180 /// LOC100499484-C9ORF174"))),levels=c("0", "1"))
plot(roc_data7, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data7), 2)), col = "blue", cex = 1.2)
Data4.GSE19188 <- cbind(data.GSE19188.1[,1],data.GSE19188$"LOC101928605 /// OR2A1-AS1 /// OR2A20P /// OR2A9P",data.GSE19188$`OVCH1-AS1`,data.GSE19188$`SNHG1 /// SNORD22 /// SNORD25 /// SNORD26 /// SNORD27 /// SNORD28 /// SNORD29 /// SNORD30 /// SNORD31`,
                        data.GSE19188$"CCDC180 /// LOC100499484-C9ORF174",data.GSE19188$`TRG-AS1`,data.GSE19188$`TTC28-AS1`)
write.csv(Data4.GSE19188,"E:/github/R/Data4.GSE19188.csv")

data9 <- read.table("E:/github/data/GSE19804/GSE19804_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T)
data9 <- data9[c(-1:-8,-10:-34),]
data9.1 <- read.table("E:/github/data/GSE19188/GPL570-55999.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data.GSE19804 <- merge(data9.1,data9,by.x="ID",by.y="X.Sample_title")
data.GSE19804 <- data.GSE19804[which(!duplicated(data.GSE19804$Gene.Symbol)),]
rownames(data.GSE19804) <- data.GSE19804[,2]
data.GSE19804 <- data.GSE19804[,c(-1,-2)]
data9.2 <- t(data9)
data9.2 <- data9.2[-1,]
data.GSE19804 <- t(data.GSE19804)
data.GSE19804 <- cbind(data9.2[,1],data.GSE19804)
data.GSE19804 <- as.data.frame(data.GSE19804)

data.GSE19804[,1] <- ifelse(data.GSE19804[,1] == "tissue: paired normal adjacent",0,1)
roc_data9 <- roc(data.GSE19804[,1],as.numeric((unlist(data.GSE19804$"CCDC180 /// LOC100499484-C9ORF174"))),levels=c("0", "1"))
plot(roc_data9, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data9), 2)), col = "blue", cex = 1.2)
data.GSE19804.1 <- data.GSE19804
data.GSE19804.1[,1] <- ifelse(data.GSE19804[,1] == "0","healthy","lung cancer")
bp9 <- boxplot(as.numeric((unlist(data.GSE19804.1$"CCDC180 /// LOC100499484-C9ORF174"))) ~ data.GSE19804.1[,1],
               data = data.GSE19804.1,
               col = c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2"),
               xlab = "label",ylab = "Gene Expression",
               main = "SUGT1P4-STRA6LP-CCDC180")
bp9$stats[3,]
means9 <- tapply(as.numeric((unlist(data.GSE19804$"CCDC180 /// LOC100499484-C9ORF174"))),data.GSE19804[,1],mean)
means9

Data7.GSE19804 <- cbind(data.GSE19804[,1],data.GSE19804$"LOC101928605 /// OR2A1-AS1 /// OR2A20P /// OR2A9P",data.GSE19804$`OVCH1-AS1`,data.GSE19804$`SNHG1 /// SNORD22 /// SNORD25 /// SNORD26 /// SNORD27 /// SNORD28 /// SNORD29 /// SNORD30 /// SNORD31`,
                        data.GSE19804$"CCDC180 /// LOC100499484-C9ORF174",data.GSE19804$`TRG-AS1`,data.GSE19804$`TTC28-AS1`)
write.csv(Data7.GSE19804,"E:/github/R/Data1.GSE19804.csv")

data10 <- read.table("E:/github/data/GSE18842/GSE18842_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T)
data10 <- data10[c(-1:-9,-11:-32),]
data10.1 <- read.table("E:/github/data/GSE19188/GPL570-55999.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data.GSE18842 <- merge(data10.1,data10,by.x="ID",by.y="X.Sample_title")
data.GSE18842 <- data.GSE18842[which(!duplicated(data.GSE18842$Gene.Symbol)),]
rownames(data.GSE18842) <- data.GSE18842[,2]
data.GSE18842 <- data.GSE18842[,c(-1,-2)]
data10.2 <- t(data10)
data10.2 <- data10.2[-1,]
data.GSE18842 <- t(data.GSE18842)
data.GSE18842 <- cbind(data10.2[,1],data.GSE18842)
data.GSE18842 <- as.data.frame(data.GSE18842)

roc_data10 <- roc(data.GSE18842[,1],as.numeric((unlist(data.GSE18842$"CCDC180 /// LOC100499484-C9ORF174"))),levels=c("sample type: control", "sample type: tumor"))
plot(roc_data10, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data10), 2)), col = "blue", cex = 1.2)
bp10 <- boxplot(as.numeric((unlist(data.GSE18842$"CCDC180 /// LOC100499484-C9ORF174"))) ~ data.GSE18842[,1],
                data = data.GSE18842,
                col = c("#8ECFC9","#FFBE7A"),
                xlab = "label",ylab = "Gene Expression",
                main = "SUGT1P4-STRA6LP-CCDC180")
bp10$stats[3,]
means10 <- tapply(as.numeric((unlist(data.GSE18842$"CCDC180 /// LOC100499484-C9ORF174"))),data.GSE18842[,1],mean)
means10
data.GSE18842.1[,1] <- ifelse(data.GSE18842[,1] == "sample type: control",0,1)

Data5.GSE18842 <- cbind(data.GSE18842.1[,1],data.GSE18842$"LOC101928605 /// OR2A1-AS1 /// OR2A20P /// OR2A9P",data.GSE18842$`OVCH1-AS1`,data.GSE18842$`SNHG1 /// SNORD22 /// SNORD25 /// SNORD26 /// SNORD27 /// SNORD28 /// SNORD29 /// SNORD30 /// SNORD31`,
                        data.GSE18842$"CCDC180 /// LOC100499484-C9ORF174",data.GSE18842$`TRG-AS1`,data.GSE18842$`TTC28-AS1`)
write.csv(Data5.GSE18842,"E:/github/R/Data5.GSE18842.csv")

data12 <- read.table("E:/github/data/GSE33532/GSE33532_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data12 <- data12[c(-1:-11,-13:-46),]
data12.1 <- read.table("E:/github/data/GSE19188/GPL570-55999.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data.GSE33532 <- merge(data12.1,data12,by.x="ID",by.y="X.Sample_title")
data.GSE33532<- data.GSE33532[which(!duplicated(data.GSE33532$Gene.Symbol)),]
rownames(data.GSE33532) <- data.GSE33532[,2]
data.GSE33532 <- data.GSE33532[,c(-1,-2)]
data12.2 <- t(data12)
data12.2 <- data12.2[-1,]
data.GSE33532 <- t(data.GSE33532)
data.GSE33532 <- cbind(data12.2[,1],data.GSE33532)
data.GSE33532 <- as.data.frame(data.GSE33532)

data.GSE33532.1 <- data.GSE33532
data.GSE33532.1[,1] <- ifelse(data.GSE33532.1[,1] == "histology: normal lung",0,1)
roc_data12 <- roc(data.GSE33532.1[,1],as.numeric((unlist(data.GSE33532.1$`CCDC180 /// LOC100499484-C9ORF174`))),levels=c("0", "1"))
plot(roc_data12, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data12), 2)), col = "blue", cex = 1.2)
bp12 <- boxplot(as.numeric((unlist(data.GSE33532$`CCDC180 /// LOC100499484-C9ORF174`))) ~ data.GSE33532[,1],
                data = data.GSE33532,
                col = c("#FFBE7A","#FA7F6F","#8ECFC9","#82B0D2"),
                xlab = "label",ylab = "Gene Expression",
                main = "SUGT1P4-STRA6LP-CCDC180")
bp12
means12 <- tapply(as.numeric((unlist(data.GSE33532$PCSK7))),data.GSE33532[,1],mean)
means12

Data6.GSE33532 <- cbind(data.GSE33532.1[,1],data.GSE33532$`RP11-295G20.2`,data.GSE33532$"LOC101928605 /// OR2A1-AS1 /// OR2A20P /// OR2A9P",
                        data.GSE33532$"CCDC180 /// LOC100499484-C9ORF174",data.GSE33532$`TTC28-AS1`)
write.csv(Data6.GSE33532,"E:/github/R/Data6.GSE33532.csv")

data16 <- read.table("E:/github/data/GSE27262/GSE27262_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data16 <- data16[c(-1:-8,-10:-37),]
data16.1 <- read.table("E:/github/GSE19188/data/GPL570-55999.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data.GSE27262 <- merge(data16.1,data16,by.x="ID",by.y="X.Sample_title")
data.GSE27262 <- data.GSE27262[which(!duplicated(data.GSE27262[,2])),]
rownames(data.GSE27262) <- data.GSE27262[,2]
data.GSE27262 <- data.GSE27262[,c(-1,-2)]
data16.2 <- t(data16)
data16.2 <- data16.2[-1,]
data.GSE27262 <- t(data.GSE27262)
data.GSE27262 <- cbind(data16.2[,1],data.GSE27262)
data.GSE27262 <- as.data.frame(data.GSE27262)
#data.GSE27262 <- subset(data.GSE27262,data.GSE27262[,1]=="egfr mutation: NA" | data.GSE27262[,1]=="egfr mutation: Yes")
data.GSE27262.1 <- data.GSE27262
means16 <- tapply(as.numeric((unlist(data.GSE27262$"CCDC180 /// LOC100499484-C9ORF174"))),data.GSE27262[,1],mean)
means16
bp16 <- boxplot(as.numeric((unlist(data.GSE27262$"CCDC180 /// LOC100499484-C9ORF174"))) ~ data.GSE27262[,1],
                data = data.GSE27262,
                col = c("#8ECFC9","#FFBE7A","#FA7F6F"),
                xlab = "label",ylab = "Gene Expression",
                main = "CCDC180 /// LOC100499484-C9ORF174")
bp16$stats[3,]
data.GSE27262.1[,1] <- ifelse(data.GSE27262.1[,1] == "egfr mutation: NA",0,1)
roc_data16 <- roc(data.GSE27262.1[,1],as.numeric((unlist(data.GSE27262.1$"CCDC180 /// LOC100499484-C9ORF174"))),levels=c("0", "1"))
plot(roc_data16, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data16), 2)), col = "blue", cex = 1.2)
Data2.GSE27262 <- cbind(data.GSE27262.1[,1],data.GSE27262$"LOC101928605 /// OR2A1-AS1 /// OR2A20P /// OR2A9P",data.GSE27262$`OVCH1-AS1`,data.GSE27262$`SNHG1 /// SNORD22 /// SNORD25 /// SNORD26 /// SNORD27 /// SNORD28 /// SNORD29 /// SNORD30 /// SNORD31`,
                        data.GSE27262$"CCDC180 /// LOC100499484-C9ORF174",data.GSE27262$`TRG-AS1`,data.GSE27262$`TTC28-AS1`)
write.csv(Data2.GSE27262,"E:/github/R/Data2.GSE27262.csv")

data26 <- read.table("E:/github/data/GSE102287/GSE102287-GPL570_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data26 <- data26[c(-1:-14,-16:-59),]
data26.1 <- read.table("E:/github/GSE19188/data/GPL570-55999.csv", comment.char = "#",header = T,sep = ",",fill = T)
data.GSE102287 <- merge(data26.1,data26,by.x="ID",by.y="X.Sample_title")
data.GSE102287 <- data.GSE102287[which(!duplicated(data.GSE102287$Gene.Symbol)),]
rownames(data.GSE102287) <- data.GSE102287[,2]
data.GSE102287 <- data.GSE102287[,c(-1:-2)]
data26.2 <- t(data26)
data26.2 <- data26.2[-1,]
data.GSE102287 <- t(data.GSE102287)
data.GSE102287 <- cbind(data26.2[,1],data.GSE102287)
data.GSE102287 <- as.data.frame(data.GSE102287)
data.GSE102287[,1] <- ifelse(data.GSE102287[,1] == "tumor_normal status: N",0,1)
roc_data32 <- roc(data.GSE102287[,1],as.numeric((unlist(data.GSE102287$"CCDC180 /// LOC100499484-C9ORF174"))),levels=c("0", "1"))
plot(roc_data32, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data32), 2)), col = "blue", cex = 1.2)
data.GSE102287.1 <- data.GSE102287
data.GSE102287.1[,1] <- ifelse(data.GSE102287[,1] == "0","healthy","NSCLC")
bp32 <- boxplot(as.numeric((unlist(data.GSE102287.1$"CCDC180 /// LOC100499484-C9ORF174"))) ~ data.GSE102287.1[,1],
                data = data.GSE102287.1,
                col = c("#8ECFC9","#FFBE7A"),
                xlab = "label",ylab = "Gene Expression",
                main = "SUGT1P4-STRA6LP-CCDC180")
bp32$stats[3,]
means32 <- tapply(as.numeric((unlist(data.GSE102287$"CCDC180 /// LOC100499484-C9ORF174"))),data.GSE102287[,1],mean)
means32

Data3.GSE102287 <- cbind(data.GSE102287[,1],data.GSE102287$"LOC101928605 /// OR2A1-AS1 /// OR2A20P /// OR2A9P",data.GSE102287$`OVCH1-AS1`,data.GSE102287$`SNHG1 /// SNORD22 /// SNORD25 /// SNORD26 /// SNORD27 /// SNORD28 /// SNORD29 /// SNORD30 /// SNORD31`,
                         data.GSE102287$"CCDC180 /// LOC100499484-C9ORF174",data.GSE102287$`TRG-AS1`,data.GSE102287$`TTC28-AS1`)
write.csv(Data3.GSE102287,"E:/github/R/Data3.GSE102287.csv")

data21 <- read.table("E:/github/data/GSE176348/GSE176348_series_matrix.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data21 <- data21[c(-1:-21,-23:-40),]
data21.1 <- read.table("E:/github/data/GSE176348/GPL25759_family.csv", comment.char = "#",header = T,sep = ",",fill = T,quote = "")
data.GSE176348 <- merge(data21.1,data21,by.x="ID",by.y="X.Sample_title")
data.GSE176348 <- data.GSE176348[which(!duplicated(data.GSE176348$GENE.SYMBOL)),]
rownames(data.GSE176348) <- data.GSE176348[,4]
data.GSE176348 <- data.GSE176348[,c(-1:-14)]
data21.2 <- t(data21)
data21.2 <- data21.2[-1,]
data.GSE176348 <- t(data.GSE176348)
data.GSE176348 <- cbind(data21.2[,1],data.GSE176348)
data.GSE176348 <- as.data.frame(data.GSE176348)

data.GSE176348[,1] <- ifelse(data.GSE176348[,1] == "adjacent non-tumor tissues",0,1)
roc_data27 <- roc(data.GSE176348[,1],as.numeric((unlist(data.GSE176348$"RP11-23J9.4"))),levels=c("0", "1"))
plot(roc_data21, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data27), 2)), col = "blue", cex = 1.2)
data.GSE176348.1 <- data.GSE176348
data.GSE176348.1[,1] <- ifelse(data.GSE176348[,1] == "0" ,"healthy","LUAD")
bp27 <- boxplot(as.numeric((unlist(data.GSE176348.1$"RP11-23J9.4"))) ~ data.GSE176348.1[,1],
                data = data.GSE176348.1,
                col = c("#8ECFC9","#FFBE7A"),
                xlab = "label",ylab = "Gene Expression",
                main = "SUGT1P4-STRA6LP-CCDC180")
bp27$stats[3,]
means27 <- tapply(as.numeric((unlist(data.GSE176348$CDKN3))),data.GSE176348[,1],mean)
means27
Data1.GSE176348 <- cbind(data.GSE176348[,1],data.GSE176348$SNHG1,data.GSE176348$"RP11-23J9.4",data.GSE176348$`TRG-AS1`)
write.csv(Data1.GSE176348,"E:/github/R/Data7.GSE176348.csv")

#TCGA
setwd("E:/github/data/TCGA")
path<-"E:/github/data/TCGA"
install.packages("future")
install.packages("future.apply")
install.packages("pbapply")
install.packages("parallel")
install.packages("doParallel")
library(limma)
library(future)
library(future.apply)
library(pbapply)
library(parallel)
library(doParallel)
#加载所需函数
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM=merge_TCGA(metadata=metaMatrix.RNA, 
                   path="RNAseq", 
                   data.type="RNAseq",
                   mRNA_expr_type="TPM",
                   symbol = T,
                   RNA_type=T
)
#提取所有表达数据
RNA_TPM_all=as.matrix(RNA_TPM)
rownames(RNA_TPM_all)=RNA_TPM_all[,1]
RNA_TPM_all_exp=RNA_TPM_all[,3:ncol(RNA_TPM_all)]
RNA_TPM_all_dimnames=list(rownames(RNA_TPM_all_exp),colnames(RNA_TPM_all_exp))
RNA_TPM_all_data=matrix(as.numeric(as.matrix(RNA_TPM_all_exp)),nrow=nrow(RNA_TPM_all_exp),dimnames=RNA_TPM_all_dimnames)
RNA_TPM_all_data=avereps(RNA_TPM_all_data)
write.table(file="TCGA_all_TPM.txt",RNA_TPM_all_data,sep="\t",quote=F)

#提取lncRNA
RNA_TPM_lnc=RNA_TPM[RNA_TPM$type=="lncRNA",]
RNA_TPM_lnc=as.matrix(RNA_TPM_lnc)
rownames(RNA_TPM_lnc)=RNA_TPM_lnc[,1]
RNA_TPM_lnc_exp=RNA_TPM_lnc[,3:ncol(RNA_TPM_lnc)]
RNA_TPM_lnc_dimnames=list(rownames(RNA_TPM_lnc_exp),colnames(RNA_TPM_lnc_exp))
RNA_TPM_lnc_data=matrix(as.numeric(as.matrix(RNA_TPM_lnc_exp)),nrow=nrow(RNA_TPM_lnc_exp),dimnames=RNA_TPM_lnc_dimnames)
RNA_TPM_lnc_data=avereps(RNA_TPM_lnc_data)
write.table(file="TCGA_lnc_TPM.txt",RNA_TPM_lnc_data,sep="\t",quote=F)

RNA_TPM_all_data_t <- t(RNA_TPM_all_data)
metaMatrix.RNA1 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1 <- cbind(rownames(RNA_TPM_all_data_t),RNA_TPM_all_data_t)
RNA_TPM_all_data1 <- merge(metaMatrix.RNA1,RNA_TPM_all_data1,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1) <- RNA_TPM_all_data1[,1]
RNA_TPM_all_data1 <- RNA_TPM_all_data1[,-1]
RNA_TPM_all_data1[,1] <- ifelse(RNA_TPM_all_data1[,1] == "SolidTissueNormal",0,1)
#roc_data21 <- roc(RNA_TPM_all_data1[,1],as.numeric((unlist(RNA_TPM_all_data1$"STRA6LP"))),levels=c("0", "1"))
#plot(roc_data21, main = "SUGT1P4-STRA6LP-CCDC180", col = "red", lwd = 2)
#text(0.8, 0.2, paste("AUC =", round(auc(roc_data21), 2)), col = "blue", cex = 1.2)
means21 <- tapply(as.numeric((unlist(RNA_TPM_all_data1$"STRA6LP"))),RNA_TPM_all_data1[,1],mean)
means21
RNA_TPM_all_data1.1 <-RNA_TPM_all_data1
RNA_TPM_all_data1.1[,1] <- ifelse(RNA_TPM_all_data1[,1] == "0","healthy","LUAD")
bp21 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1.1$"STRA6LP"))) ~ RNA_TPM_all_data1.1[,1],
                data = RNA_TPM_all_data1.1,
                col = c("#8ECFC9","#FFBE7A"),
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA-LUAD")
bp21$stats[3,]
bp21

RNA_TPM_lnc_data_t <- t(RNA_TPM_lnc_data)
RNA_TPM_lnc_data1 <- cbind(rownames(RNA_TPM_lnc_data_t),RNA_TPM_lnc_data_t)
RNA_TPM_lnc_data1 <- merge(metaMatrix.RNA1,RNA_TPM_lnc_data1,by.x="sample",by.y="V1")
rownames(RNA_TPM_lnc_data1) <- RNA_TPM_lnc_data1[,1]
RNA_TPM_lnc_data1 <- RNA_TPM_lnc_data1[,-1]
RNA_TPM_lnc_data1[,1] <- ifelse(RNA_TPM_lnc_data1[,1] == "SolidTissueNormal",0,1)
roc_data22 <- roc(RNA_TPM_lnc_data1[,1],as.numeric((unlist(RNA_TPM_lnc_data1$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data22, main = "TCGA-LUAD", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data22), 2)), col = "blue", cex = 1.2)
write.csv(RNA_TPM_lnc_data1,"E:/github/R/TCGA-LUAD.csv")

#TCGA-COAD 结肠癌
setwd("E:/github/data/TCGA-COAD")
path<-"E:/github/data/TCGA-COAD"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_coad = merge_TCGA(metadata=metaMatrix.RNA, 
                   path="RNAseq", 
                   data.type="RNAseq",
                   mRNA_expr_type="TPM",
                   symbol = T,
                   RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_coad = as.matrix(RNA_TPM_coad)
rownames(RNA_TPM_all_coad)=RNA_TPM_all_coad[,1]
RNA_TPM_all_exp_coad=RNA_TPM_all_coad[,3:ncol(RNA_TPM_all_coad)]
RNA_TPM_all_dimnames_coad=list(rownames(RNA_TPM_all_exp_coad),colnames(RNA_TPM_all_exp_coad))
RNA_TPM_all_data_coad=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_coad)),nrow=nrow(RNA_TPM_all_exp_coad),dimnames=RNA_TPM_all_dimnames_coad)
RNA_TPM_all_data_coad=avereps(RNA_TPM_all_data_coad)

RNA_TPM_all_data_t_coad <- t(RNA_TPM_all_data_coad)
metaMatrix.RNA2 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_coad <- cbind(rownames(RNA_TPM_all_data_t_coad),RNA_TPM_all_data_t_coad)
RNA_TPM_all_data1_coad <- merge(metaMatrix.RNA2,RNA_TPM_all_data1_coad,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_coad) <- RNA_TPM_all_data1_coad[,1]
RNA_TPM_all_data1_coad <- RNA_TPM_all_data1_coad[,-1]
RNA_TPM_all_data1_coad[,1] <- ifelse(RNA_TPM_all_data1_coad[,1] == "SolidTissueNormal",0,1)
roc_data23 <- roc(RNA_TPM_all_data1_coad[,1],as.numeric((unlist(RNA_TPM_all_data1_coad$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data23, main = "TCGA COAD", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data23), 2)), col = "blue", cex = 1.2)
means23 <- tapply(as.numeric((unlist(RNA_TPM_all_data1_coad$CD226))),RNA_TPM_all_data1_coad[,1],mean)
means23
RNA_TPM_all_data1_coad.1 <- RNA_TPM_all_data1_coad
RNA_TPM_all_data1_coad.1[,1] <- ifelse(RNA_TPM_all_data1_coad[,1] == "0","SolidTissueNormal","COAD")
bp23 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_coad.1$"STRA6LP"))) ~ RNA_TPM_all_data1_coad.1[,1],
                data = RNA_TPM_all_data1_coad.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA COAD")
bp23$stats[3,]
bp23

RNA_TPM2 <- t(RNA_TPM_coad)
RNA_TPM2 <- as.data.frame(RNA_TPM2)

#TCGA-LIHC 肝癌
setwd("E:/github/data/TCGA-LIHC")
path<-"E:/github/data/TCGA-LIHC"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_lihc = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_lihc = as.matrix(RNA_TPM_lihc)
rownames(RNA_TPM_all_lihc)=RNA_TPM_all_lihc[,1]
RNA_TPM_all_exp_lihc=RNA_TPM_all_lihc[,3:ncol(RNA_TPM_all_lihc)]
RNA_TPM_all_dimnames_lihc=list(rownames(RNA_TPM_all_exp_lihc),colnames(RNA_TPM_all_exp_lihc))
RNA_TPM_all_data_lihc=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_lihc)),nrow=nrow(RNA_TPM_all_exp_lihc),dimnames=RNA_TPM_all_dimnames_lihc)
RNA_TPM_all_data_lihc=avereps(RNA_TPM_all_data_lihc)

RNA_TPM_all_data_t_lihc <- t(RNA_TPM_all_data_lihc)
metaMatrix.RNA3 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_lihc <- cbind(rownames(RNA_TPM_all_data_t_lihc),RNA_TPM_all_data_t_lihc)
RNA_TPM_all_data1_lihc <- merge(metaMatrix.RNA3,RNA_TPM_all_data1_lihc,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_lihc) <- RNA_TPM_all_data1_lihc[,1]
RNA_TPM_all_data1_lihc <- RNA_TPM_all_data1_lihc[,-1]
RNA_TPM_all_data1_lihc[,1] <- ifelse(RNA_TPM_all_data1_lihc[,1] == "SolidTissueNormal",0,1)
#AC022400.3 = "BMS1P4-AGAP5"
roc_data24 <- roc(RNA_TPM_all_data1_lihc[,1],as.numeric((unlist(RNA_TPM_all_data1_lihc$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data24, main = "TCGA LIHC", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data24), 2)), col = "blue", cex = 1.2)
means24 <- tapply(as.numeric((unlist(RNA_TPM_all_data1_lihc$CD226))),RNA_TPM_all_data1_lihc[,1],mean)
means24
RNA_TPM_all_data1_lihc.1 <- RNA_TPM_all_data1_lihc
RNA_TPM_all_data1_lihc.1[,1] <- ifelse(RNA_TPM_all_data1_lihc[,1] == "0","SolidTissueNormal","LIHC")
bp24 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_lihc.1$"STRA6LP"))) ~ RNA_TPM_all_data1_lihc.1[,1],
                data = RNA_TPM_all_data1_lihc.1,
                col = c("#8ECFC9","#FFBE7A"),
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA LIHC")
bp24$stats[3,]
bp24

RNA_TPM2 <- t(RNA_TPM_coad)
RNA_TPM2 <- as.data.frame(RNA_TPM2)

#TCGA-LUSC 肺鳞癌
setwd("E:/github/data/TCGA-LUSC")
path<-"E:/github/data/TCGA-LUSC"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_lucs = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_lucs = as.matrix(RNA_TPM_lucs)
rownames(RNA_TPM_all_lucs)=RNA_TPM_all_lucs[,1]
RNA_TPM_all_exp_lucs=RNA_TPM_all_lucs[,3:ncol(RNA_TPM_all_lucs)]
RNA_TPM_all_dimnames_lucs=list(rownames(RNA_TPM_all_exp_lucs),colnames(RNA_TPM_all_exp_lucs))
RNA_TPM_all_data_lucs=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_lucs)),nrow=nrow(RNA_TPM_all_exp_lucs),dimnames=RNA_TPM_all_dimnames_lucs)
RNA_TPM_all_data_lucs=avereps(RNA_TPM_all_data_lucs)

RNA_TPM_all_data_t_lucs <- t(RNA_TPM_all_data_lucs)
metaMatrix.RNA4 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_lucs <- cbind(rownames(RNA_TPM_all_data_t_lucs),RNA_TPM_all_data_t_lucs)
RNA_TPM_all_data1_lucs <- merge(metaMatrix.RNA4,RNA_TPM_all_data1_lucs,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_lucs) <- RNA_TPM_all_data1_lucs[,1]
RNA_TPM_all_data1_lucs <- RNA_TPM_all_data1_lucs[,-1]
RNA_TPM_all_data1_lucs[,1] <- ifelse(RNA_TPM_all_data1_lucs[,1] == "SolidTissueNormal",0,1)
write.csv(RNA_TPM_all_data1_lucs,"E:/github/R/TCGA-LUSC.csv")
roc_data25 <- roc(RNA_TPM_all_data1_lucs[,1],as.numeric((unlist(RNA_TPM_all_data1_lucs$"CD274"))),levels=c("0", "1"))
plot(roc_data25, main = "AC108010.1", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data25), 2)), col = "blue", cex = 1.2)
means25 <- tapply(as.numeric((unlist(RNA_TPM_all_data1_lucs$"STRA6LP"))),RNA_TPM_all_data1_lucs[,1],mean)
means25
RNA_TPM_all_data1_lucs.1 <- RNA_TPM_all_data1_lucs
RNA_TPM_all_data1_lucs.1[,1] <- ifelse(RNA_TPM_all_data1_lucs[,1] == "0","healthy","LUSC")
bp25 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_lucs.1$"STRA6LP"))) ~ RNA_TPM_all_data1_lucs.1[,1],
                data = RNA_TPM_all_data1_lucs.1,
                col = c("#8ECFC9","#FFBE7A"),
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA-LUSC")
bp25$stats[3,]
bp25

#TCGA-STAD 胃癌
setwd("E:/github/data/TCGA-STAD")
path<-"E:/github/data/TCGA-STAD"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_stad = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_stad = as.matrix(RNA_TPM_stad)
rownames(RNA_TPM_all_stad)=RNA_TPM_all_stad[,1]
RNA_TPM_all_exp_stad=RNA_TPM_all_stad[,3:ncol(RNA_TPM_all_stad)]
RNA_TPM_all_dimnames_stad=list(rownames(RNA_TPM_all_exp_stad),colnames(RNA_TPM_all_exp_stad))
RNA_TPM_all_data_stad=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_stad)),nrow=nrow(RNA_TPM_all_exp_stad),dimnames=RNA_TPM_all_dimnames_stad)
RNA_TPM_all_data_stad=avereps(RNA_TPM_all_data_stad)

RNA_TPM_all_data_t_stad <- t(RNA_TPM_all_data_stad)
metaMatrix.RNA6 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_stad <- cbind(rownames(RNA_TPM_all_data_t_stad),RNA_TPM_all_data_t_stad)
RNA_TPM_all_data1_stad <- merge(metaMatrix.RNA6,RNA_TPM_all_data1_stad,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_stad) <- RNA_TPM_all_data1_stad[,1]
RNA_TPM_all_data1_stad <- RNA_TPM_all_data1_stad[,-1]
RNA_TPM_all_data1_stad[,1] <- ifelse(RNA_TPM_all_data1_stad[,1] == "SolidTissueNormal",0,1)
#AC022400.3 = "BMS1P4-AGAP5"
roc_data27 <- roc(RNA_TPM_all_data1_stad[,1],as.numeric((unlist(RNA_TPM_all_data1_stad$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data27, main = "TCGA STAD", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data27), 2)), col = "blue", cex = 1.2)
means27 <- tapply(as.numeric((unlist(RNA_TPM_all_data1_stad$"STRA6LP"))),RNA_TPM_all_data1_stad[,1],mean)
means27
RNA_TPM_all_data1_stad.1 <- RNA_TPM_all_data1_stad
RNA_TPM_all_data1_stad.1[,1] <- ifelse(RNA_TPM_all_data1_stad[,1] == "0","SolidTissueNormal","STAD")
bp27 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_stad.1$"STRA6LP"))) ~ RNA_TPM_all_data1_stad.1[,1],
                data = RNA_TPM_all_data1_stad.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA STAD")
bp27$stats[3,]
bp27

#TCGA-PAAD 胰腺癌
setwd("E:/github/data/TCGA-PAAD")
path<-"E:/github/data/TCGA-PAAD"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_paad = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_paad = as.matrix(RNA_TPM_paad)
rownames(RNA_TPM_all_paad)=RNA_TPM_all_paad[,1]
RNA_TPM_all_exp_paad=RNA_TPM_all_paad[,3:ncol(RNA_TPM_all_paad)]
RNA_TPM_all_dimnames_paad=list(rownames(RNA_TPM_all_exp_paad),colnames(RNA_TPM_all_exp_paad))
RNA_TPM_all_data_paad=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_paad)),nrow=nrow(RNA_TPM_all_exp_paad),dimnames=RNA_TPM_all_dimnames_paad)
RNA_TPM_all_data_paad=avereps(RNA_TPM_all_data_paad)

RNA_TPM_all_data_t_paad <- t(RNA_TPM_all_data_paad)
metaMatrix.RNA5 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_paad <- cbind(rownames(RNA_TPM_all_data_t_paad),RNA_TPM_all_data_t_paad)
RNA_TPM_all_data1_paad <- merge(metaMatrix.RNA5,RNA_TPM_all_data1_paad,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_paad) <- RNA_TPM_all_data1_paad[,1]
RNA_TPM_all_data1_paad <- RNA_TPM_all_data1_paad[,-1]
RNA_TPM_all_data1_paad[,1] <- ifelse(RNA_TPM_all_data1_paad[,1] == "SolidTissueNormal",0,1)
#AC022400.3 = "BMS1P4-AGAP5"
roc_data26 <- roc(RNA_TPM_all_data1_paad[,1],as.numeric((unlist(RNA_TPM_all_data1_paad$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data26, main = "TCGA PAAD", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data26), 2)), col = "blue", cex = 1.2)
means26 <- tapply(as.numeric((unlist(RNA_TPM_all_data1_paad$"BMS1P4-AGAP5"))),RNA_TPM_all_data1_paad[,1],mean)
means26
RNA_TPM_all_data1_paad.1 <- RNA_TPM_all_data1_paad
RNA_TPM_all_data1_paad.1[,1] <- ifelse(RNA_TPM_all_data1_paad[,1] == "0","SolidTissueNormal","PAAD")
bp26 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_paad.1$"STRA6LP"))) ~ RNA_TPM_all_data1_paad.1[,1],
                data = RNA_TPM_all_data1_paad.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA PAAD")
bp26$stats[3,]
bp26


#TCGA-GBM 多形成性胶质细胞瘤
setwd("E:/github/data/TCGA-GBM")
path<-"E:/github/data/TCGA-GBM"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_GBM = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_GBM = as.matrix(RNA_TPM_GBM)
rownames(RNA_TPM_all_GBM)=RNA_TPM_all_GBM[,1]
RNA_TPM_all_exp_GBM=RNA_TPM_all_GBM[,3:ncol(RNA_TPM_all_GBM)]
RNA_TPM_all_dimnames_GBM=list(rownames(RNA_TPM_all_exp_GBM),colnames(RNA_TPM_all_exp_GBM))
RNA_TPM_all_data_GBM=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_GBM)),nrow=nrow(RNA_TPM_all_exp_GBM),dimnames=RNA_TPM_all_dimnames_GBM)
RNA_TPM_all_data_GBM=avereps(RNA_TPM_all_data_GBM)

RNA_TPM_all_data_t_GBM <- t(RNA_TPM_all_data_GBM)
metaMatrix.RNA7 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_GBM <- cbind(rownames(RNA_TPM_all_data_t_GBM),RNA_TPM_all_data_t_GBM)
RNA_TPM_all_data1_GBM <- merge(metaMatrix.RNA7,RNA_TPM_all_data1_GBM,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_GBM) <- RNA_TPM_all_data1_GBM[,1]
RNA_TPM_all_data1_GBM <- RNA_TPM_all_data1_GBM[,-1]
RNA_TPM_all_data1_GBM[,1] <- ifelse(RNA_TPM_all_data1_GBM[,1] == "SolidTissueNormal",0,1)
roc_data33 <- roc(RNA_TPM_all_data1_GBM[,1],as.numeric((unlist(RNA_TPM_all_data1_GBM$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data33, main = "TCGA GBM", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data33), 2)), col = "blue", cex = 1.2)
means33 <- tapply(as.numeric((unlist(RNA_TPM_all_data1_GBM$"BMS1P4-AGAP5"))),RNA_TPM_all_data1_GBM[,1],mean)
means33
RNA_TPM_all_data1_GBM.1 <- RNA_TPM_all_data1_GBM
RNA_TPM_all_data1_GBM.1[,1] <- ifelse(RNA_TPM_all_data1_GBM[,1] == "0","SolidTissueNormal","GBM")
bp33 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_GBM.1$"STRA6LP"))) ~ RNA_TPM_all_data1_GBM.1[,1],
                data = RNA_TPM_all_data1_GBM.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA GBM")
bp33$stats[3,]
bp33

#TCGA-KIRC 肾透明细胞瘤
setwd("E:/github/data/TCGA-KIRC")
path<-"E:/github/data/TCGA-KIRC"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_KIRC = merge_TCGA(metadata=metaMatrix.RNA, 
                         path="RNAseq", 
                         data.type="RNAseq",
                         mRNA_expr_type="TPM",
                         symbol = T,
                         RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_KIRC = as.matrix(RNA_TPM_KIRC)
rownames(RNA_TPM_all_KIRC)=RNA_TPM_all_KIRC[,1]
RNA_TPM_all_exp_KIRC=RNA_TPM_all_KIRC[,3:ncol(RNA_TPM_all_KIRC)]
RNA_TPM_all_dimnames_KIRC=list(rownames(RNA_TPM_all_exp_KIRC),colnames(RNA_TPM_all_exp_KIRC))
RNA_TPM_all_data_KIRC=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_KIRC)),nrow=nrow(RNA_TPM_all_exp_KIRC),dimnames=RNA_TPM_all_dimnames_KIRC)
RNA_TPM_all_data_KIRC=avereps(RNA_TPM_all_data_KIRC)

RNA_TPM_all_data_t_KIRC <- t(RNA_TPM_all_data_KIRC)
metaMatrix.RNA8 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_KIRC <- cbind(rownames(RNA_TPM_all_data_t_KIRC),RNA_TPM_all_data_t_KIRC)
RNA_TPM_all_data1_KIRC <- merge(metaMatrix.RNA8,RNA_TPM_all_data1_KIRC,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_KIRC) <- RNA_TPM_all_data1_KIRC[,1]
RNA_TPM_all_data1_KIRC <- RNA_TPM_all_data1_KIRC[,-1]
RNA_TPM_all_data1_KIRC[,1] <- ifelse(RNA_TPM_all_data1_KIRC[,1] == "SolidTissueNormal",0,1)
roc_data34 <- roc(RNA_TPM_all_data1_KIRC[,1],as.numeric((unlist(RNA_TPM_all_data1_KIRC$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data34, main = "TCGA KIRC", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data34), 2)), col = "blue", cex = 1.2)
RNA_TPM_all_data1_KIRC.1 <- RNA_TPM_all_data1_KIRC
RNA_TPM_all_data1_KIRC.1[,1] <- ifelse(RNA_TPM_all_data1_KIRC[,1] == "0","SolidTissueNormal","KIRC")
bp34 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_KIRC.1$"STRA6LP"))) ~ RNA_TPM_all_data1_KIRC.1[,1],
                data = RNA_TPM_all_data1_KIRC.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA KIRC")
bp34$stats[3,]
bp34

#TCGA-HNSC 头颈鳞状细胞瘤
setwd("E:/github/data/TCGA-HNSC")
path<-"E:/github/data/TCGA-HNSC"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_HNSC = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_HNSC = as.matrix(RNA_TPM_HNSC)
rownames(RNA_TPM_all_HNSC)=RNA_TPM_all_HNSC[,1]
RNA_TPM_all_exp_HNSC=RNA_TPM_all_HNSC[,3:ncol(RNA_TPM_all_HNSC)]
RNA_TPM_all_dimnames_HNSC=list(rownames(RNA_TPM_all_exp_HNSC),colnames(RNA_TPM_all_exp_HNSC))
RNA_TPM_all_data_HNSC=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_HNSC)),nrow=nrow(RNA_TPM_all_exp_HNSC),dimnames=RNA_TPM_all_dimnames_HNSC)
RNA_TPM_all_data_HNSC=avereps(RNA_TPM_all_data_HNSC)

RNA_TPM_all_data_t_HNSC <- t(RNA_TPM_all_data_HNSC)
metaMatrix.RNA9 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_HNSC <- cbind(rownames(RNA_TPM_all_data_t_HNSC),RNA_TPM_all_data_t_HNSC)
RNA_TPM_all_data1_HNSC <- merge(metaMatrix.RNA9,RNA_TPM_all_data1_HNSC,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_HNSC) <- RNA_TPM_all_data1_HNSC[,1]
RNA_TPM_all_data1_HNSC <- RNA_TPM_all_data1_HNSC[,-1]
RNA_TPM_all_data1_HNSC[,1] <- ifelse(RNA_TPM_all_data1_HNSC[,1] == "SolidTissueNormal",0,1)
roc_data35 <- roc(RNA_TPM_all_data1_HNSC[,1],as.numeric((unlist(RNA_TPM_all_data1_HNSC$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data35, main = "TCGA HNSC", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data35), 2)), col = "blue", cex = 1.2)
RNA_TPM_all_data1_HNSC.1 <- RNA_TPM_all_data1_HNSC
RNA_TPM_all_data1_HNSC.1[,1] <- ifelse(RNA_TPM_all_data1_HNSC[,1] == "0","SolidTissueNormal","HNSC")
bp35 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_HNSC.1$"STRA6LP"))) ~ RNA_TPM_all_data1_HNSC.1[,1],
                data = RNA_TPM_all_data1_HNSC,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA HNSC")
bp35$stats[3,]
bp35

#TCGA-BLCA 脑低级别胶质瘤
setwd("E:/github/data/TCGA-BLCA")
path<-"E:/github/data/TCGA-BLCA"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_BLCA = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_BLCA = as.matrix(RNA_TPM_BLCA)
rownames(RNA_TPM_all_BLCA)=RNA_TPM_all_BLCA[,1]
RNA_TPM_all_exp_BLCA=RNA_TPM_all_BLCA[,3:ncol(RNA_TPM_all_BLCA)]
RNA_TPM_all_dimnames_BLCA=list(rownames(RNA_TPM_all_exp_BLCA),colnames(RNA_TPM_all_exp_BLCA))
RNA_TPM_all_data_BLCA=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_BLCA)),nrow=nrow(RNA_TPM_all_exp_BLCA),dimnames=RNA_TPM_all_dimnames_BLCA)
RNA_TPM_all_data_BLCA=avereps(RNA_TPM_all_data_BLCA)

RNA_TPM_all_data_t_BLCA <- t(RNA_TPM_all_data_BLCA)
metaMatrix.RNA13 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_BLCA <- cbind(rownames(RNA_TPM_all_data_t_BLCA),RNA_TPM_all_data_t_BLCA)
RNA_TPM_all_data1_BLCA <- merge(metaMatrix.RNA13,RNA_TPM_all_data1_BLCA,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_BLCA) <- RNA_TPM_all_data1_BLCA[,1]
RNA_TPM_all_data1_BLCA <- RNA_TPM_all_data1_BLCA[,-1]
RNA_TPM_all_data1_BLCA[,1] <- ifelse(RNA_TPM_all_data1_BLCA[,1] == "SolidTissueNormal",0,1)
roc_data38 <- roc(RNA_TPM_all_data1_BLCA[,1],as.numeric((unlist(RNA_TPM_all_data1_BLCA$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data38, main = "TCGA BLCA", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data38), 2)), col = "blue", cex = 1.2)
RNA_TPM_all_data1_BLCA.1 <- RNA_TPM_all_data1_BLCA
RNA_TPM_all_data1_BLCA.1[,1] <- ifelse(RNA_TPM_all_data1_BLCA[,1] == "0","SolidTissueNormal","BLCA")
bp38 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_BLCA.1$"STRA6LP"))) ~ RNA_TPM_all_data1_BLCA.1[,1],
                data = RNA_TPM_all_data1_BLCA.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA BLCA")
bp38$stats[3,]
bp38

#TCGA-UCEC 子宫内膜癌
setwd("E:/github/data/TCGA-UCEC")
path<-"E:/github/data/TCGA-UCEC"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_UCEC = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_UCEC = as.matrix(RNA_TPM_UCEC)
rownames(RNA_TPM_all_UCEC)=RNA_TPM_all_UCEC[,1]
RNA_TPM_all_exp_UCEC=RNA_TPM_all_UCEC[,3:ncol(RNA_TPM_all_UCEC)]
RNA_TPM_all_dimnames_UCEC=list(rownames(RNA_TPM_all_exp_UCEC),colnames(RNA_TPM_all_exp_UCEC))
RNA_TPM_all_data_UCEC=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_UCEC)),nrow=nrow(RNA_TPM_all_exp_UCEC),dimnames=RNA_TPM_all_dimnames_UCEC)
RNA_TPM_all_data_UCEC=avereps(RNA_TPM_all_data_UCEC)

RNA_TPM_all_data_t_UCEC <- t(RNA_TPM_all_data_UCEC)
metaMatrix.RNA14 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_UCEC <- cbind(rownames(RNA_TPM_all_data_t_UCEC),RNA_TPM_all_data_t_UCEC)
RNA_TPM_all_data1_UCEC <- merge(metaMatrix.RNA14,RNA_TPM_all_data1_UCEC,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_UCEC) <- RNA_TPM_all_data1_UCEC[,1]
RNA_TPM_all_data1_UCEC <- RNA_TPM_all_data1_UCEC[,-1]
RNA_TPM_all_data1_UCEC[,1] <- ifelse(RNA_TPM_all_data1_UCEC[,1] == "SolidTissueNormal",0,1)
roc_data39 <- roc(RNA_TPM_all_data1_UCEC[,1],as.numeric((unlist(RNA_TPM_all_data1_UCEC$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data39, main = "TCGA UCEC", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data39), 2)), col = "blue", cex = 1.2)
RNA_TPM_all_data1_UCEC.1 <- RNA_TPM_all_data1_UCEC
RNA_TPM_all_data1_UCEC.1[,1] <- ifelse(RNA_TPM_all_data1_UCEC[,1] == "0","SolidTissueNormal","UCEC")
bp39 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_UCEC.1$"STRA6LP"))) ~ RNA_TPM_all_data1_UCEC.1[,1],
                data = RNA_TPM_all_data1_UCEC.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA UCEC")
bp39$stats[3,]
bp39

#TCGA-PRAD 前列腺癌
setwd("E:/github/data/TCGA-PRAD")
path<-"E:/github/data/TCGA-PRAD"
{
  merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){
    
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
    if (data.type=='RNAseq') {
      message ('###############    正在进行提取，请稍后   ################')
      if(mRNA_expr_type=="STAR"){
        column=4
      }else if(mRNA_expr_type=="TPM"){
        column=7
      }else if(mRNA_expr_type=="FPKM"){
        column=8
      }else if(mRNA_expr_type=="FPKM_UQ"){
        column=9
      }
      plan(multisession)
      rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
        read.table(fl,skip=6,sep="\t")[,column]))
      ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
      gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
      type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
      index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
      rnaMatrix=rnaMatrix[index,]
      rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
      gene_symbol=gene_symbol[index]
      type=type[index]
      colnames(rnaMatrix) <- metadata$sample
      nSamples = ncol(rnaMatrix)
      nGenes = nrow(rnaMatrix)
      if(RNA_type){
        rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      
      if(symbol){
        rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
      }
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of genes: ', nGenes, '\n', sep=''))
      #返回最后的基因表达矩阵
      return (rnaMatrix)
      
    }else if (data.type=='miRNAs') { 
      message ('############### Merging miRNAs data ###############\n')
      mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
      mirs <- mirbase$V1
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      rownames(mirMatrix) <- mirbase$V2
      colnames(mirMatrix) <- metadata$sample
      mirMatrix[is.na(mirMatrix)] <- 0
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      return (mirMatrix)
    }else{ 
      stop('data type error!')
    }
  }
  
  filtermir <- function(fl) {
    expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"),]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
    expr <- expr[,-1]
    names(expr) <- mirs
    return(expr)
  }
  
  FilterDuplicate <- function(metadata) {
    filter <- which(duplicated(metadata[,'sample']))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  FilterSampleType <- function(metadata) {
    filter <- which(! metadata$sample_type %in% 
                      c('PrimaryTumor', 'SolidTissueNormal'))
    if (length(filter) != 0) {
      metadata <- metadata[-filter,]
    }
    message (paste('Removed', length(filter), 'samples', sep=' '))
    return (metadata)
  }
  metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
  names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
  metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
  metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
  metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}
#TPM
RNA_TPM_PRAD = merge_TCGA(metadata=metaMatrix.RNA, 
                          path="RNAseq", 
                          data.type="RNAseq",
                          mRNA_expr_type="TPM",
                          symbol = T,
                          RNA_type=T
)
#提取所有表达数据
RNA_TPM_all_PRAD = as.matrix(RNA_TPM_PRAD)
rownames(RNA_TPM_all_PRAD)=RNA_TPM_all_PRAD[,1]
RNA_TPM_all_exp_PRAD=RNA_TPM_all_PRAD[,3:ncol(RNA_TPM_all_PRAD)]
RNA_TPM_all_dimnames_PRAD=list(rownames(RNA_TPM_all_exp_PRAD),colnames(RNA_TPM_all_exp_PRAD))
RNA_TPM_all_data_PRAD=matrix(as.numeric(as.matrix(RNA_TPM_all_exp_PRAD)),nrow=nrow(RNA_TPM_all_exp_PRAD),dimnames=RNA_TPM_all_dimnames_PRAD)
RNA_TPM_all_data_PRAD=avereps(RNA_TPM_all_data_PRAD)

RNA_TPM_all_data_t_PRAD <- t(RNA_TPM_all_data_PRAD)
metaMatrix.RNA15 <- metaMatrix.RNA[,c(7,8)]
RNA_TPM_all_data1_PRAD <- cbind(rownames(RNA_TPM_all_data_t_PRAD),RNA_TPM_all_data_t_PRAD)
RNA_TPM_all_data1_PRAD <- merge(metaMatrix.RNA15,RNA_TPM_all_data1_PRAD,by.x="sample",by.y="V1")
rownames(RNA_TPM_all_data1_PRAD) <- RNA_TPM_all_data1_PRAD[,1]
RNA_TPM_all_data1_PRAD <- RNA_TPM_all_data1_PRAD[,-1]
RNA_TPM_all_data1_PRAD[,1] <- ifelse(RNA_TPM_all_data1_PRAD[,1] == "SolidTissueNormal",0,1)
roc_data40 <- roc(RNA_TPM_all_data1_PRAD[,1],as.numeric((unlist(RNA_TPM_all_data1_PRAD$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data40, main = "TCGA PRAD", col = "red", lwd = 2)
text(0.8, 0.2, paste("AUC =", round(auc(roc_data40), 2)), col = "blue", cex = 1.2)
RNA_TPM_all_data1_PRAD.1 <- RNA_TPM_all_data1_PRAD
RNA_TPM_all_data1_PRAD.1[,1] <- ifelse(RNA_TPM_all_data1_PRAD[,1] == "0","SolidTissueNormal","PRAD")
bp40 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_PRAD.1$"STRA6LP"))) ~ RNA_TPM_all_data1_PRAD.1[,1],
                data = RNA_TPM_all_data1_PRAD.1,
                col = "pink",
                xlab = "label",ylab = "Gene Expression",
                main = "TCGA PRAD")
bp40$stats[3,]
bp40

library(RColorBrewer)
display.brewer.all()
br_pal <- brewer.pal(12,"Set3")
br_pal
plot(roc_data40, main = "ROC", col = "#8DD3C7", lwd = 2)
plot(roc_data39, add = T, col = "#FDB462", lwd = 2)
plot(roc_data38, add = T, col = "#BEBADA", lwd = 2)
plot(roc_data35, add = T, col = "#FB8072", lwd = 2)
plot(roc_data34, add = T, col = "#80B1D3", lwd = 2)
plot(roc_data33, add = T, col = "#FFED6F", lwd = 2)
plot(roc_data25, add = T, col = "#B3DE69", lwd = 2)
plot(roc_data27, add = T, col = "#FCCDE5", lwd = 2)
plot(roc_data28, add = T, col = "#D9D9D9", lwd = 2)
plot(roc_data24, add = T, col = "#BC80BD", lwd = 2)
plot(roc_data23, add = T, col = "#CCEBC5", lwd = 2)
plot(roc_data22, add = T, col = "red", lwd = 2)
legend("bottomright",legend = c("LUSC","LUAD","PRAD","UCEC","BLCA","HNSC","KIRC","GBM","STAD","PAAD","LIHC","COAD"),col = c("#B3DE69","red","#8DD3C7","#FDB462","#BEBADA","#FB8072","#80B1D3","#FFED6F","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5"),lwd = 2)

par(mfrow = c(2,6))
bp41 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_lucs.1$"STRA6LP"))) ~ RNA_TPM_all_data1_lucs.1[,1],
                col = c("#8DD3C7","pink"),
                xlab = "",ylab = "",
                main = "LUSC",
                axes = F)
bp42 <- boxplot(as.numeric((unlist(RNA_TPM_lnc_data1$"STRA6LP"))) ~ RNA_TPM_lnc_data1[,1],
                col = c("#8DD3C7","pink"),
                xlab = "",ylab = "",
                main = "LUAD",
                axes = F)
bp43 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_GBM.1$"STRA6LP"))) ~ RNA_TPM_all_data1_GBM.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "GBM",
                axes = F)
bp44 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_PRAD.1$"STRA6LP"))) ~ RNA_TPM_all_data1_PRAD.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "PRAD",
                axes = F)
bp45 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_lihc.1$"STRA6LP"))) ~ RNA_TPM_all_data1_lihc.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "LIHC",
                axes = F)
bp46 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_UCEC.1$"STRA6LP"))) ~ RNA_TPM_all_data1_UCEC.1[,1],
                col = c("#8DD3C7","pink"),
                xlab = "",ylab = "",
                main = "UCEC",
                axes = F)
bp47 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_KIRC.1$"STRA6LP"))) ~ RNA_TPM_all_data1_KIRC.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "KIRC",
                axes = F)
bp48 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_paad.1$"STRA6LP"))) ~ RNA_TPM_all_data1_paad.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "PAAD",
                axes = F)
bp49 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_coad.1$"STRA6LP"))) ~ RNA_TPM_all_data1_coad.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "COAD",
                axes = F)
bp50 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_stad.1$"STRA6LP"))) ~ RNA_TPM_all_data1_stad.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "STAD",
                axes = F)
bp51 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_HNSC.1$"STRA6LP"))) ~ RNA_TPM_all_data1_HNSC.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "HNSC",
                axes = F)
bp52 <- boxplot(as.numeric((unlist(RNA_TPM_all_data1_BLCA.1$"STRA6LP"))) ~ RNA_TPM_all_data1_BLCA.1[,1],
                col = c("pink","#8DD3C7"),
                xlab = "",ylab = "",
                main = "BLCA",
                axes = F)

par(mfrow = c(1,1))
roc_data25 <- roc(RNA_TPM_all_data1_lucs[,1],as.numeric((unlist(RNA_TPM_all_data1_lucs$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data25, main = "LUSC", col = "blue", lwd = 2)
roc_data22 <- roc(RNA_TPM_lnc_data1[,1],as.numeric((unlist(RNA_TPM_lnc_data1$"STRA6LP"))),levels=c("0", "1"))
plot(roc_data22, main = "luad", col = "blue", lwd = 2)


