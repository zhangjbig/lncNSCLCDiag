library(plyr)
library(readr)
setwd("E:/github")
path<-"E:/github"

Data.lnc <- read_csv("E:/github/R/exosome.lncRNA.csv",col_names=TRUE)
count <- rowSums(Data.lnc!=0)
count <- as.data.frame(count)
count <- count-1
Data.lnc <- cbind(Data.lnc,count)
Data.lnc <- Data.lnc[which(Data.lnc$count >= 9),]
Data.lnc <- as.data.frame(Data.lnc)
Data.lnc <- Data.lnc[!duplicated(Data.lnc[,1]),]
rownames(Data.lnc) <- Data.lnc[,1]
Data.lnc<- Data.lnc[,-1]
Data.lnc.2 <- read_csv("E:/github/R/0.58_exosome.lncRNA.csv") 
Data.lnc.2 <- as.data.frame(Data.lnc.2)
rownames(Data.lnc.2) <- Data.lnc.2[,1]
Data.lnc.2 <- Data.lnc.2[,-1]
Data.lnc.3 <- Data.lnc.2[,c(1,4,5)]
Data <- merge(Data.lnc.3,Data.lnc,by = "row.names",all = F)
Data[,1] <- substring(Data[,1],1,nchar(Data[,1])-4)
Data <- as.data.frame(Data)
Data <- Data[!duplicated(Data[,1]),]
write.csv(Data,"E:/github/R/0.58.exo.lnc.csv")

LUAD <- read_csv("E:/github/R/TCGA-LUAD.csv")
LUAD <- as.data.frame(LUAD)
rownames(LUAD) <- LUAD[,1]
LUAD <- LUAD[,-1]
LUAD.1 <- LUAD[,c(-1,-2)]
LUAD.1 <- t(LUAD.1)
count.1 <- rowSums(LUAD.1!=0)
count.1 <- as.data.frame(count.1)
LUAD.2 <- cbind(count.1,LUAD.1)
LUAD.2 <- LUAD.2[which(LUAD.2$count.1 >= 293),]
LUAD.2 <- LUAD.2[,-1]

LUSC <- read_csv("E:/github/R/TCGA-LUSC.csv")
LUSC <- as.data.frame(LUSC)
rownames(LUSC) <- LUSC[,1]
LUSC.1 <- LUSC[,c(-1,-2)]
LUSC.1 <- t(LUSC.1)
count.2 <- rowSums(LUSC.1!=0)
count.2 <- as.data.frame(count.2)
LUSC.2 <- cbind(count.2,LUSC.1)
LUSC.2 <- LUSC.2[which(LUSC.2$count.2 >= 276),]
LUSC.2 <- LUSC.2[,-1]

suppressMessages(library(limma))
LUAD.3 <- t(LUAD.2)
LUAD.3 <- cbind(LUAD[,1],LUAD.3)
group.1 <- LUAD[,1]
LUSC.3 <- t(LUSC.2)
LUSC.3 <- cbind(LUSC[,2],LUSC.3)
group.2 <- LUSC[,2]

design.1 <- model.matrix(~0+factor(group.1)) #对样本实验设置进行计算。设计矩阵
colnames(design.1) <- factor(c("Healthy","LUAD"))
rownames(design.1) <- colnames(LUAD.2)
design.2 <- model.matrix(~0+factor(group.2)) 
colnames(design.2) <- factor(c("Healthy","LUSC"))
rownames(design.2) <- colnames(LUSC.2)
contrasts.matrix.1 <- makeContrasts(LUAD-Healthy,
                                  levels = design.1)
contrasts.matrix.2 <- makeContrasts(LUSC-Healthy,
                                    levels = design.2)

fit_1 <- lmFit(LUAD.2,design.1) 
fit2_1 <- contrasts.fit(fit_1,contrasts.matrix.1) 
fit2_1 <- eBayes(fit2_1)
tempOutput1 <- topTable(fit2_1,coef = 1,n=Inf)
nrDEG1 = na.omit(tempOutput1)
nrDEG1 <- cbind(rownames(nrDEG1),nrDEG1)

fit_2 <- lmFit(LUSC.2,design.2) 
fit2_2 <- contrasts.fit(fit_2,contrasts.matrix.2) 
fit2_2 <- eBayes(fit2_2)               
tempOutput2 <- topTable(fit2_2,coef = 1,n=Inf)
nrDEG2 = na.omit(tempOutput2)
nrDEG2 <- cbind(rownames(nrDEG2),nrDEG2)

res_healthy_vs_NSCLC1<- topTable(fit2_1,coef = 1,n=Inf)
res_healthy_vs_NSCLC2<- topTable(fit2_2,coef = 1,n=Inf)
diffsig_healthy_vs_NSCLC1 <- res_healthy_vs_NSCLC1[with(res_healthy_vs_NSCLC1,(abs(logFC)>0.58& adj.P.Val<0.05)),]
diffsig_healthy_vs_NSCLC2 <- res_healthy_vs_NSCLC2[with(res_healthy_vs_NSCLC2,(abs(logFC)>0.58& adj.P.Val<0.05)),]

LUAD.4 <- diffsig_healthy_vs_NSCLC1[,c(1,5)]
LUAD.4 <- cbind(rownames(LUAD.4),LUAD.4)
overlap.1 <- merge(LUAD.4,Data,by.x = "rownames(LUAD.4)", by.y = "Row.names")
LUSC.4ealthy_vs_NSCLC2[,c(1,5)]
LUSC.4 <- cbind(rownames(LUSC.4),LUSC.4)
overlap.2 <- merge(LUSC.4,Data,by.x = "rownames(LUSC.4)", by.y = "Row.names")

res_healthy_vs_NSCLC1.1 <- cbind(rownames(res_healthy_vs_NSCLC1),res_healthy_vs_NSCLC1)
res_healthy_vs_NSCLC2.1 <- cbind(rownames(res_healthy_vs_NSCLC2),res_healthy_vs_NSCLC2) 

Data.lnc1 <- t(Data.lnc)
Data.lnc1 <- as.data.frame(Data.lnc1)
exosome <- cbind(Data.lnc1$`AL445524.1-201`,Data.lnc1$`MMP25-AS1-201`,Data.lnc1$`OR2A1-AS1-201`,Data.lnc1$`OVCH1-AS1-201`,
                 Data.lnc1$`RTCA-AS1-201`,Data.lnc1$`SNHG1-204`,Data.lnc1$`STRA6LP-201`,Data.lnc1$`TRG-AS1-201`,Data.lnc1$`TTC28-AS1-201`)
write.csv(exosome,"E:/github/R/exosome.csv") 

LUAD.exo <- cbind(LUAD[,1],LUAD$`AL445524.1`,LUAD$"MMP25-AS1",LUAD$"OR2A1-AS1",LUAD$"OVCH1-AS1",LUAD$"RTCA-AS1",LUAD$"SNHG1",LUAD$"STRA6LP",
                  LUAD$"TRG-AS1",LUAD$`TTC28-AS1`)
write.table(LUAD.exo,"E:/github/R/LUAD.csv",quote=F,sep = ",")
LUSC.exo <- cbind(LUSC[,2],LUSC$`AL445524.1`,LUSC$"MMP25-AS1",LUSC$"OR2A1-AS1",LUSC$"OVCH1-AS1",LUSC$"RTCA-AS1",LUSC$"SNHG1",LUSC$"STRA6LP",
                  LUSC$"TRG-AS1",LUSC$`TTC28-AS1`)
write.table(LUSC.exo,"E:/github/R/LUSC.csv",quote=F,sep = ",")


diffsig_healthy_vs_NSCLC1$Row.names <- rownames(diffsig_healthy_vs_NSCLC1)
diffsig_healthy_vs_NSCLC2$Row.names <- rownames(diffsig_healthy_vs_NSCLC2)
TCGA <- merge(diffsig_healthy_vs_NSCLC1,diffsig_healthy_vs_NSCLC2,by= "Row.names")
TCGA <- TCGA[,c(-3:-5,-7,-9:-11,-13)]
rownames(TCGA) <- TCGA[,1]
TCGA <- TCGA[,-1]
TCGA <- as.data.frame(TCGA)
TCGA1 <- subset(TCGA,(abs(logFC.x)>1 & abs(logFC.y)>1))
list <- rownames(TCGA1)
list <- as.data.frame(list)

LUAD.1.1 <- cbind(rownames(LUAD.1),LUAD.1)
LUSC.1.1 <- cbind(rownames(LUSC.1),LUSC.1)
list.1 <- merge(list,LUAD.1.1,by.x = "list",by.y = "V1")
list.1 <- merge(list.1,LUSC.1.1,by.x = "list",by.y = "V1" )
rownames(list.1) <- list.1[,1]
list.1 <- list.1[,-1]
list.2 <- t(list.1)
list.2 <- cbind(rownames(list.2),list.2)
LUAD.1.2 <- LUAD
LUAD.1.2 <- cbind(rownames(LUAD),LUAD)
LUAD.1.2 <- LUAD.1.2[,c(1,2)]
LUSC.1.2 <- LUSC
LUSC.1.2 <- LUSC.1.2[,c(1,2)]
colnames(LUSC.1.2)[1] <- "rownames(LUAD)"
TCGA2 <- rbind(LUAD.1.2,LUSC.1.2)
list.2 <- merge(TCGA2,list.2,by.x= "rownames(LUAD)",by.y ="V1")
rownames(list.2) <- list.2[,1]
list.2 <- list.2[,-1]
list.2 <- as.data.frame(list.2)
