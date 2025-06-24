setwd("E:/github")
path <- "E:/github"

install.packages("ggplot2")
install.packages("tidyverse")
install.packages("ggvenn")
install.packages("glmnet")
install.packages("survival")
install.packages("randomForestSRC") 

library(randomForestSRC)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(glmnet)
library(survival)
library(dplyr)

#读取生存信息的基因表达矩阵
train <- read.csv("E:/github/R/exosome.lncRNA.csv",header=T,row.names = 1)
count <- rowSums(train!=0)
count <- as.data.frame(count)
per <- count/19
per_num <- as.numeric(unlist(per))
percent <- sprintf("%.2f%%", per_num *100)
train <- cbind(train,count,percent)
train <- train[which(train$count >= 12),]

train <- train[,c(-20:-21)]
train <- t(train)
train <- as.data.frame(train)
#scaled_train <- scale(train)
os <- c(1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
os <- as.data.frame(os)
train <- cbind(train,os)
#提取基因表达矩阵
x <- as.matrix(train[, 1:136])
y <- as.matrix(train[, 137])
mod <- glmnet(x,y,family = "binomial",alpha = 1) 
set.seed(123)
cvmod <- cv.glmnet(x,y,family = "binomial",alpha = 1,nfolds = 18)   # 交叉验证
cvmod
plot(mod,label = T,lwd=2)
plot(mod,xvar = "lambda",label = T,lwd=2)
print(mod)
plot(cvmod)
#基于该图选择最佳的λ，一般可以采用两个内置函数实现cvmod$lambda.min和 cvmod$lambda.1se
cvmod$lambda.min
cvmod$lambda.1se
##基因筛选，采用coef函数即可，有相应参数的gene则被保留，采用λ使用的是lambda.min
coef.min <- coef(cvmod,s="lambda.min")
coef.min
codf.min <- as.data.frame.matrix(coef.min)
cvmod <- as.data.frame.matrix(cvmod)
  
write.table(codf.min,"E:/github/R/exo.lncRNA.csv",quote=F,sep = ",")

codf.min1 <- tibble::rownames_to_column(codf.min, var = "genename")
result <- codf.min1[codf.min1$s1 != 0,]
result <- as.data.frame(result)
result <- result[-1,-2]
result <- as.data.frame(result)
result1 <- t(result)
train1 <- t(train)
train1 <- as.data.frame.matrix(train1)
train1 <- tibble::rownames_to_column(train1,var = "result")
merge_data <- merge(train1,result,by="result")
merge_data <- t(merge_data)
merge_data <- as.data.frame(merge_data)
colnames(merge_data) <- merge_data[1,]
merge_data <- merge_data[-1,]
#os <- c(1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
#os <- as.data.frame(os)
merge_data <- cbind(os,merge_data)
write.table(merge_data,"E:/github/R/exo.lasso.csv",quote=F,sep = ",")

