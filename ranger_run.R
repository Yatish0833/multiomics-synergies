
#install.packages('optparse')

#!/usr/bin/env Rscript
library(optparse)
library(rpart)
library(ggplot2)
library(reshape2)
library(ranger)
#library(getTreeranger)


nTrees <- 1000
setwd("/Users/reg032/workspace/procan")
set.seed(13)
X_df <- read.csv("tmp/data.csv")


# Run ranger
ra <-ranger(label ~ ., X_df, num.trees = nTrees, importance = 'impurity')

#Saving importances
write.csv(ra$variable.importance, 'output/importances.csv',row.names = FALSE)

# Saving the tree info 
for (i in 1:nTrees){
  write.csv(treeInfo(ra,i),paste0("output/tree",i,".csv"), row.names = FALSE)
}

# print('a')
