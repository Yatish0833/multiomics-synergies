
#install.packages('optparse')

#!/usr/bin/env Rscript
library(optparse)
library(rpart)
library(ggplot2)
library(reshape2)
library(ranger)
#library(getTreeranger)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

nTrees <- 1000
#DIR
setwd("/Users/reg032/workspace/procan")
set.seed(13)
#DIR
X_df <- read.csv(paste0("tmp/tmp",args[1],"/data.csv"))


# Run ranger
ra <-ranger(label ~ ., X_df, num.trees = nTrees, importance = 'impurity')

#Saving importances
#DIR
write.csv(ra$variable.importance, paste0("tmp/tmp",args[1],'/importances.csv'),row.names = FALSE)

# Saving the tree info 
for (i in 1:nTrees){
  #DIR
  write.csv(treeInfo(ra,i),paste0("tmp/tmp",args[1],"/tree",i,".csv"), row.names = FALSE)
}

# print('a')
