#!/usr/bin/env Rscript

#install.packages('optparse')
library(optparse)
library(rpart)
library(ggplot2)
library(reshape2)
library(ranger)
#library(getTreeranger)

# option_list = list(
#   make_option(c("-n", "--n-trees"), dest = "nTrees", type = "integer", default = 1000,
#               help = "input the number of trees to train the model [default= %default]", metavar = "character"),

#   make_option(c("-t", "--mtry"), dest = "mtry", type = "integer", default = NULL,
#               help = "input the number of variables to possibly split at in each node [default= %default]", metavar = "character"),

#   make_option(c("-s", "--min-node-size"), dest = "minNodeSize", type = "integer", default = NULL,
#               help = "min node size [default= %default]", metavar = "numeric"),

#   make_option(c("-d", "--max-depth"), dest = "maxDepth", type = "integer", default = NULL,
#               help = "max depth [default= %default]", metavar = "numeric"),

#   make_option(c("-w", "--work-dir"), dest = "workDir", type = "character", default = "./",
#               help = "working directory [default= %default]", metavar = "character"),

#   make_option(c("-c", "--split"), dest = "split", type = "integer", default = 0,
#               help = "split number for the temporal working directory [default= %default]", metavar = "character")
# );

# opt_parser = OptionParser(option_list = option_list)
# args = parse_args(opt_parser)

#args$nTrees = 1000
#args$mtry = 100
#args$minNodeSize = 5
#args$maxDepth = 0
#args$split = 1
#args$workDir = "/Users/XXXX/workspace/procan/tmp/"

setwd(snakemake@params[["workDir"]])
set.seed(13)
X_df <- read.csv(snakemake@input[["train"]])
test_df <- read.csv(snakemake@input[["test"]])


# Run ranger
ra <-ranger(label ~ ., X_df, num.trees = snakemake@params[["nTrees"]], mtry=snakemake@params[["mtry"]], min.node.size=snakemake@params[["minNodeSize"]], max.depth= snakemake@params[["maxDepth"]], importance = 'impurity')

#Saving importances
print("extracting trees:")
print(Sys.time())
# write.csv(ra$variable.importance, file=snakemake@output[["imp"]], row.names = FALSE)


# Saving the tree info 
final_df<-treeInfo(ra,1)
final_df$tree <- 1
for (i in 2:snakemake@params[["nTrees"]]){
  #write.csv(treeInfo(ra,i),paste0("tmp",args$split,"/tree",i,".csv"), row.names = FALSE)
  tr<-treeInfo(ra,i)
  tr$tree <- i
  final_df<-rbind(final_df, tr)
}
print("saving trees")
print(Sys.time())
write.csv(final_df, file=snakemake@output[["trees"]], row.names = FALSE)

# MSE = ra$prediction.error
# OOB = ra$r.squared (R squared OOB by ranger definition)
#ra

test_y <- predict(ra, test_df)


# print("create data frame")
df <- data.frame (train_MSE  = c(ra$prediction.error),
                  train_OOB = c(ra$r.squared),
                  train_pearsonR = cor(ra$predictions, X_df$label, method='pearson'),
                  MSE = mean((test_df$label - test_y$predictions)^2),
                  pearsonR = cor(test_y$predictions, test_df$label, method = 'pearson'),
                  drug = snakemake@params[["drug"]]
)
print("saving performances")
print(Sys.time())
write.table(df, file=snakemake@output[["score"]], quote=FALSE, sep='\t')


