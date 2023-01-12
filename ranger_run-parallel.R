
#install.packages('optparse')

#!/usr/bin/env Rscript
library(optparse)
library(rpart)
library(ggplot2)
library(reshape2)
library(ranger)
#library(getTreeranger)

option_list = list(
  make_option(c("-n", "--n-trees"), dest = "nTrees", type = "integer", default = 1000,
              help = "input the number of trees to train the model [default= %default]", metavar = "character"),

  make_option(c("-t", "--mtry"), dest = "mtry", type = "integer", default = NULL,
              help = "input the number of variables to possibly split at in each node [default= %default]", metavar = "character"),
  
  make_option(c("-s", "--min-node-size"), dest = "minNodeSize", type = "integer", default = NULL,
              help = "min node size [default= %default]", metavar = "numeric"),

  make_option(c("-d", "--max-depth"), dest = "maxDepth", type = "integer", default = NULL,
              help = "max depth [default= %default]", metavar = "numeric"),

  make_option(c("-w", "--work-dir"), dest = "workDir", type = "character", default = "./",
              help = "working directory [default= %default]", metavar = "character"),

  make_option(c("-c", "--split"), dest = "split", type = "integer", default = 0,
              help = "split number for the temporal working directory [default= %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)

#nTrees <- 1000
#DIR
#setwd("/Users/reg032/workspace/procan")
setwd(args$workDir)
set.seed(13)
#DIR
X_df <- read.csv(paste0("tmp",args$split,"/data.csv"))


# Run ranger
ra <-ranger(label ~ ., X_df, num.trees = args$nTrees, mtry=args$mtry, min.node.size=args$minNodeSize, max.depth= args$maxDepth, importance = 'impurity')

#Saving importances
#DIR
write.csv(ra$variable.importance, paste0("tmp",args$split,'/importances.csv'),row.names = FALSE)

# Saving the tree info 
for (i in 1:args$nTrees){
  #DIR
  write.csv(treeInfo(ra,i),paste0("tmp",args$split,"/tree",i,".csv"), row.names = FALSE)
}

# print('a')
