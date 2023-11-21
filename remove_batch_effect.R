library(sva)
library(ggplot2)
library(dplyr)
library(ggfortify)
setwd("C:/Users/wen068/OneDrive - CSIRO/Documents/multiomics-synergies")


# need to read in pc and bc data, combine, set batch columns (transpose)
# do one pca before and after (prcomp) - allows na?



pc <- read.csv('data/pc_protein.csv', row.names="Cell_Line") 
bc <- read.csv('data/bc_protein.csv', row.names="Cell_Line") 

pc_bc <- cbind(t(pc), t(bc))
pc_bc <- replace(pc_bc, is.na(pc_bc), 0)
pc_bc <- 2 ** pc_bc


pc_batch <- rep(0, 784)
bc_batch <- rep(1, 76)

batch <- c(pc_batch, bc_batch)


bc_pc_pca <- prcomp(t(pc_bc),scale=TRUE)



bc_pc_pca <- data.frame(PC1 = bc_pc_pca$x[,1], PC2 = bc_pc_pca$x[,2], Group = c(rep('ProCan', 784), rep('Westlake', 76)))

ggplot(bc_pc_pca,aes(x=PC1,y=PC2,col=Group))+
  geom_point(size=3,alpha=0.7)+ #Size and alpha just for fun
  scale_color_manual(values = c('#007377', '#6D2077'))+ #your colors here
  theme_classic()


adjusted <- ComBat(pc_bc , batch=batch)
# these leads to errors   covar_mod = NULL, model.matrix(~factor(meta.oxford$Condition1)),full_mod=TRUE)

adjusted_pca <- prcomp(t(adjusted),scale=TRUE)

adjusted_pca <- data.frame(PC1 = adjusted_pca$x[,1], PC2 = adjusted_pca$x[,2], Group = c(rep('ProCan', 784), rep('Westlake', 76)))

ggplot(adjusted_pca,aes(x=PC1,y=PC2,col=Group))+
  geom_point(size=3,alpha=0.7)+ #Size and alpha just for fun
  scale_color_manual(values = c('#007377', '#6D2077'))+ #your colors here
  theme_classic()

adjusted <- log2(adjusted)

write.csv(t(adjusted[,1:784]), 'data/pc_protein_nbatch.csv')

write.csv(t(adjusted[,785:860]), 'data/bc_protein_nbatch.csv')


