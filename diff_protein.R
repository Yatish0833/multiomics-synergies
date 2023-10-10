library("data.table")
library(ggplot2)
library(dplyr)    # alternatively, this also loads %>%
library(edgeR)
library(Biobase)
setwd("C:/Users/wen068/OneDrive - CSIRO/Documents/multiomics-synergies")


meta <- read.csv('data/pc_drug_response.csv', sep = ',', header=TRUE)
meta <- meta[,c('tissue', 'cell_line_name')]
meta <- meta[!duplicated(meta),]
rownames(meta) <- meta$cell_line_name
meta$cell_line_name <- NULL
meta$breast <- as.integer(meta$tissue == 'Breast')


count <- read.csv('data/pc_protein.csv', row.names="Cell_Line") 
count <- t(count)

count_filtered <- count
#count_filtered <- count[rowSums(is.na(count)) < 784/2,];
#count_filtered <- count_filtered %>% na.omit();


all(rownames(meta) %in% colnames(count_filtered))
all(rownames(meta) == colnames(count_filtered))

meta <- meta[colnames(count_filtered),]
all(rownames(meta) == colnames(count_filtered))

count_pos <- 2 ** count_filtered  # some values are negative for some reason 

count_pos <- replace(count_pos, is.na(count_pos), 0) # I am filtering completely at the moment

boxplot(count_pos[1,] ~ meta[, 'breast'], main=rownames(count_pos)[1])





eset <- ExpressionSet(assayData = count_pos, 
                      phenoData = AnnotatedDataFrame(meta))

exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))

plotDensities(eset, legend=FALSE)


design <- model.matrix(~breast, data = pData(eset))


library(limma)

fit <- lmFit(eset, design)

fit <- eBayes(fit)

results <- decideTests(fit[, "breast"])
summary(results)

stats <- topTable(fit, sort.by = "P", n = Inf)

hist(stats[, 'P.Value']) #a lot of 0.

volcanoplot(fit, highlight = 10, names = rownames(fit))

data <- as.data.frame(stats)
data$proteins <- rownames(data)

data$diffexpressed <- "NO"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 
data$diffexpressed[data$logFC > 1 & data$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.05, set as "DOWN"
data$diffexpressed[data$logFC < -1 & data$adj.P.Val < 0.05] <- "DOWN"

data$delabel <- NA
data <- data[order(data$adj.P.Val, decreasing=FALSE),]
top20 <- Reduce(rbind, by(data, data$diffexpressed, head, 15))
data$delabel[data$proteins %in% top20$proteins[top20$diffexpressed != "NO"]] <- data$proteins[data$proteins %in% top20$proteins[top20$diffexpressed != "NO"]]

library(ggrepel)
# plot adding
ggplot(data=data, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "green")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
