#!/usr/bin/env Rscript


library(optparse)
library(ranger)
library(data.table)
library(parallel)

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
  make_option(c("-p", "--protein-matrix"), type="character", default=NULL, 
              help="Path to protein matrix CSV file", metavar="FILE"),
  make_option(c("-y", "--phenotype"), type="character", default=NULL, 
              help="Path to phenotype CSV file", metavar="FILE"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="Model name", metavar="character"),
  make_option(c("-r", "--her"), type="character", default=NULL, 
              help="HER value", metavar="character"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="Version number", metavar="character"),
  make_option(c("-o", "--order"), type="character", default=NULL, 
              help="Order value", metavar="character")
)


opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)


setwd(args$workDir)

output_dir_name <- sprintf("%s_order_%s_phenotypes_ranger_results_her%s_V%s", args$model, args$order, args$her, args$version)
args$output_dir <- file.path("results", "Ranger", output_dir_name)


cat("Creating output directory:", args$output_dir, "\n")
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(args$output_dir)) {
  stop("Failed to create output directory: ", args$output_dir)
}

num_threads <- max(1, detectCores() - 1)
cat("Using", num_threads, "threads for parallel processing\n")


cat("Reading input files...\n")
protein_matrix <- read.csv(args$`protein-matrix`, header=T, row.names = 1)
phenotype <- read.csv(args$phenotype, header = T, row.names = 1)


phenotype_file_name <- basename(args$phenotype)
phenotype_name <- tools::file_path_sans_ext(phenotype_file_name)


interaction_model <- args$model
order <- args$order
her <- args$her 
version <- args$version

# Check if the number of rows match
if (nrow(protein_matrix) != nrow(phenotype)) {
  stop("Error: The number of rows in the protein matrix does not match the number of rows in the phenotype data.")
}

# Check if the row names match
if (!all(rownames(protein_matrix) == rownames(phenotype))) {
  stop("Error: The row identifiers in the protein matrix do not match those in the phenotype data.")
}

cat("Data validation passed: Protein matrix and phenotype data are correctly aligned.\n")

# Combine protein_matrix and phenotype for Ranger function
combined_df <- cbind(ln_IC50 = phenotype, protein_matrix)

cat("Dimensions of combined_df:", dim(combined_df), "\n")

# Run ranger with parallel processing
cat("Training random forest model...\n")
ra <- ranger(
  ln_IC50 ~ .,
  data = combined_df,
  num.trees = args$nTrees,
  mtry = args$mtry,
  min.node.size = args$minNodeSize,
  max.depth = args$maxDepth,
  importance = 'impurity',
  seed = 14,
  node.stats = TRUE,
  num.threads = num_threads
)


rm(protein_matrix, phenotype, combined_df)
gc()


cat("Current working directory:", getwd(), "\n")
cat("Output directory path:", args$output_dir, "\n")


tryCatch({
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
}, error = function(e) {
  cat("Error creating directory:", conditionMessage(e), "\n")
  cat("Full path attempted:", file.path(getwd(), args$output_dir), "\n")
  stop("Execution halted due to directory creation error")
})


write_output <- function(data, filename) {
  setnames(data, old = names(data), new = gsub("\\.", ";", names(data)))
  fwrite(data, file.path(args$output_dir, filename))
}

cat("Writing variable importances...\n")
write_output(data.table(variable = gsub("\\.", ";", names(ra$variable.importance)), 
                        importance = ra$variable.importance), 
             "variable_importances.csv")

print(ra)

# Function to get only parent-child pairs - exactly matching original script's approach
get_parent_child_pairs <- function(tree_info) {
  setDT(tree_info)
  
  tree_info[, splitvarName := gsub("\\.", ";", splitvarName)]
  
  parent_child_pairs <- tree_info[!is.na(splitvarName), .(
    var1 = splitvarName,
    var2 = c(tree_info$splitvarName[match(leftChild, nodeID)],
             tree_info$splitvarName[match(rightChild, nodeID)])
  )]
  
  parent_child_pairs <- parent_child_pairs[!is.na(var1) & !is.na(var2) & var1 != var2]
  
  return(parent_child_pairs)
}

# Get tree information, count parent-child pairs, and aggregate tree info
cat("Processing trees...\n")
num_trees <- ra$num.trees
parent_child_pairs <- data.table()
aggregated_tree <- data.table()

for (i in 1:num_trees) {
  if (i %% 500 == 0) {
    cat("Processing tree", i, "of", num_trees, "\n")
  }

  tree_info <- data.table(treeInfo(ra, i))
  
  pairs <- get_parent_child_pairs(tree_info)
  
  if (nrow(pairs) > 0) {
    pairs[, tree := i]
    parent_child_pairs <- rbind(parent_child_pairs, pairs)
  }
  

  tree_info[, tree := i]
  aggregated_tree <- rbind(aggregated_tree, tree_info)
  
  if (i %% 200 == 0) {
    gc()
  }
}

# Count pair occurrences
cat("Calculating pair statistics...\n")
pair_counts <- parent_child_pairs[, .(
  count = .N
), by = .(
  var1 = pmin(var1, var2),
  var2 = pmax(var1, var2)
)][order(-count)]
cat("Top 10 most frequent parent-child variable pairs:\n")
print(head(pair_counts, 10))

write_output(pair_counts, "variable_pair_counts.csv")
write_output(parent_child_pairs, "detailed_parent_child_pairs.csv")

cat("Aggregating tree statistics...\n")
aggregated_stats <- aggregated_tree[, .(
  usage_count = sum(!is.na(splitvarName)),
  avg_split_value = mean(splitval, na.rm = TRUE),
  min_split_value = min(splitval, na.rm = TRUE),
  max_split_value = max(splitval, na.rm = TRUE),
  avg_node_samples = mean(numSamples, na.rm = TRUE),
  avg_split_stat = mean(splitStat, na.rm = TRUE)
), by = .(splitvarName)]

aggregated_stats <- aggregated_stats[order(-usage_count)]

cat("Top most frequently used variables:\n")
print(head(aggregated_stats, 10))

write_output(aggregated_stats, "aggregated_tree_stats.csv")

cat("Analyzing variable usage by tree level...\n")
level_usage <- aggregated_tree[!is.na(splitvarName), .(
  usage_count = .N
), by = .(splitvarName, level = floor(log2(nodeID + 1)))]

top_vars_by_level <- level_usage[order(level, -usage_count), head(.SD, 5), by = level]

cat("Top variables used at each tree level:\n")
print(top_vars_by_level)

write_output(top_vars_by_level, "variable_usage_by_level.csv")

cat("Analysis complete. Output files written to:", args$output_dir, "\n")
cat("Files generated: variable_importances.csv, variable_pair_counts.csv, aggregated_tree_stats.csv, variable_usage_by_level.csv, detailed_parent_child_pairs.csv\n")