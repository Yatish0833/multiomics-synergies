#!/bin/bash


if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <interaction_model> <heritability> <version> <order> <n_trees> <min_node_size> <max_depth> <mtry>"
    exit 1
fi

model=$1
her=$2
version=$3
order=$4
n_trees=$5
min_node_size=$6
max_depth=$7
mtry=$8

# paths
BASEDIR=$(pwd)
INPUT_DIR="${BASEDIR}/results/simulated-phenotypes-her-${her}_V${version}"
PROTEIN_MATRIX="${BASEDIR}/results/cleaned_ProCan-DepMapSanger_protein_matrix_6692_averaged.csv"
PHENOTYPE_FILE="${INPUT_DIR}/${model}_order_${order}_phenotypes.csv"
OUTPUT_DIR="${BASEDIR}/results/Ranger/${model}_order_${order}_phenotypes_ranger_results_her${her}_V${version}"

mkdir -p "$OUTPUT_DIR"

#Rscript
Rscript scripts/Run-ranger-analysis.R \
   --model "$model" \
   --her "$her" \
   --version "$version" \
   --order "$order" \
   --n-trees "$n_trees" \
   --min-node-size "$min_node_size" \
   --max-depth "$max_depth" \
   --mtry "$mtry" \
   --protein-matrix "${PROTEIN_MATRIX}" \
   --phenotype "${PHENOTYPE_FILE}" \
   --work-dir "${BASEDIR}"


if [ $? -ne 0 ]; then
   echo "Error: Ranger analysis failed"
   exit 1
fi

echo "Ranger analysis complete for ${model} model with heritability ${her}, version ${version}, order ${order}"

# Check true-variables occurance frequencies under each model 
JSON_FILE="${INPUT_DIR}/interacting_proteins.json"
if [ "$order" -eq 2 ]; then
    RESULTS_FILE="${OUTPUT_DIR}/variable_pair_counts.csv"
    python scripts/check_model_frequencies.py "$JSON_FILE" "$RESULTS_FILE" "$model" "$order"
else
    TRIPLET_FILE="${OUTPUT_DIR}/variable_triplet_counts.csv"
    PAIR_FILE="${OUTPUT_DIR}/variable_pair_counts.csv"
    python scripts/check_model_frequencies.py "$JSON_FILE" "$TRIPLET_FILE" "$model" "$order" "$PAIR_FILE"
fi

echo "Frequency check complete for ${model} model with heritability ${her}, version ${version}, order ${order}"
