#!/bin/bash


heritability_values=(10 30 50 70 100)
# Number of data versions to generate for each heritability value
num_versions=100 #2

for her in "${heritability_values[@]}"; do
    for version in $(seq 1 $num_versions); do
        echo "Generating data for heritability: $her, version: V$version"
        python Complex-phenotypes-generation.py --percent_explainability $her --version V$version
    done
done

echo "Data generation complete."