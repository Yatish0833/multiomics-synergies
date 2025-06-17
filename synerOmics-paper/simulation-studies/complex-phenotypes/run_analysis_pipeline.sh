#!/bin/bash

# Generate data for different heritability values
bash scripts/generate_data.sh

# Array of interaction models
#interaction_models=("Interaction_only" "Modifier_protein" "No_Interaction" "Dominant" "Redundant" "Synergistic" "Cubic" "Sigmoidal" "Exponential")
interaction_models=("complex")

# Array of heritability values
heritability_values=(10 30 50 70 100)
#heritability_values=(100)

# Number of trees for Ranger
n_trees=1500

# Min node size and max depth for Ranger
min_node_size=5 # 5 seems reasonable
max_depth=82  # 0 means no limit, sqrt of number of predictors

# mtry for Ranger (set to sqrt(number of predictors) by default)
mtry=500  #

# Loop through heritability values and interaction models
for her in "${heritability_values[@]}"; do
    for model in "${interaction_models[@]}"; do
        for version in $(seq 1 100); do
            for order in 2; do
                echo "Running analysis for model: $model, heritability: $her, version: V$version, order: $order"
                bash Run-ranger-R-script.sh "$model" "$her" "$version" "$order" "$n_trees" "$min_node_size" "$max_depth" "$mtry"
            done
        done
    done
done

# Now after completing the above step run parse_frequency_information_and_plot.py



