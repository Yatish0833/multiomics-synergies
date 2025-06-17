#!/bin/bash

### Workflow steps ####
# Step 1: Generate data for different heritability values using the data generation script 
# Step 2: Set interaction models and heritability values 
# Step 3: Set parameters for the Ranger analysis (number of trees, min node size, max depth, mtry)
# Step 4: For each her-interaction-version-order combination, run the Ranger analysis script with the appropriate arguments
# Step 5:  Optionally, run the script to parse frequency information and plot results



# Generate data for different heritability values - needs the protein expression data
bash scripts/generate_data.sh

# interaction models considered
interaction_models=("Interaction_only" "Modifier_protein" "No_Interaction" "Dominant" "Redundant" "Synergistic" "Cubic" "Sigmoidal" "Exponential")
#interaction_models=("Interaction_only")

# Heritability values to test
heritability_values=(10 30 50 70 100)
#heritability_values=(100)

# tree params
n_trees=1500 # Number of trees for Ranger
min_node_size=5 # 5 seems reasonable
max_depth=82  # 0 means no limit, sqrt of number of predictors
mtry=500  # mtry for Ranger (set to sqrt(number of predictors) by default)

# Loop through heritability values and interaction models
for her in "${heritability_values[@]}"; do
    for model in "${interaction_models[@]}"; do
        for version in $(seq 1 10); do
            for order in 2; do #change accordingly for order 3
                echo "Running analysis for model: $model, heritability: $her, version: V$version, order: $order"
                bash scripts/Run-ranger-R-script.sh "$model" "$her" "$version" "$order" "$n_trees" "$min_node_size" "$max_depth" "$mtry"
            done
        done
    done
done

# Now after completing the above step run parse_frequency_information_and_plot.py

python scripts/parse_frequency_information_and_plot.py 


