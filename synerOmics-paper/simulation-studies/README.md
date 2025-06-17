# Simulation Pipeline

This repository contains a simulation pipeline designed to test Syneromics assumptions. Syneromics assumes that variables co-occur more frequently as parent-child nodes across the regression trees when they are interacting than when they are not. The pipeline generates synthetic drug IC50 phenotype data with assignable interaction patterns (model) using real protein expression data, and analyses them using Random Forest model (through Ranger R package) to create frequency count tables for the true variable combinations. In addition to standard phenotype generation (labelled simple phenotypes generation as they are generated under one randomly assigned interaction model), the pipeline also supports creating complex phenotypes, where the phenotype is influenced by multiple sets of proteins (each set has at most 2 proteins) and each set can be modeled with a different assumed interaction model for added complexity. The scripts for generating these complex phenotypes can be found inside the `complex-phenotypes` folder.

## Directory Structure

```
simulation/
|-- Syneromics-simulation-avoiding-correlated-proteins.py  # Main simulation script for simple phenotype
|-- Simulation-pipeline.sh                                 # workflow pipeline (complete)
|-- Run-ranger-R-script.sh                                  # Ranger analysis pipeline 
|--generate_data.sh                                   # script to generate datasets 
|-- Run-ranger-analysis.R                           #Rscript for running Ranger
|-- parse_frequency_information_and_plot.py #plotting helper functions 
|-- check_model_frequencies.py         # Helper function to generate occurance count
|-- complex-phenotypes/
    |-- Complex-phenotypes-generation.py             # Complex pheno generation script
    |-- Simulation-pipeline.sh                       # workflow pipeline
    |-- Run-ranger-R-script.sh                          # Ranger analysis pipeline 
    |-- generate_data.sh                                # script to generate datasets 
    |-- Run-ranger-analysis.R                           #Rscript for running Ranger
    |-- parse_frequency_information_and_plot.py #plotting helper functions 
    |-- check_model_frequencies.py       # Helper function to generate occurance count
```

## Setup requirements

- Python 3.11.5
- R (4.4.1) with the following packages:
  - ranger (0.16.0)
  - data.table (1.16.0)
  - ggplot2 (3.5.1)
- Required Python packages:
  - numpy (1.25.2)
  - pandas (2.1.1)
  - matplotlib (3.8.0)
  - seaborn (0.12.2)

## How to use the folder

1. **If we only want to generate simulated pheno data, run**:

   ```bash
   bash generate_data.sh
   ```

2. **If we want to run the complete pipeline, use**:

   ```bash
   bash Simulation-pipeline.sh #make sure directory paths are proper
   ```

3. **Run individual analysis (for each using each script)**:

   ```bash
   bash Run-ranger-R-script.sh <model> <heritability> <version> <order> <n_trees> <min_node_size> <max_depth> <mtry> #example run
   ```

### Interaction Models Considered

- Interaction_only
- Modifier_protein
- No_Interaction
- Dominant
- Redundant
- Synergistic
- Cubic
- Sigmoidal
- Exponential
- Complex (different protein sets and acting through different model stated above)

### Tested heritability Values

- 10%
- 30%
- 50%
- 70%
- 100%

### Random Forest Parameters (all changeable within the scripts)

- Number of trees: 1500
- Minimum node size: 5
- Maximum depth: 82
- mtry: 500

## Notes

- The simulation avoids correlated proteins (Pearson's correlation > 0.70) to ensure clean interaction detection
- Each model can be (and is) run with multiple versions (can be changed) for statistical robustness, this can be set using versions parameter
