# Analysis from the paper #

## Simulation Pipeline - Figure 1

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

## Running synerOmics on real world datasets (ProCan and validation datasets) - Figure 2 and 3

For result 2 and 3 in the paper, synerOmics was run using the following parameters (same parameters as simulation tests):

```bash
./batch_array.sh --minRep = 2, --maxOrder = 2, --nTrees = 1500, --mtry = 500, --maxDepth = 82, --minNode = 5
```

All input data can be found in synerOmics/paper_data. For the ProCan-DepMapSanger (Primary) dataset, the following inputs were used:

```bash
x_input = “synerOmics/paper_data/ProCan-DepMapSanger_protein_matrix_6692_averaged.csv”
y_input = “synerOmics/paper_data/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv”
```

For the Sun et al., (breast cancer validation) dataset, the following inputs were used:

```bash
x_input = “synerOmics/paper_data/validation_protein_matrix.csv”
y_input = “synerOmics/paper_data/validation_drug_response.csv”
```

Source of datasets used in paper:
Primary dataset (ProCan-DepMapSanger): Proteomic data and cell line information was derived from https://depmap.org/portal/download/all/. IC50 data was derived from: https://doi.org/10.1016/j.ccell.2022.06.010 
Validation: Proteomic data, cell line information and phenotype data was derived from Sun, R. et al. Proteomic Dynamics of Breast Cancer Cell Lines Identifies Potential Therapeutic Protein Targets. Mol Cell Proteomics 22, 100602 (2023)
