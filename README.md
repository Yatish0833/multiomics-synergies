# **synerOmics**

_synerOmics_ is a tool that aims to identify interactions between continuous features that may have a greater-than-additive effect on a continuous phenotype. We have demonstrated its use case on the ProCan-DepMapSanger dataset (Goncalves _et al.,_ 2022), and a validation dataset (Sun _et al.,_ 2023), however it can be applied to any continous ‘omic data.

![synerOmics_readme_figure](https://github.com/user-attachments/assets/3fbdf255-6bf4-414f-b003-7ce82fe435b3)



_synerOmics step 1_ utilises the R package _ranger_ to build regression trees and generates a count table for how frequently two features occur as parent-child nodes. _synerOmics step 2_ then takes a user-supplied cutoff of occurrences and performs OLS regression through python’s _statsmodels_, identifying features that show “synergistic” effects on phenotypes. 


## Running synerOmics

Clone the repository and navigate to `/synerOmics` directory.

```bash
git clone https://github.com/csiro-internal/synerOmics
cd synerOmics 
```

## Option 1: bash script for array job

In `batch_array.sh`, change working directory "workingDir".

You can specify synerOmics and ranger parameters as follows:

```bash
--minRep (default =2), --maxOrder (default = 2), --nTrees (default = 1000), --mtry (default = 100), --maxDepth (default = 0), --minNode (default = 5)
```

_maxOrder_ is set to 2 by default. synerOmics has the ability parse parent-child1-child2 (order = 3) up to parent-child1-child_n_ (order = _n_) for _n_th order, however it has not been sufficiently tested through our synthetic data to guarantee reliable results. In data_prep_parallel.py, change x_input for your input feature matrix, and y_input for your phenotype file.

Run `batch_array.sh`:

```bash
chmod +x ./batch_array.sh
sbatch batch_array.sh
```
## Option 2: snakemake

Navigate to `/synerOmics/snakemake` directory. Change parameters in config file `config.json`:

```bash
train: train_subset_matrix.csv
test: test_subset_matrix.csv
y: pc_drug_response.csv
test_y: pc_drug_response.csv
nan: 0
nTrees: 6000
mtry: 500
minNodeSize: 10
maxDepth: 15
maxOrder: 2
minReps: 2
exp_mem: 10000
alpha: 0.05  # currently parsed, can be turned into list but then explainability.py needs to be updated
filter: 0.2  # currently parsed, can be turned into list but then explainability.py needs to be updated
regularized: 0
```

You can run snakemake locally, but it can be also run on a server as a job through `snake_server.sh`. 

# Simulation Pipeline

In our paper, we also included a simulation pipeline for synthetic data. If you would like to run your own simulation experiments, the instructions have been detailed in the  `synerOmics-paper/` folder for Figure 1.
