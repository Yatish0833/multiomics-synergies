import pandas as pd
import numpy as np

import os


for synergie in snakemake.input["synergies"]:
    if os.path.exists(snakemake.output["synergies"]):
        pd.read_csv(synergie, sep='\t').to_csv(snakemake.output["synergies"], index=False, sep='\t', mode='a', header=False)
    else:
        pd.read_csv(synergie, sep='\t').to_csv(snakemake.output["synergies"], index=False, sep='\t',mode='a')
        

for tree_score in snakemake.input["rf_score"]:
    if os.path.exists(snakemake.output["rf_scores"]):
        pd.read_csv(tree_score, sep='\t').to_csv(snakemake.output["rf_scores"], index=False, sep='\t', mode='a', header=False)
    else:
        pd.read_csv(tree_score, sep='\t').to_csv(snakemake.output["rf_scores"], index=False, sep='\t',mode='a')


for model in snakemake.input["explained"]:
    df = pd.read_csv(model, sep='\t')
    if os.path.exists(snakemake.output["exp_models"]):
        df.to_csv(snakemake.output["exp_models"], index=False, sep='\t', mode='a', header=False)
    else:
        df.to_csv(snakemake.output["exp_models"], index=False, sep='\t',mode='a')

    df = df[["train_pseudo_r2", "train_adj_r2", "train_MSE", "MSE", "train_pearsonR", "pearsonR", "drug", "order", "version"]].drop_duplicates(ignore_index=True)
    if os.path.exists(snakemake.output["exp_scores"]):
        df.to_csv(snakemake.output["exp_scores"], index=False, sep='\t', mode='a', header=False)
    else:
        df.to_csv(snakemake.output["exp_scores"], index=False, sep='\t',mode='a')

del df