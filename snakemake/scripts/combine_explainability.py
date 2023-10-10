import pandas as pd
import numpy as np

import os



for model in snakemake.input["explained"]:
    df = pd.read_csv(model, sep='\t')
    if df.shape[0] == 0: continue   # only works if one model has a result, otherwise we dont get an output file (can just create file first, but empty?)
    if os.path.exists(snakemake.output["exp_models"]):
        df.to_csv(snakemake.output["exp_models"], index=False, sep='\t', mode='a', header=False)
    else:
        df.to_csv(snakemake.output["exp_models"], index=False, sep='\t',mode='a')
        
    if df.fit[1] == 'normal':  # don't want to use the header(don't wanna test how it is done), i think 0 should be fine too
        df = df[["train_pseudo_r2", "train_adj_r2", "train_MSE", "MSE", "train_pearsonR", "pearsonR", "drug", "order", "n_prot", "n_obs", "n_feat", "n_syn", "config", "fit"]].drop_duplicates(ignore_index=True)
    else:
        df = df[["train_MSE", "MSE", "train_pearsonR", "pearsonR", "drug", "order", "n_prot", "n_obs", "n_feat", "n_syn", "config", "fit"]].drop_duplicates(ignore_index=True)
        
    if os.path.exists(snakemake.output["exp_scores"]):
        df.to_csv(snakemake.output["exp_scores"], index=False, sep='\t', mode='a', header=False)
    else:
        df.to_csv(snakemake.output["exp_scores"], index=False, sep='\t',mode='a')

        
for model in snakemake.input["proteins"]:
    df = pd.read_csv(model)
    if os.path.exists(snakemake.output["exp_proteins"]):
        df.to_csv(snakemake.output["exp_proteins"], index=False, sep='\t', mode='a', header=False)
    else:
        df.to_csv(snakemake.output["exp_proteins"], index=False, sep='\t', mode='a')
        
        
del df