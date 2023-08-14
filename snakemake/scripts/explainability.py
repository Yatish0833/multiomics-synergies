import pandas as pd
import numpy as np
import os
import sys
import shutil

import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import pearsonr

from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns

from utils_explainability import explain_regression, undo
from utils_reg_explainability import regularized_regression



sys.setrecursionlimit(50000)

data = pd.read_csv(snakemake.input["train"])
test = pd.read_csv(snakemake.input["test"])

synergies = pd.read_csv(snakemake.input["synergies"], sep='\t')

results = []
for alpha in snakemake.params["alpha"]:
    for filter in snakemake.params["filter"]:
        if snakemake.params["regularized"]:
            result = regularized_regression(data, test, synergies, snakemake.wildcards["drug"], alpha, filter, 'sqrt_lasso')
        else: 
            result = explain_regression(data, test, synergies, snakemake.wildcards["drug"], alpha, filter)
        if result:
            df = pd.concat(result)
            df["config"] = f'a: {alpha} f: {filter}'
            df["Proteins"] = ['+'.join([y if y in ["Intercept", "maxscreeningconc"] else undo(x) for y in x.split(':')]) for x in df.coef_id]
            results.append(df)
if results:        
    results = pd.concat(results)
    results.to_csv(snakemake.output[0], sep='\t', index=False)
else:
    pd.DataFrame(columns = ["train_pseudo_r2", "train_adj_r2", "train_MSE", "MSE", "train_pearsonR", "pearsonR", "drug", "order", "config"]).to_csv(snakemake.output[0], sep='\t', index=False)