import pandas as pd
import numpy as np
import os
import sys
import shutil

import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection, multipletests
from scipy.stats import pearsonr

from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns

from utils_explainability import explain_regression, undo
from utils_reg_explainability import regularized_regression

import argparse


sys.setrecursionlimit(50000)

parser = argparse.ArgumentParser()


parser.add_argument('train')
parser.add_argument('test')
parser.add_argument('synergies')
parser.add_argument('proteins')
parser.add_argument('output')
parser.add_argument('drug')
parser.add_argument('-a', nargs=1, help='alpha', type=float)
parser.add_argument('-f', nargs=2, help='filter', type=float)
parser.add_argument('-r', type=int, default=0)

args = parser.parse_args()


data = pd.read_csv(args.train)
test = pd.read_csv(args.test)

correction = 'fdr_bh'

synergies = pd.read_csv(args.synergies, sep='\t')
synergies['p_corrected'] = multipletests(synergies['P>|z|'], method=correction)[1]  # make into config option as well

proteins = pd.read_csv(args.proteins)
proteins = proteins[proteins.coef_id.str.contains('HUMAN')].copy().dropna()
proteins['p_corrected'] = multipletests(proteins['P>|z|'], method=correction)[1]
print(f'P value correction done using {correction}')

results = []
for alpha in args.a:
    for filter in args.f:
        if args.r == 1:
            print('regularization applied')
            result = regularized_regression(data, test, synergies, proteins, args.drug, alpha, float(filter))  #    DEFAULT is lasso
        else: 
            result = explain_regression(data, test, synergies, proteins, args.drug, alpha, float(filter))
        if result:
            df = pd.concat(result)
            df["config"] = f'a: {alpha} f: {filter}'
            df["Proteins"] = ['+'.join([y if y in ["Intercept", "maxscreeningconc"] else undo(x) for y in x.split(':')]) for x in df.coef_id]
            results.append(df)
        print(f'config a: {alpha} f: {filter} is finished')
if results:        
    results = pd.concat(results)
    results.to_csv(args.output, sep='\t', index=False)
else:
    pd.DataFrame(columns = ["train_pseudo_r2", "train_adj_r2", "train_MSE", "MSE", "train_pearsonR", "pearsonR",'n_prot', 'n_obs ', 'n_feat', 'n_syn', "drug", "order", "config", "fit"]).to_csv(args.output, sep='\t', index=False)