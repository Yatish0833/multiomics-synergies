#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Importing libraries
import pandas as pd
import numpy as np
import os
import sys
import shutil

import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from itertools import combinations

import matplotlib.pyplot as plt

from scipy.stats import pearsonr
import seaborn as sns
from random import sample, seed

import warnings
warnings.filterwarnings("ignore")
# In[10]:


working_dir = 'tmp/'

# File inputs
# x_input = 'data/ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx'
x_input = 'data/train_protein_matrix.csv'
y_input = '../data/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv'

test_input = 'data/test_protein_matrix.csv'


# In[36]:


x = pd.read_csv(x_input)
c = [a.replace('.','').replace('_','') for a in x.columns]
x.columns = c

x_test = pd.read_csv(test_input)
x_test.columns = c

y = pd.read_csv(y_input)[['drug_id','cell_line_name','ln_IC50', 'max_screening_conc']]


# In[37]:


drug_list = np.unique(y.drug_id)
protein_list = c[1:]


# In[39]:


def results_fit_to_df(results, ols, y, test_data):
    coeffs = results.params.tolist()
    pvals = results.pvalues.tolist()
    pseudo_r2 = results.rsquared
    # adj_r2 = results.rsquared_adj
    tvals = results.tvalues.tolist()
    cint_low = results.conf_int()[0].tolist()
    cint_high = results.conf_int()[1].tolist()

    pred_y = results.predict()
    train_pear_R = pearsonr(pred_y, y).statistic
    train_mse = np.square(y - pred_y).mean()

    if len(test_data.ln_IC50) < 8:
        test_pear_R = None
        test_mse = None
    else:    
        test_y = test_data.ln_IC50.to_list()
        pred_y = results.predict(test_data)
        test_pear_R = pearsonr(pred_y, test_y).statistic   # want 20 entries for this
        test_mse = np.square(test_y - pred_y).mean()

    
    try:
        results = results.summary()
    except:
        #ValueError: resids must contain at least 2 elements
        r = pd.DataFrame([1,2,3]) #dirty...
        r['z']='nan'
        return r
    converged = results.tables[0].data[5][1].strip()
    results = results.tables[1].data
    results = pd.DataFrame(results[1:], columns=['coef_id', 'coef', 'std err', 'z', 'P>|z|', '[0.025', '0.975]'])
    results['P>|z|'] = pvals
    results['z'] = tvals 
    results['coef'] = coeffs
    results['converged'] = converged
    results['pseudo_r2'] = pseudo_r2
    # results['adj_r2'] = adj_r2
    results["train_MSE"] = train_mse
    results["MSE"] = test_mse
    results['[0.025'] = cint_low
    results['0.975]'] = cint_high
    results["train_pearsonR"] = train_pear_R
    results["pearsonR"] = test_pear_R
    return results


regression_results = []
for d in drug_list:
    print(f'currently looking at drug {d}')
    c = 0
    xy = x.merge(y[y["drug_id"]==d], left_on='CellLine', right_on='cell_line_name')
    test_xy = x_test.merge(y[y["drug_id"] == d], left_on='CellLine', right_on='cell_line_name')
    if len(xy.CellLine) < 300*0.8: continue

    for p in protein_list:

        
        data = xy[["ln_IC50", "max_screening_conc", p]].dropna()
        test_data = test_xy[["ln_IC50", "max_screening_conc", p]].dropna()

        # print(f'For drug {d} and {p} there are {len(data[p])} and {len(test_data[p])} cell lines')

        if len(test_data[p]) < 8 or len(data[p]) < 8:
            # print(f'Not enough test cell lines for protein {p} and drug {d}, n = {len(test_data[p])}')
            c += 1
            continue

        formula = "ln_IC50 ~ max_screening_conc + " + p
        ols = smf.ols(formula,data=data)
        results = ols.fit(disp=False, maxiter=1000)
        results = results_fit_to_df(results, ols, data["ln_IC50"].to_list(), test_data)
        results["drug"] = d
        results["n_train"] = len(data[p])
        results["n_test"] = len(test_data[p])
        regression_results.append(results.loc[[2]])
        del data, test_data   
    print(f'For {c} proteins we did not have enough data')
    del xy, test_xy
regression_results = pd.concat(regression_results)


regression_results.to_csv("protein_regression.csv", index=False)

pd.pivot_table(regression_results, values='pearsonR', index=['drug'], columns=['coef_id'], aggfunc=np.sum, fill_value=0).to_csv("protein_regression_pearsonR.csv")





