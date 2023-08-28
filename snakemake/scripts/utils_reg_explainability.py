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


def results_fit_to_df(results, ols, y, test_data):
    coeffs = results.params.tolist()
    ids = results.params.index.tolist()

    pred_y = results.predict()
    train_pear_R = pearsonr(pred_y, y)
    train_mse = np.square(y - pred_y).mean()

    test_y = test_data.label.to_list()
    pred_y = results.predict(test_data) 
    test_pear_R = pearsonr(pred_y, test_y)
    test_mse = np.square(test_y - pred_y).mean()

    results = pd.DataFrame(columns=['coef_id', 'coef', 'train_MSE', 'MSE', 'train_pearsonR', 'pearsonR'])
    results['coef_id'] = ids
    results['coef'] = coeffs
    results["train_MSE"] = train_mse
    results["MSE"] = test_mse
    results["train_pearsonR"] = train_pear_R.statistic
    results["pearsonR"] = test_pear_R.statistic
    return results


def regularized_regression(data, test_data, synergies, proteins, drug, alpha=0.05, filter=0, reg='sqrt_lasso'):
    # sig = fdrcorrection(synergies["P>|z|"], alpha=alpha, method='indep', is_sorted=False)[0]
    # synergies = synergies[sig].reset_index()
    sig_synergies = synergies[(np.abs(synergies.coef)>=filter) & (synergies.p_corrected < alpha)][["coef_id", "order"]]
    # protein_list = np.unique([y for x in sig_synergies.coef_id for y in x.split(":")]).tolist()  # all proteins used in significant synergies of that drug
    protein_list = proteins.coef_id[proteins.p_corrected < alpha].to_list()  # get this from single protein regression
    
    final_results = []
    if len(sig_synergies.order) == 0: 
        print("No sig synergies")
        return final_results

    
    # with each interration add more interaction information to the regression formula
    for o in [1,2,4]:
        included = list(set(sig_synergies.coef_id[sig_synergies.order<=o])) + protein_list
        
        formular = "label" + " ~ maxscreeningconc + "
        formular = formular + " + ".join(included)
        
        try:
            ols = smf.ols(formular,data=data)
        except Exception as inst:
            print(type(inst))
            print('error in Lasso, drug', drug, "order: ", o, "number of variables (+)", formular.count('+') )
            print(sig_synergies.shape, len(protein_list))
            break  # we do not need to check any higher order, if the lower order already fails due to recursion depth
            
        ols.raise_on_perfect_prediction = False #preventing the perfect separation error
        try:
            results = ols.fit_regularized(method=reg)
        except Exception as inst:
            print(type(inst))
            print('error in Lasso, drug', drug, "order: ", o, "number of variables (+)", formular.count('+') )
            print(sig_synergies.shape, len(protein_list))
            break
            
        results = results_fit_to_df(results, ols, data["label"].to_list(), test_data)
        results["order"] = o
        results["drug"] = drug
        results["n_prot"] = len(protein_list) 
        results["n_obs"] = len(data.label)
        results["n_feat"] = len(included)
        results["n_syn"] = formular.count(':')  # this value can be wrong since we don't only have 2nd order synergies
        results["fit"] = reg
        final_results.append(results)
    return final_results
    
    