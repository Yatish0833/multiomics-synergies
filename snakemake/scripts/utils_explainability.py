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
    pvals = results.pvalues.tolist()
    pseudo_r2 = results.rsquared
    adj_r2 = results.rsquared_adj
    tvals = results.tvalues.tolist()
    cint_low = results.conf_int()[0].tolist()
    cint_high = results.conf_int()[1].tolist()

    pred_y = results.predict()
    train_pear_R = pearsonr(pred_y, y)
    train_mse = np.square(y - pred_y).mean()

    test_y = test_data.label.to_list()
    pred_y = results.predict(test_data)
    test_pear_R = pearsonr(pred_y, test_y)
    test_mse = np.square(test_y - pred_y).mean()

    # print(type(train_pear_R), train_pear_R)
    
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
    results['train_pseudo_r2'] = pseudo_r2
    results['train_adj_r2'] = adj_r2
    results["train_MSE"] = train_mse
    results["MSE"] = test_mse
    results['[0.025'] = cint_low
    results['0.975]'] = cint_high
    results["train_pearsonR"] = train_pear_R.statistic
    results["pearsonR"] = test_pear_R.statistic
    return results


def explain_regression(xy, test_xy, synergies, drug, alpha=0.05, filter=0):
    sig = fdrcorrection(synergies["P>|z|"], alpha=alpha, method='indep', is_sorted=False)[0]
    synergies = synergies[sig].reset_index()
    sig_synergies = synergies[np.abs(synergies.coef)>=filter][["coef_id", "order"]]
    protein_list = np.unique([y for x in sig_synergies.coef_id for y in x.split(":")]).tolist()  # all proteins used in significant synergies of that drug
    data = xy[["label", "maxscreeningconc"] + protein_list]
    test_data = test_xy[["label", "maxscreeningconc"] + protein_list]
    del xy, test_xy
    # print(protein_list)
    final_results = []
    # with each interration add more interaction information to the regression formula
    for o in [1,2,4]:
        excluded = sig_interactions.coef_id[sig_interactions.order<=o].drop_duplicates(ignore_index=True).to_list() + protein_list
        included = []
        print("number of possible variables:", len(excluded))

        while excluded:
            p_vals = []
            # forward
            for variable in excluded:
                formular = "lnICfa" + " ~ maxscreeningconc + "
                formular = formular + " + ".join(included + [variable])
        
                results = smf.ols(formular, data=data).fit(disp=False, maxiter=1000)
                p_vals.append(results.pvalues[variable])
            best_p = np.min(p_vals)
            if best_p < alpha:
                i = np.argmin(p_vals)
                included.append(excluded.pop(i))
            else:
                print(f'Terminated with best p-value: {best_p}')
                break
            # backward
            formular = "lnICfa" + " ~ maxscreeningconc + "
            formular = formular + " + ".join(included)
            results = smf.ols(formular, data=data).fit(disp=False, maxiter=1000)
            worst_p = np.max(results.pvalues)
            i = np.argmax(results.pvalues) - 2  #  (intercept and maxscreen)
            if i <= 0: continue
            if worst_p > alpha:
                del included[i]
        
        formular = "lnICfa" + " ~ maxscreeningconc + "
        formular = formular + " + ".join(included)
        
        try:
            ols = smf.ols(formular,data=data)
        except Exception as inst:
            print(type(inst))
            print('error in OLS, drug', d, "order: ", o, "number of variables (+)", formular.count('+') )
            print(sig_interactions.shape, len(protein_list))
            j += 1
            
            break  # we do not need to check any higher order, if the lower order already fails due to recursion depth
        ols.raise_on_perfect_prediction = False #preventing the perfect separation error
        results = ols.fit(maxiter=100) #method prevents singular matrix
        results = results_fit_to_df(results, ols, data["lnICfa"].to_list(), test_data)
        results["order"] = o
        results["drug"] = d
        results["n_prot"] = len(protein_list) 
        results["n_obs"] = len(data.lnICfa)
        results["n_feat"] = len(included)
        results["n_syn"] = formular.count(':')
        final_results.append(results)
    return final_results
    
    