# Importing libraries
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

working_dir = 'julika/tmp/'

print("loading data")

# File inputs
x_input = 'data/ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx'
y_input = 'data/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv'
interaction_input = 'downstream_analysis/final_results_all.tsv'

#One row per cell line
x = pd.read_excel(x_input, engine='openpyxl').drop(columns=['Project_Identifier'])
c = [a.replace(';','.') for a in x.columns]
x.columns = c

y = pd.read_csv(y_input)[['drug_id','cell_line_name','ln_IC50']]

tested_interactions = pd.read_csv(interaction_input, sep="\t", index_col=False,  names=["coef_id", "coef" ,"std err" ,"z", "P>|z|", "[0.025", 	"0.975]", "converged", "pseudo_r2", "standard_fitting", "snps", "order", "drug"])

proteins = tested_interactions[tested_interactions["coef_id"].str.contains("HUMAN")].copy()  # remove Intercept from data
drug_list = np.unique(tested_interactions["drug"])


## functions for regression

# to save regression results
def results_fit_to_df(results):
    coeffs = results.params.tolist()
    pvals = results.pvalues.tolist()
    pseudo_r2 = results.rsquared
    adj_r2 = results.rsquared_adj
    tvals = results.tvalues.tolist()
    cint_low = results.conf_int()[0].tolist()
    cint_high = results.conf_int()[1].tolist()

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
    results['adj_r2'] = adj_r2
    results['[0.025'] = cint_low
    results['0.975]'] = cint_high
    return results


def regression_per_drug(proteins, drug_list, alpha=0.05, filter=0):
    final_results = []
    j = 0
    for i,d in enumerate(drug_list[7:9]):  # [[45,76]]
        xy = x.merge(y[y["drug_id"]==d], left_on='Cell_Line', right_on='cell_line_name')
        xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns] 
        
        interactions = proteins[proteins["drug"]==d]
        sig = fdrcorrection(interactions["P>|z|"], alpha=alpha, method='indep', is_sorted=False)[0]
        interactions = interactions[sig]
    
        sig_inter = interactions[np.abs(interactions.coef)>=filter][["coef_id", "order"]] # all sig interaction independet from order
        sig_inter["coef_id"] = sig_inter["coef_id"].str.strip("\[01\]")  # remove certain characters
        sig_inter.loc[~sig_inter["coef_id"].str.contains(":"), 'order'] = 1  # denote single proteins as first interaction ( by it self)
        sig_inter = sig_inter.drop_duplicates(ignore_index=True)
        prot_list = [str(y).rstrip("\[01\]") for x in sig_inter["coef_id"] for y in x.split(":")]  # all proteins connected to any sig information (needed for the regression)
        prot_list = np.unique(prot_list).tolist()
        n_prot = len(sig_inter["coef_id"][sig_inter["order"]<=1])

        data = xy[["lnICfa"] + prot_list].fillna(0)
        # with each interation add more interaction information to the regression formula
        for o in range(1,5):
            # print("Order: ", o)
            formular = "lnICfa" + " ~ "
            formular = formular + " + ".join(sig_inter["coef_id"][sig_inter["order"]<=o].drop_duplicates(ignore_index=True))

            try:
                ols = smf.ols(formular,data=data)
            except:
                # print('error in OLS, drug', d, "number of variables (+)", formular.count('+') )
                j += 1
                
                break  # we do not need to check any higher order, if the lower order already fails due to recursion depth
            ols.raise_on_perfect_prediction = False #preventing the perfect separation error
            results = ols.fit(disp=False, maxiter=1000) #method prevents singular matrix
            pred_y = ols.predict(results.params)
            pear_R = pearsonr(pred_y, data["lnICfa"].to_list())
            results = results_fit_to_df(results)
            results["order"] = o
            results["drug"] = d
            results["pearsonR"] = pear_R.statistic
            results["n_prot"] = n_prot
            results["MSE"] = np.square(data["lnICfa"].to_list() - pred_y).mean()
            final_results.append(results)
    final_results = pd.concat(final_results)
    print(j/(len(drug_list)), " % of drugs failed OLS at some point")        
    print("The largest amount of proteins used for a regression was: ", final_results.n_prot.max())
    return final_results



# main part
# 
print("baseline regression")
baseline = regression_per_drug(proteins, drug_list)
baseline.to_csv(working_dir+"baseline_results.csv", index=False)

print("|coef| > 0.5 regression")
coef_05 = regression_per_drug(proteins, drug_list, filter=0.5)
coef_05.to_csv(working_dir+"coef_05_results.csv", index=False)

print("combinig tables")
r2_5["batch"] = "|Coef| > 0.5"
r2_05["batch"] = "baseline"

combined = pd.concat([r2_5, r2_05], ignore_index=True)
combined = r_combined[["MSE", "pseudo_r2", "adj_r2", "pearsonR", "order", "drug", "batch"]].drop_duplicates(ignore_index=True)
combined.to_csv(working_dir+"combined_performance.csv", index=False)

del baseline, coef_05, combined