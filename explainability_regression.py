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

sys.setrecursionlimit(50000)

working_dir = 'julika/tmp/'

print("loading data")

# File inputs
x_input = 'data/ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx'
y_input = 'data/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv'
interaction_input = 'downstream_analysis/final_results_all.tsv'
# interaction_input = 'julika/tmp/final3_aux.tsv'

print("loading data", interaction_input)

#One row per cell line
x = pd.read_excel(x_input, engine='openpyxl').drop(columns=['Project_Identifier'])
c = [a.replace(';','.') for a in x.columns]
x.columns = c

y = pd.read_csv(y_input, usecols=['drug_id','cell_line_name','ln_IC50'])

tested_regression = pd.read_csv(interaction_input, sep="\t", index_col=False,  names=["coef_id", "coef" ,"std err" ,"z", "P>|z|", "[0.025", 	"0.975]", "converged", "pseudo_r2", "standard_fitting", "snps", "order", "drug"])

# filter for interactions that were activly tested in syneromics 
tested_mask = tested_regression.order == tested_regression.coef_id.str.count(":")+1
interactions = tested_regression[tested_mask]

drug_list = np.unique(interactions.drug)

print(len(drug_list), "drugs are considered")

## functions for regression

# to save regression results
def results_fit_to_df(results, ols, y):
    coeffs = results.params.tolist()
    pvals = results.pvalues.tolist()
    pseudo_r2 = results.rsquared
    adj_r2 = results.rsquared_adj
    tvals = results.tvalues.tolist()
    cint_low = results.conf_int()[0].tolist()
    cint_high = results.conf_int()[1].tolist()

    pred_y = ols.predict(results.params)
    pear_R = pearsonr(pred_y, y)
    mse = np.square(y - pred_y).mean()
    
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
    results["MSE"] = mse
    results['[0.025'] = cint_low
    results['0.975]'] = cint_high
    results["pearsonR"] = pear_R.statistic
    return results


def regression_per_drug(interactions, drug_list, alpha=0.05, filter=0):
    final_results = []
    j = 0
    for i,d in enumerate(drug_list):
        xy = x.merge(y[y["drug_id"] == d], left_on='Cell_Line', right_on='cell_line_name')
        xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns] 

        drug_interactions = interactions[interactions.drug == d].drop_duplicates(ignore_index=True)
        sig = fdrcorrection(drug_interactions["P>|z|"], alpha=alpha, method='indep', is_sorted=False)[0]
        drug_interactions = drug_interactions[sig]
        sig_interactions = drug_interactions[np.abs(drug_interactions.coef)>=filter][["coef_id", "order"]]
        protein_list = np.unique([y for x in sig_interactions.coef_id for y in x.split(":")]).tolist()  # all proteins used in significant synergies of that drug
        
        data = xy[["lnICfa"] + protein_list].fillna(0)
        # with each interration add more interaction information to the regression formula
        for o in range(1,5):
            # print("Order: ", o)
            formular = "lnICfa" + " ~ "
            formular = formular + " + ".join(protein_list + sig_interactions.coef_id[sig_interactions.order<=o].to_list())

            try:
                ols = smf.ols(formular,data=data)
            except:
                print('error in OLS, drug', d, "order: ", o, "number of variables (+)", formular.count('+') )
                print(sig_interactions.shape, len(protein_list))
                j += 1        
                break  # we do not need to check any higher order, if the lower order already fails due to recursion depth
            ols.raise_on_perfect_prediction = False #preventing the perfect separation error
            results = ols.fit(disp=False, maxiter=1000) #method prevents singular matrix
            results = results_fit_to_df(results, ols, data["lnICfa"].to_list())
            results["order"] = o
            results["drug"] = d
            results["n_prot"] = len(protein_list) 
            final_results.append(results)
    final_results = pd.concat(final_results)
    print(j/(len(drug_list)), " % of drugs failed OLS at some point")        
    print("The largest amount of proteins used for a regression was: ", final_results.n_prot.max())
    return final_results


# main part
# 
print("baseline regression")
baseline = regression_per_drug(interactions, drug_list)
baseline.to_csv(working_dir+"baseline_results.csv", index=False)

print("|coef| > 0.5 regression")
coef_05 = regression_per_drug(interactions, drug_list, filter=0.5)
coef_05.to_csv(working_dir+"coef_05_results.csv", index=False)

print("combinig tables")
coef_05["version"] = "|Coef| > 0.5"
baseline["version"] = "baseline"

combined = pd.concat([baseline, coef_05], ignore_index=True) # TODO
combined = combined[["MSE", "pseudo_r2", "adj_r2", "pearsonR", "order", "drug", "version"]].drop_duplicates(ignore_index=True)
combined.to_csv(working_dir+"combined_performance.csv", index=False)

sns.violinplot(data=combined, x="MSE", y="order", hue="version", cut=0, orient="h" )
plt.savefig(working_dir+"img/mse_combined.png")

del baseline, coef_05, combined
