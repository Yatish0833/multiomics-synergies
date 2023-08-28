import pandas as pd
import numpy as np

import os
import statsmodels.formula.api as smf
from scipy.stats import pearsonr

from utils_explainability import undo


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

    if len(test_data.label) < 8:
        test_pear_R = None
        test_mse = None
    else:    
        test_y = test_data.label.to_list()
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


xy = pd.read_csv(snakemake.input["train"])
test_xy = pd.read_csv(snakemake.input["test"])
protein_list = xy.columns[:-2].to_list() + ['DOSIS']

drop_nans = True

regression_results = []
for p in protein_list:
    if p == 'DOSIS':
        data = xy[["label", "maxscreeningconc"]]
        test_data = test_xy[["label", "maxscreeningconc"]]
        formula = "label ~ maxscreeningconc"
    else:
        data = xy[["label", "maxscreeningconc", p]]
        test_data = test_xy[["label", "maxscreeningconc", p]].dropna()    
        formula = "label ~ maxscreeningconc + " + p

    if drop_nans:
        data = data.loc[~(data==0.0).all(axis=1)]  # is this the correct call?
        test_data = test_data.loc[~(test_data==0.0).all(axis=1)]
    if (data.shape[0] < 8) or (test_data.shape[0] < 8): continue
    
    ols = smf.ols(formula,data=data)
    results = ols.fit(disp=False, maxiter=1000)
    results = results_fit_to_df(results, ols, data["label"].to_list(), test_data)
    results["drug"] = snakemake.params['drug']
    results["n_train"] = len(data['label'])
    results["n_test"] = len(test_data['label'])
    results['protein'] = undo(p)
    regression_results.append(results)
    del data, test_data   

regression_results = pd.concat(regression_results)

regression_results.to_csv(snakemake.output[0], index=False)
del xy, test_xy

