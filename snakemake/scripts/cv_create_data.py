import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

x = pd.read_csv(snakemake.input["x"]).set_index('Cell_Line')
y = pd.read_csv(snakemake.input["response"], usecols=['CHEMBL','cell_line_name', "max_screening_conc",'ln_IC50'])

x.columns = [''.join([chr(int(y)+97) if b.isnumeric() else b for b in a.replace('_','').replace('.','')]) for a in x.columns]
y.columns = [''.join([chr(int(y)+97) if b.isnumeric() else b for b in a.replace('_','').replace('.','')]) for a in y.columns]

data_files = snakemake.output["folds"]

if snakemake.params["nan_removed"] != 0:
    t = snakemake.params["nan_removed"]
    n = x.shape[0]
    na_count = pd.concat([x, test]).isna().sum(axis=0)
    mask = na_count.iloc[:].apply(lambda x: x > t * n)  # ture for all proteins, that have more nan than the threshold t
    x = x.drop(columns=na_count.index[mask])
    test = test.drop(columns=na_count.index[mask])

    
def shift(x, oldLow=-5, oldHigh=16, newLow=0, newHigh=21):
    return ((x - oldLow) / (oldHigh - oldLow) * (newHigh - newLow) + newLow )

x = x.apply(shift, raw=True)

kf = KFold(n_splits=5, random_state=42, shuffle=False)

folds = []

for (_, test_index) in kf.split(X):
    folds.append(x.iloc[test_index,:])

drug_list = y.CHEMBL.unique()
drug_list = np.sort(drug_list)

for j, x in enumerate(folds):
  # drugs must be sorted in the same way in snakemake and here
    for i, d in enumerate(drug_list[:snakemake.params["maxdrug"]]):
        xy = x.merge(y[y.CHEMBL==d], left_on='Cell_Line', right_on='celllinename')
        xy = xy.fillna(0).rename(columns={'lnICfa':'label'}).drop(columns=['celllinename','CHEMBL'])
        
        xy.to_csv(datafiles[5*i + j], index=False)  # need to have the correct expand in snakemake




