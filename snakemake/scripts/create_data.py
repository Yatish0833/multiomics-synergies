import pandas as pd
import numpy as np

x = pd.read_csv(snakemake.input["train"])
test = pd.read_csv(snakemake.input["test"])
y = pd.read_csv(snakemake.input["response"], usecols=['drug_id','cell_line_name', "max_screening_conc",'ln_IC50'])

data_file = snakemake.output["train"]
test_file = snakemake.output["test"]

if snakemake.params["nan_removed"] != 0:
    t = snakemake.params["nan_removed"]
    n = len(x.Cell_Line) + len(test.Cell_Line)  
    # count missing over full data set?
    na_count = pd.concat([x, test]).set_index("Cell_Line").isna().sum(axis=0)
    mask = na_count.iloc[:].apply(lambda x: x > t * n)  # ture for all proteins, that have more nan than the threshold t
    x = x.drop(columns=na_count.index[mask])
    test = test.drop(columns=na_count.index[mask])



drug_list = y.drug_id.unique()
drug_list = np.sort(drug_list)
  # drugs must be sorted in the same way in snakemake and here
for i, d in enumerate(drug_list[:snakemake.params["max_drug"]]):
    xy = x.merge(y[y.drug_id==d], left_on='Cell_Line', right_on='cell_line_name')
    
    xy = xy.rename(columns={'ln_IC50':'label'}).drop(columns=['Cell_Line', 'cell_line_name','drug_id'])
    xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns]
    xy.to_csv(data_file[i], index=False)
    
    testXY = test.merge(y[y.drug_id==d], left_on='Cell_Line', right_on='cell_line_name')
    testXY = testXY.rename(columns={'ln_IC50':'label'}).drop(columns=['Cell_Line', 'cell_line_name','drug_id'])
    testXY.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in testXY.columns]
    testXY.to_csv(test_file[i], index=False)

del data_file, test_file



