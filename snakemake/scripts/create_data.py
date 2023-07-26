import pandas as pd
import numpy as np

x = pd.read_csv(snakemake.input["train"])
test = pd.read_csv(snakemake.input["test"])
y = pd.read_csv(snakemake.input["response"], usecols=['drug_id','cell_line_name', "max_screening_conc",'ln_IC50'])

data_file = snakemake.output["train"]
test_file = snakemake.output["test"]

drug_list = y.drug_id.unique()
drug_list = np.sort(drug_list)
  # drugs must be sorted in the same way in snakemake and here
for i, d in enumerate(drug_list[:snakemake.params["max_drug"]]):
    xy = x.merge(y[y.drug_id==d], left_on='Cell_Line', right_on='cell_line_name')
    xy = xy.fillna(0).rename(columns={'ln_IC50':'label'}).drop(columns=['Cell_Line', 'cell_line_name','drug_id'])
    xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns]
    xy.to_csv(data_file[i], index=False)
    
    testXY = test.merge(y[y.drug_id==d], left_on='Cell_Line', right_on='cell_line_name')
    testXY = testXY.fillna(0).rename(columns={'ln_IC50':'label'}).drop(columns=['Cell_Line', 'cell_line_name','drug_id'])
    testXY.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in testXY.columns]
    testXY.to_csv(test_file[i], index=False)

del data_file, test_file



