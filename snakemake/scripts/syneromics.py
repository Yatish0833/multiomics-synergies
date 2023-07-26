import pandas as pd
import numpy as np
import os
import sys
import argparse
import shutil

import statsmodels.formula.api as smf
from itertools import combinations

from utils_syneromics import get_interactions, from_table_to_json, test_interactions_high

data = pd.read_csv(snakemake.input["train"])
aggregated_trees_df = pd.read_csv(snakemake.input["trees"])

synergies = {}

aggregated_trees_list = aggregated_trees_df.tree.unique()  # shouldn't this just be 1:numberOfTrees?
for tree_nr in aggregated_trees_list:
    tree_df = aggregated_trees_df[aggregated_trees_df.tree==tree_nr]  # get rows to current tree
    tree_json = from_table_to_json(tree_df)        # change format
    get_interactions(tree_json, [], synergies)   

df = pd.DataFrame({'variants':synergies.keys(), 'repetitions':synergies.values()})
df['order'] = df.variants.apply(lambda x: x.count('+')+1)

tested_synergies = test_interactions_high(df, data, max_order=snakemake.params["maxOrder"], repetitions_threshold=snakemake.params["minReps"], drop_nans=False) #here you define which order of interactions you want to compute
    
if tested_synergies:    
    tested_synergies = pd.concat(tested_synergies)
    tested_synergies['drug'] = snakemake.params["drug"]
    
    tested_mask = tested_synergies.order == tested_synergies.coef_id.str.count(":")+1
    tested_synergies = tested_synergies[tested_mask]
    
    tested_synergies.to_csv(snakemake.output[0], sep="\t", index=False)
else: 
    pd.DataFrame(columns = ["coef_id", "coef", "std err", "z", "P>|z|", "[0.025", "0.975]", "converged", "pseudo_r2", "standard_fitting", "snps", "order", "drug"]).to_csv(snakemake.output[0], sep='\t', index=False)


del tested_synergies, aggregated_trees_list, aggregated_trees_df