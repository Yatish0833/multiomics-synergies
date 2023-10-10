import pandas as pd
import numpy as np

import os


for synergie in snakemake.input["synergies"]:
    if os.path.exists(snakemake.output["synergies"]):
        pd.read_csv(synergie, sep='\t').to_csv(snakemake.output["synergies"], index=False, sep='\t', mode='a', header=False)
    else:
        pd.read_csv(synergie, sep='\t').to_csv(snakemake.output["synergies"], index=False, sep='\t',mode='a')
        

for tree_score in snakemake.input["rf_score"]:
    if os.path.exists(snakemake.output["rf_scores"]):
        pd.read_csv(tree_score, sep='\t').to_csv(snakemake.output["rf_scores"], index=False, sep='\t', mode='a', header=False)
    else:
        pd.read_csv(tree_score, sep='\t').to_csv(snakemake.output["rf_scores"], index=False, sep='\t',mode='a')
        