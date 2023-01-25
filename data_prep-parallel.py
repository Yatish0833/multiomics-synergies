#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# You may need to remove the '#' from the following commands to install the required dependencies

#these two are to read the excel
#! pip install xlrd
#! pip install install openpyxl


# In[ ]:


## ***** THESE ARE THE DEFAULTS UNLESS THEY ARE CHANGED WHEN YOU RUN THE CODE!!! *****

min_repetitions = 2 #Number of repetitions an interaction appears in the trees
max_order = 4 
working_dir = 'tmp/' #make sure it is empty

#Ranger parameters
n_trees = 1000 #Number of trees
mtry = 100 # Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables. 
max_depth = 0# Maximal tree depth. A value of NULL or 0 (the default) corresponds to unlimited depth, 1 to tree stumps (1 split per tree).
min_node= 5 # Minimal node size. default ranger: 5

# File inputs
x_input = 'data/ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx'
y_input = 'data/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv'


# In[ ]:


# Importing libraries
import pandas as pd
import numpy as np
import os
import sys
import argparse
import shutil

import statsmodels.formula.api as smf
from itertools import combinations


# In[ ]:


# Defining the number of splits in the data
try:
    parser = argparse.ArgumentParser()
    parser.add_argument('--split', default=1)
    parser.add_argument('--nSplits', default=1) #This makes it not parallel
    parser.add_argument('--minRep', default=min_repetitions)
    parser.add_argument('--maxOrder', default=max_order)
    parser.add_argument('--nTrees', default=n_trees)
    parser.add_argument('--mtry', default=mtry)
    parser.add_argument('--maxDepth', default=max_depth)
    parser.add_argument('--minNode', default=min_node)
    parser.add_argument('--workingDir', default=working_dir)
    
    
    parser_args = parser.parse_args()
    
    
    split_nr = int(parser_args.split)
    n_splits = int(parser_args.nSplits)
    
    
    min_repetitions = int(parser_args.minRep) 
    max_order = int(parser_args.maxOrder)
    working_dir = parser_args.workingDir

    #Ranger parameters
    n_trees = int(parser_args.nTrees)
    mtry = int(parser_args.mtry)
    max_depth = int(parser_args.maxDepth)
    min_node= int(parser_args.minNode)


except Exception as e: 
    print(e)
    split_nr = 1
    n_splits = 2
    print('Did you forget to add the split number and the number of splits?')
    #exit()

    


# In[ ]:


#split_nr = 1
#n_splits = 2


# In[ ]:


# If this split finalized successfully, it does not re-run
if os.path.exists(working_dir+f"final_results{split_nr}.tsv") and os.path.exists(working_dir+f"tree_performances{split_nr}.tsv"):
    print('Job previously run successfully!\nExiting')
    exit()


# In[ ]:


# Easing things for other systems
dir_trees_tmp = working_dir+"tmp" # temportal directory where each of the separate trees are going to be saved
#final_results_dir ="tmp/" # directory where the final results are going to be saved
path_to_R = 'Rscript' # Path to run R
path_to_ranger_script = 'ranger_run-parallel.R' # Path to the ranger script


# In[ ]:


#One row per cell line
#DIR
x = pd.read_excel(x_input, engine='openpyxl').drop(columns=['Project_Identifier'])
c = [a.replace(';','.') for a in x.columns]
x.columns = c
x.columns


# In[ ]:


#DIR
y = pd.read_csv(y_input)[['drug_id','cell_line_name','ln_IC50']]
y.columns


# In[ ]:


# This cell is the function to go from the table to a JSON file (variantSpark format and structure)

def merge_tree_node(tree, node):
    
    # Empty tree
    if type(tree)==float: return tree
    if len(tree)==0: return node

    # Direct children
    if tree['right'] == node['nodeID']:
        tree['right'] = node
        return tree
    elif tree['left'] == node['nodeID']:
        tree['left'] = node
        return tree

    # Create
    right = merge_tree_node(tree['right'], node)
    left = merge_tree_node(tree['left'], node)
    tree['right'] = right
    tree['left'] = left
    return tree
            

def from_table_to_json(m):
    tree = {}
    for _id,row in m.iterrows():
        current_node = {'nodeID': row['nodeID'], 
                        'splitvarID':row['splitvarID'],
                        'splitVar':row['splitvarName'],
                        'splitval':row['splitval'], 
                        'terminal':row['terminal'], 
                        'prediction':row['prediction'], 
                        'left':row['leftChild'], 
                        'right':row['rightChild'] }
        tree = merge_tree_node(tree, current_node)
    return tree



# Test
#m = pd.read_csv('output/tree1.csv')
#from_table_to_json(m)


# In[ ]:


# Algorithm to get the interacting nodes (no testing done yet)

def get_interactions(tree, current_list, interactions):
    if not 'splitVar' in tree.keys():
        return 0
    if str(tree['splitVar']) == 'nan': return 0 #ranger adds a fake predicting node at the end
    
    # Adding the interaction
    current_list.append(tree['splitVar'])
    if len(current_list) >= 2:
        for i in range(2,len(current_list)+1):
            aux = '+'.join(sorted(current_list[-i:]))
            if aux in interactions.keys():
                interactions[aux] +=1
            else:
                interactions[aux] = 1
                    
    if 'left' in tree.keys():
        get_interactions(tree['left'], current_list, interactions)
    if 'right' in tree.keys():
        get_interactions(tree['right'], current_list, interactions)
        
    _ = current_list.pop()
    


# In[ ]:


# Testing all the interactions

def results_fit_to_df(results):
    coeffs = results.params.tolist()
    pvals = results.pvalues.tolist()
    pseudo_r2 = results.rsquared
    tvals = results.tvalues.tolist()
    cint_low = results.conf_int()[0].tolist()
    cint_high = results.conf_int()[1].tolist()

    results = results.summary()
    converged = results.tables[0].data[5][1].strip()
    results = results.tables[1].data
    results = pd.DataFrame(results[1:], columns=['coef_id', 'coef', 'std err', 'z', 'P>|z|', '[0.025', '0.975]'])
    results['P>|z|'] = pvals
    results['z'] = tvals 
    results['coef'] = coeffs
    results['converged'] = converged
    results['pseudo_r2'] = pseudo_r2
    results['[0.025'] = cint_low
    results['0.975]'] = cint_high
    return results
    
def test_interactions_high(df, data, max_order=4, repetitions_threshold=2, min_cell_lines=20):
    """
    I use GLM because:
    The main difference between the two approaches is that the general linear model strictly assumes that
    the residuals will follow a conditionally normal distribution, while the GLM loosens this assumption 
    and allows for a variety of other distributions from the exponential family for the residuals.
    """
    final_results = []
    counts = 0

    for m_or in range(2,max_order+1):
        #print('current order',m_or)
        
        for v in df[(df.repetitions>=2) & (df.order==m_or)].variants.tolist():
            #preparing the input
            sp=v.split('+')
            xy = data[sp+['ln_IC50']].dropna()#.fillna(-1)
            if len(xy) <min_cell_lines: continue 
            if len(xy.columns) <3: continue
            sp=v.replace('_','').split('+')
            xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns]
            formula = xy.columns[-1]+' ~ '
            for i in range(1,len(xy.columns)):
                formula = formula + ' + '.join(['*'.join(o) for o in list(combinations(xy.columns[:-1],i))])
                formula = formula + ' + '
            formula = formula.rstrip(' + ')
            
            #Recreating the formula
            if m_or>2:
                #gathering all interactions
                fs = formula.split(' + ')
                formula = ' + '.join([a for a in fs if '*' not in a]+[a for a in fs if a.count('*')== m_or-1])
                all_interactions = [a.replace('*',':') for a in fs if '*' in a]
                final_results = pd.concat(final_results)
                subset = final_results[final_results.coef_id.apply(lambda a: a in all_interactions)].reset_index(drop=True)
                final_results = [final_results]
                if len(subset)>0:
                    max_idx = subset['coef'].astype(float).abs().idxmax()
                    coef_id = subset.loc[max_idx].coef_id
                    formula = formula +' + '+coef_id.replace(':','*')
                else:
                    #pass
                    continue # bc i dont think it is a valid tree form (interaction-wise)
                    #There is no sub epistasis (P>Q>O>P, tree 503, first compound)

            # Standard fitting
            try:
                ols = smf.ols(formula.replace('*',':'),data=xy)
                # "*" vs ":" #https://stackoverflow.com/questions/33050104/difference-between-the-interaction-and-term-for-formulas-in-statsmodels-ols
            except:
                print('error in OLS')
                print('coef_id',coef_id)
                print('formula OLS',type(formula),formula)
                #return pd.concat(final_results)
                continue
            ols.raise_on_perfect_prediction = False #preventing the perfect separation error
            results = ols.fit(disp=False, maxiter=1000) #mehtod prevents singular matrix
#            return results
            results = results_fit_to_df(results)
            results['standard_fitting'] = True

            #If nan means no convergence bc singular matrix
            #adding regularization
            if 'nan' == pd.DataFrame(results)['z'].astype(str).iloc[2].strip():
                try:
                    results = ols.fit_regularized(method='l1', disp=False, maxiter=1000, alpha=0.3) #mehtod prevents singular matrix
                    results = results_fit_to_df(results)
                    results['standard_fitting'] = False        
                except:
                    #crashed the regularized
                    counts +=1
                    continue


            results['snps'] = v
            results['order'] = len(sp)
            final_results.append(results)

    final_results = pd.concat(final_results)
    return final_results


# In[ ]:


def undo(string):
    
    string = ''.join([ x if ord(x)<90 else str(ord(x)-97) for x in string ])
    string = string[:6]+'.'+string[6:].replace('HUMAN', '_HUMAN') #not sure these 6
    return string
    
#undo('PacfbbCRYABHUMAN')


# In[ ]:


#TODO: By compound

#Looping over all drugs
# drug_id
# Other options:
# - drug_name
# - CHEMBL = Chemical compound ID
#for compound_name, group in x.merge(y, left_on='Cell_Line', right_on='cell_line_name').groupby('drug_id'): # may require too much memory

#Making a temp file to run all R stuff
#DIR
os.system("mkdir -p "+dir_trees_tmp+f"{split_nr}")


column_to_group = 'drug_id'
drugs_list = y[column_to_group].sort_values().unique()
drugs_list = [drugs_list[i] for i in range(len(drugs_list)) if i%n_splits==split_nr-1]
i = -1
all_drug_results = []
tree_performances = []
for elm in drugs_list:
    print(f"{path_to_R} {path_to_ranger_script}  -w {working_dir} -c {split_nr} -n {n_trees} -t {mtry} -s {min_node} -d {max_depth}")
    i+=1
    
    if i%10==0 or i<10: print(i,'out of',len(drugs_list), 'drugs in split nr', split_nr)

    xy = x.merge(y[y[column_to_group]==elm], left_on='Cell_Line', right_on='cell_line_name')
    #Enhancement: Remove peptides that are all zero 
    
    # saving csv for R df
    # file name is generic but we could personalize it
    #DIR
    xy.drop(columns=['Cell_Line', 'cell_line_name','drug_id']).fillna(0).rename(columns={'ln_IC50':'label'}).to_csv(dir_trees_tmp+f"{split_nr}/data.csv", index=False)

    #Run the R script to generate the outputs
    #os.system(f"{path_to_R} {path_to_ranger_script}  -w {working_dir} -c {split_nr} -n {n_trees} -t {mtry} -s {min_node} -d {max_depth}")
    print("R output (in case there is an error or something")
    #os.popen(f"{path_to_R} {path_to_ranger_script}  -w {working_dir} -c {split_nr} -n {n_trees} -t {mtry} -s {min_node} -d {max_depth}").read()    
    print(os.popen(f"{path_to_R} {path_to_ranger_script}  -w {working_dir} -c {split_nr} -n {n_trees} -t {mtry} -s {min_node} -d {max_depth}").read())
    
    #load the R outputs (the trees, one file each), and convert it to VS look-alike and get interactions
    interactions = {}
    #DIR
#    trees = os.listdir(dir_trees_tmp+f"{split_nr}")
    #files = [x for x in files if 'tree' in x]
#    for tree in trees:
#        if 'tree' not in tree: continue #if it is not a tree file ignore
        #DIR
#        tree_json = from_table_to_json(pd.read_csv(os.path.join(dir_trees_tmp+str(split_nr),tree)))        
#        get_interactions(tree_json,[],interactions) #the interactions are found in "interactions"
        
    aggregated_trees_df = pd.read_csv(os.path.join(dir_trees_tmp+str(split_nr),'aggregated_trees.csv'))
    aggregated_trees_list = aggregated_trees_df.tree.unique()
    for tree_nr in aggregated_trees_list:
        tree_df = aggregated_trees_df[aggregated_trees_df.tree==tree_nr]
        tree_json = from_table_to_json(tree_df)        
        get_interactions(tree_json,[],interactions)
    
    # Creating a df out of the interactions
    df = pd.DataFrame({'variants':interactions.keys(),'repetitions':interactions.values()})
    df['order'] = df.variants.apply(lambda x: x.count('+')+1)
    
    #get tree performances
    aux_performances = pd.read_csv(os.path.join(dir_trees_tmp+str(split_nr),"performance.tsv"), sep='\t')
    aux_performances['drug'] = elm
    tree_performances.append(aux_performances)
    
    tested_interactions = test_interactions_high(df, xy, max_order=max_order, repetitions_threshold=min_repetitions) #here you define which order of interactions you want to compute
    tested_interactions['drug'] = elm
    all_drug_results.append(tested_interactions)
    #break
    

    

final_results = pd.concat(all_drug_results)
final_results.to_csv(working_dir+f"final_results{split_nr}.tsv", index=False, sep='\t')

tree_performances = pd.concat(tree_performances)
tree_performances.to_csv(working_dir+f"tree_performances{split_nr}.tsv", index=False, sep='\t')

# deliting temp folder from raneger
shutil.rmtree(dir_trees_tmp+f"{split_nr}")
#Option B: os.system("rm -rf _path_to_dir")


# In[ ]:


print('Split ',split_nr, 'out of ',n_splits,' is DONE')


# In[ ]:


fr = [x for x in os.listdir(working_dir) if 'final_results' in x and 'final_results_all.tsv' not in x]
if len(fr) == n_splits:
    df = pd.concat([pd.read_csv(os.path.join(working_dir,final_result), sep='\t') for final_result in fr])
    #df = df[df.coef_id.apply(lambda x: x.count(':'))==df.snps.apply(lambda x: x.count('+'))] #this keeps only the interactions
    df.to_csv(working_dir+"final_results_all.tsv", index=False, sep='\t')
    
    #Removing temp final results
    for final_result in fr:
        #os.remove(os.path.join(working_dir,final_result))
        pass
        
    # now the same for performances    
    pr = [x for x in os.listdir(working_dir) if 'tree_performances' in x and 'tree_performances_all.tsv' not in x]
    df = pd.concat([pd.read_csv(os.path.join(working_dir,tree_performance), sep='\t') for tree_performance in pr])
    df.to_csv(working_dir+"tree_performances_all.tsv", index=False, sep='\t')
    
    #Removing temp final results
    for final_result in pr:
        #os.remove(os.path.join(working_dir,final_result))
        pass

    print('All jobs finished successfully!\n final_results_all.tsv has all the aggregated output')

