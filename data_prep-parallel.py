#!/usr/bin/env python
# coding: utf-8

# In[1]:


# You may need to remove the '#' from the following commands to install the required dependencies

#these two are to read the excel
#! pip install xlrd
#! pip install install openpyxl


# In[2]:


# Importing libraries
import pandas as pd
import numpy as np
import os
import sys

import statsmodels.formula.api as smf
from itertools import combinations


# In[3]:


min_repetitions=2


# In[4]:


# Defining the number of splits in the data
try:
    split_nr = int(sys.argv[1])
    n_splits = int(sys.argv[2])
except:
    split_nr = 1
    n_splits = 10
    print('Did you forget to add the split number and the number of splits?')
    #exit()


# In[5]:


#One row per cell line
#DIR
x = pd.read_excel('data/ProCan-DepMapSanger_protein_matrix_6692_averaged.xlsx', engine='openpyxl').fillna(0).drop(columns=['Project_Identifier'])
c = [a.replace(';','.') for a in x.columns]
x.columns = c
x.columns


# In[6]:


#DIR
y = pd.read_csv('data/DrugResponse_PANCANCER_GDSC1_GDSC2_20200602.csv')[['drug_id','cell_line_name','ln_IC50']]
y.columns


# In[7]:


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


# In[8]:


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
    


# In[9]:


# Testing all the interactions

def test_interactions(df, data):
    """
    I use GLM because:
    The main difference between the two approaches is that the general linear model strictly assumes that
    the residuals will follow a conditionally normal distribution, while the GLM loosens this assumption 
    and allows for a variety of other distributions from the exponential family for the residuals.
    """
    final_results = []
    counts = 0

    for v in df[(df.repetitions>=2) & (df.order ==2)].variants.tolist():
        
        #preparing the input
        sp=v.split('+')
        xy = data[sp+['ln_IC50']].fillna(-1)
        sp=v.replace('_','').split('+')
        xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns]
        formula = xy.columns[-1]+' ~ '
        for i in range(1,len(xy.columns)):
            formula = formula + ' + '.join(['*'.join(o) for o in list(combinations(xy.columns[:-1],i))])
            formula = formula + ' + '
        formula = formula.rstrip(' + ')

        # Standard fitting
        ols = smf.ols(formula,data=xy)
        ols.raise_on_perfect_prediction = False #preventing the perfect separation error
        results = ols.fit(disp=False, maxiter=1000) #mehtod prevents singular matrix
        results = results.summary()
        converged = results.tables[0].data[5][1].strip()
        pseudo_r2 = results.tables[0].data[3][3].strip()
        results = results.tables[1].data
        results = pd.DataFrame(results[1:], columns=['coef_id', 'coef', 'std err', 'z', 'P>|z|', '[0.025', '0.975]'])
        results['standard_fitting'] = True

        #If nan means no convergence bc singular matrix
        #adding regularization
        if 'nan' == pd.DataFrame(results)['z'].iloc[2].strip():
            try:
                results = ols.fit_regularized(method='l1', disp=False, maxiter=1000, alpha=0.3) #mehtod prevents singular matrix
                results = results.summary()
                converged = results.tables[0].data[5][1].strip()
                pseudo_r2 = results.tables[0].data[3][3].strip()
                results = results.tables[1].data
                results = pd.DataFrame(results[1:], columns=['coef_id', 'coef', 'std err', 'z', 'P>|z|', '[0.025', '0.975]'])
                results['standard_fitting'] = False        
            except:
                #crashed the regularized
                counts +=1

        results['converged'] = converged
        results['pseudo_r2'] = pseudo_r2
        results['snps'] = v
        results['order'] = len(sp)
        final_results.append(results)

    final_results = pd.concat(final_results)
    #print(counts)
    return final_results


# In[10]:


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
    
def test_interactions_high(df, data, max_order=4):
    """
    I use GLM because:
    The main difference between the two approaches is that the general linear model strictly assumes that
    the residuals will follow a conditionally normal distribution, while the GLM loosens this assumption 
    and allows for a variety of other distributions from the exponential family for the residuals.
    """
    final_results = []
    counts = 0

    for m_or in range(2,max_order+1):
        print('current order',m_or)
        
        for v in df[(df.repetitions>=2) & (df.order==m_or)].variants.tolist():
            #preparing the input
            sp=v.split('+')
            xy = data[sp+['ln_IC50']].fillna(-1)
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
                    max_idx = subset['coef'].astype(float).idxmax()
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


# In[11]:


def undo(string):
    
    string = ''.join([ x if ord(x)<90 else str(ord(x)-97) for x in string ])
    string = string[:6]+'.'+string[6:].replace('HUMAN', '_HUMAN') #not sure these 6
    return string
    
#undo('PacfbbCRYABHUMAN')


# In[12]:


#TODO: By compound

#Looping over all drugs
# drug_id
# Other options:
# - drug_name
# - CHEMBL = Chemical compound ID
#for compound_name, group in x.merge(y, left_on='Cell_Line', right_on='cell_line_name').groupby('drug_id'): # may require too much memory

#Making a temp file to run all R stuff
#DIR
os.system(f"mkdir -p tmp/tmp{split_nr}")


column_to_group = 'drug_id'
drugs_list = y[column_to_group].sort_values().unique()
drugs_list = [drugs_list[i] for i in range(len(drugs_list)) if i%n_splits==split_nr]
i = -1
all_drug_results = []
for elm in drugs_list[:3]:
    i+=1
    
    if i%10==0 or i<10: print(i,'out of',len(drugs_list), 'drugs in split nr', split_nr)

    xy = x.merge(y[y[column_to_group]==elm], left_on='Cell_Line', right_on='cell_line_name')
    #Enhancement: Remove peptides that are all zero 
    
    # saving csv for R df
    # file name is generic but we could personalize it
    #DIR
    xy.drop(columns=['Cell_Line', 'cell_line_name','drug_id']).rename(columns={'ln_IC50':'label'}).to_csv(f"tmp/tmp{split_nr}/data.csv", index=False)

    #Run the R script to generate the outputs
    os.system(f"Rscript ranger_run.R {split_nr}")
    
    #load the R outputs (the trees, one file each), and convert it to VS look-alike and get interactions
    interactions = {}
    #DIR
    trees = os.listdir(f"tmp/tmp{split_nr}/")
    #files = [x for x in files if 'tree' in x]
    for tree in trees:
        if 'tree' not in tree: continue #if it is not a tree file ignore
        #DIR
        tree_json = from_table_to_json(pd.read_csv(f"tmp/tmp{split_nr}/"+tree))        
        get_interactions(tree_json,[],interactions) #the interactions are found in "interactions"
        
    # Creating a df out of the interactions
    df = pd.DataFrame({'variants':interactions.keys(),'repetitions':interactions.values()})
    df['order'] = df.variants.apply(lambda x: x.count('+')+1)
    
    
    tested_interactions = test_interactions_high(df, xy) #here you define which order of interactions you want to compute
    tested_interactions['drug'] = elm
    all_drug_results.append(tested_interactions)
    

    

final_results = pd.concat(all_drug_results)
#DIR
final_results.to_csv(f"tmp/final_results{split_nr}.tsv", index=False, sep='\t')


# In[13]:


print('Done:',split_nr)

