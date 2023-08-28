import pandas as pd
import numpy as np
import os
import sys
import argparse
import shutil

import statsmodels.formula.api as smf
from itertools import combinations

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


def get_interactions(tree, current_list, interactions):
    if not 'splitVar' in tree.keys():
        return 0
    if str(tree['splitVar']) == 'nan': return 0 #ranger adds a fake predicting node at the end
    
    # Adding the interaction
    current_list.append(tree['splitVar'])
    if len(current_list) >= 2:
        for i in range(2,len(current_list)+1):
              if len(current_list[-i-1:]) == len(set(current_list[-i-1:])) and "maxscreeningconc" not in current_list:
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


def undo(string):    
    string = ''.join([ x if ord(x)<91 else str(ord(x)-97) for x in string ])
    string = string[:6]+'.'+string[6:].replace('HUMAN', '_HUMAN') #not sure these 6
    return string

def results_fit_to_df(results):
    coeffs = results.params.tolist()
    pvals = results.pvalues.tolist()
    pseudo_r2 = results.rsquared
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
    results['[0.025'] = cint_low
    results['0.975]'] = cint_high
    return results
    
def test_interactions_high(df, data, max_order=4, repetitions_threshold=2, min_samples=20, drop_nans=False):
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
            if len(sp) != len(set(sp)): continue
            xy = data[sp+["maxscreeningconc", 'label']]
            if drop_nans:
                xy = xy.loc[~(df==0).all(axis=1)]
            if len(xy.columns) <= m_or: continue #not enough columns # how would we get that case? deccision over the same protein several times?
            if len(xy) <= min_samples: continue #not enough rows
            sp=v.replace('_','').split('+')
            xy.columns = [''.join([chr(int(y)+97) if y.isnumeric() else y for y in x.replace('_','').replace('.','')]) for x in xy.columns]
            formula = xy.columns[-1]+' ~ '
            for i in range(1,len(xy.columns)):
                formula = formula + ' + '.join(['*'.join(o) for o in list(combinations(xy.columns[:-2],i))])
                formula = formula + ' + '
            formula = formula.rstrip(' + ')
            # print("We are using this formula for this case")
            # print(formula)
            
            #Recreating the formula
            if m_or>2:
                #gathering all interactions
                fs = formula.split(' + ')
                formula = ' + '.join([a for a in fs if '*' not in a]+[a for a in fs if a.count('*')== m_or-1])
                all_interactions = [a.replace('*',':') for a in fs if '*' in a]
                final_results = pd.concat(final_results)
                subset = final_results[final_results.coef_id.apply(lambda a: a in all_interactions)].reset_index(drop=True)
                # print("Before formula for this case")
                # print(formula)
                # print("this the subset: ", subset)
                final_results = [final_results]
                if len(subset)>0:
                    max_idx = subset['coef'].astype(float).abs().idxmax()
                    coef_id = subset.loc[max_idx].coef_id
                    formula = formula +' + '+coef_id.replace(':','*')
                    # print("Updated formula for this case")
                    # print(formula)
                else:
                    #pass
                    continue # bc i dont think it is a valid tree form (interaction-wise)
                    #There is no sub epistasis (P>Q>O>P, tree 503, first compound), P>Q>O and Q>O>P not counted as 2 repetitions
                    
            
            
            # Standard fitting
            try:
                ols = smf.ols(formula.replace('*',':') + " + maxscreeningconc",data=xy)  # why are we not encoding it the right way immediately?
                # "*" vs ":" #https://stackoverflow.com/questions/33050104/difference-between-the-interaction-and-term-for-formulas-in-statsmodels-ols
            except:
                print('error in OLS')
                #print('coef_id',coef_id)
                print('formula OLS',type(formula),formula)
                #return pd.concat(final_results)
                continue
            ols.raise_on_perfect_prediction = False #preventing the perfect separation error
            results = ols.fit(disp=False, maxiter=1000) #method prevents singular matrix
#            print(formula)
#            return results
            results = results_fit_to_df(results)
            results['standard_fitting'] = True

            #If nan means no convergence bc singular matrix
            #adding regularization
            if results['z'].isna().sum():
                print("found missing value")
                try:
                    results = ols.fit_regularized(method='l1', disp=False, maxiter=1000, alpha=0.3) #mehtod prevents singular matrix
                    results = results_fit_to_df(results)
                    results['standard_fitting'] = False        
                except:
                    #crashed the regularized
                    counts +=1
                    continue


            results['snps'] = '*'.join([undo(x) for x in v.split('+')])   # aren't these proteins here
            results['order'] = len(sp)
            final_results.append(results)

    return final_results