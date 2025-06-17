from collections import defaultdict
from csv import DictReader
from typing import Dict, List, Tuple
import json
import os
import itertools
import argparse

import numpy as np
import pandas as pd

BASEDIR = os.getcwd()
CSV = BASEDIR + "/data/ProCan-DepMapSanger_protein_matrix_6692_averaged.csv"
SAMPLE_COLUMN = 'Cell_Line'
PROTEIN_ID_COLUMN = 'protein_id'
PHENOTYPE_COLUMN = 'ln_IC50'

MAX_ORDER = 3
MAX_TERMS = 5

CORRELATION_THRESHOLD = 0.7

seed = 182652386547354083 + np.random.randint(0, 1000000)
rng = np.random.default_rng(seed=seed)

def sigmoid(x: float) -> float:
    return 1 / (1 + np.exp(-x))

def get_marginal_effect() -> float:
    """Randomly sample a marginal effect from a uniform distribution between -0.9 and 0.9"""
    return np.random.choice([-1, 1])

def get_interaction_effect(marginal_effects: List[float]) -> float:
    """
    Randomly sample an interaction effect from a normal distribution,
    ensuring it's smaller than the smallest absolute marginal effect.
    """
    min_abs_marginal = min(abs(effect) for effect in marginal_effects)
    scale = min_abs_marginal / 3  # Ensures most effects are smaller
    return np.random.normal(scale=scale)

def get_clean_dataframe(csv_path: str) -> pd.DataFrame:
    """Load and clean the protein data."""
    df = pd.read_csv(csv_path, index_col=0)
    df = df.fillna(df.mean())
    with open(csv_path, encoding="utf-8") as in_csv:
        reader = DictReader(in_csv)
        columns = reader.fieldnames
    seen_columns = set()
    duplicate_columns = defaultdict(list)
    for df_column, csv_column in zip(df.columns, columns):
        if csv_column in seen_columns:
            duplicate_columns[csv_column].append(df_column)
        else:
            seen_columns.add(csv_column)
    for column, duplicates in duplicate_columns.items():
        if column == '':
            continue
        if duplicates:
            df[column] = df[[column] + duplicates].mean(axis='columns')
    return df.drop(
        labels=[
            dup
            for duplicate_list in duplicate_columns.values()
            for dup in duplicate_list
        ],
        axis='columns'
    ).dropna(axis='columns')

def convert_ndarrays_to_lists(obj):
    """Recursively convert Numpy arrays and scalars in a dict or list to Python native types."""
    if isinstance(obj, dict):
        return {key: convert_ndarrays_to_lists(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_ndarrays_to_lists(element) for element in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    else:
        return obj

def identify_correlated_proteins(correlation_matrix: pd.DataFrame,
                                 threshold: float) -> Dict[str, Dict[str, float]]:
    """Identify correlated proteins using a pre-computed correlation matrix."""
    correlated_proteins = {}
    for protein in correlation_matrix.columns:
        correlated = correlation_matrix[protein][
            (correlation_matrix[protein].abs() > threshold) & 
            (correlation_matrix[protein] < 1.0)
        ]
        if not correlated.empty:
            correlated_proteins[protein] = {
                f"{partner}:{coef:.2f}": coef
                for partner, coef in correlated.items()
            }
    return correlated_proteins

def add_protein_effects(protein: str, effect: float, protein_data: pd.DataFrame) -> pd.Series:
    """ individual 'main' effects """
    effect_value = effect * protein_data[protein]
    return effect_value

def add_noise_to_phenotype(phenotype: pd.Series, percent_explainability: float) -> pd.Series:
    """Add noise to the phenotype based on percent_explainability."""
    std_dev = phenotype.std()
    if percent_explainability > 0:
        error_std = std_dev * np.sqrt(100/percent_explainability - 1)
        noisy_phenotype = phenotype + rng.normal(0, error_std, len(phenotype))
    else:
        noisy_phenotype = rng.normal(0, std_dev, len(phenotype))
    return noisy_phenotype

def get_uncorrelated_proteins(protein_data: pd.DataFrame, 
                              correlated_proteins: Dict[str, Dict[str, float]], 
                              num_proteins: int) -> List[str]:
    """Select proteins that have no correlated partners."""
    uncorrelated_proteins = [p for p in protein_data.columns if p not in correlated_proteins]
    
    if len(uncorrelated_proteins) < num_proteins:
        print(f"Warning: Not enough uncorrelated proteins. Using {len(uncorrelated_proteins)} proteins.")
        return uncorrelated_proteins
    
    return rng.choice(uncorrelated_proteins, size=num_proteins, replace=False).tolist()

def generate_phenotype(protein_data: pd.DataFrame, model: str, proteins: List[str], 
                       correlated_proteins: Dict[str, Dict[str, float]],
                       base_effect: float, marginal_effects: List[float]) -> Tuple[pd.Series, Dict]:
    """Generate continuous phenotype based on the specified interaction model"""
    phenotype = pd.Series(base_effect, index=protein_data.index)

    actual_effects = {
        'model': model,
        'base_effect': base_effect,
        'marginal_effects': {},
        'interaction_effect': {},
        'proteins': proteins,
        'correlated_proteins': {p: list(correlated_proteins.get(p, {}).keys()) for p in proteins}
    }

    protein_effects = [add_protein_effects(p, effect, protein_data)
                       for p, effect in zip(proteins, marginal_effects)]

    if model == 'Interaction_only':
        interaction_matrix = np.array([protein_data[p] for p in proteins])
        
        if len(proteins) == 2:
            pairwise_interactions = interaction_matrix[0] * interaction_matrix[1]
            interaction_effects = np.array([get_interaction_effect(marginal_effects)])
            interaction_contribution = interaction_effects[0] * pairwise_interactions
            actual_effects['interaction_effect'][f'{proteins[0]}_{proteins[1]}'] = interaction_effects[0]
        else:
            pairwise_interactions = np.einsum('i...,j...->ij...', interaction_matrix, interaction_matrix)
            i_upper, j_upper = np.triu_indices_from(pairwise_interactions[0], k=1)
            valid_indices = (i_upper < pairwise_interactions.shape[0]) & (j_upper < pairwise_interactions.shape[1])
            i_upper = i_upper[valid_indices]
            j_upper = j_upper[valid_indices]
            pairwise_interactions = pairwise_interactions[i_upper, j_upper]
            interaction_effects = np.array([get_interaction_effect(marginal_effects) for _ in range(len(pairwise_interactions))])
            interaction_contribution = np.sum(interaction_effects[:, np.newaxis] * pairwise_interactions, axis=0)
            for idx, (i, j) in enumerate(zip(i_upper, j_upper)):
                actual_effects['interaction_effect'][f'{proteins[i]}_{proteins[j]}'] = interaction_effects[idx]
        
        if len(proteins) == 3:
            three_way_interaction = np.prod(interaction_matrix, axis=0)
            three_way_interaction_effect = get_interaction_effect(marginal_effects)
            interaction_contribution += three_way_interaction_effect * three_way_interaction
            interaction_effects = np.append(interaction_effects, three_way_interaction_effect)
            actual_effects['interaction_effect']['three_way'] = three_way_interaction_effect
        
        phenotype += interaction_contribution

    elif model == 'Modifier_protein':
        phenotype += protein_effects[0]
        actual_effects['marginal_effects'][proteins[0]] = marginal_effects[0]

        if len(proteins) == 2:
            interaction_effect = get_interaction_effect(marginal_effects)
            two_way_interaction = protein_data[proteins[0]] * protein_data[proteins[1]]
            phenotype += interaction_effect * two_way_interaction
            actual_effects['interaction_effect'][f'{proteins[0]}_{proteins[1]}'] = interaction_effect
        else:
            interaction_matrix = np.array([protein_data[p] for p in proteins[1:]])
            two_way_interactions = protein_data[proteins[0]].values[:, np.newaxis] * interaction_matrix.T

            num_interactions = len(proteins) - 1
            interaction_effects = np.array([get_interaction_effect(marginal_effects) for _ in range(num_interactions)])
            
            interaction_effects = interaction_effects[:, np.newaxis].T

            phenotype += pd.Series(np.sum(interaction_effects * two_way_interactions, axis=1), index=phenotype.index)
            
            for i, protein in enumerate(proteins[1:]):
                if i <= len(interaction_effects):
                    actual_effects['interaction_effect'][f'{proteins[0]}_{protein}'] = interaction_effects[0, i]

            if len(proteins) >= 3:
                three_way_interaction = np.prod([protein_data[p] for p in proteins], axis=0)
                three_way_interaction_effect = get_interaction_effect(marginal_effects)
                phenotype += three_way_interaction_effect * three_way_interaction
                actual_effects['interaction_effect']['three_way'] = three_way_interaction_effect

    elif model == 'No_Interaction':
        for protein, effect, protein_effect in zip(proteins, marginal_effects, protein_effects):
            phenotype += protein_effect
            actual_effects['marginal_effects'][protein] = effect

    elif model == 'Dominant':
        phenotype += np.max(protein_effects, axis=0)
        for protein, effect in zip(proteins, marginal_effects):
            actual_effects['marginal_effects'][protein] = effect

    elif model == 'Redundant':
        combined_effect = np.sum(protein_effects, axis=0)
        interaction_effect = -np.min(protein_effects, axis=0)
        phenotype += combined_effect + interaction_effect

        for protein, effect in zip(proteins, marginal_effects):
            actual_effects['marginal_effects'][protein] = effect
        min_effect = np.min(marginal_effects)
        min_protein = proteins[np.argmin(marginal_effects)]
        actual_effects['interaction_effect'] = {
            'redundancy': min_effect,
            'protein': min_protein
        }

    elif model in ['Synergistic', 'Cubic', 'Sigmoidal', 'Exponential']:
        interaction_effect = get_interaction_effect(marginal_effects)

        for protein, effect, protein_effect in zip(proteins, marginal_effects, protein_effects):
            phenotype += protein_effect
            actual_effects['marginal_effects'][protein] = effect
        
        if model == 'Synergistic':
            interaction_term = np.prod([protein_data[p] for p in proteins], axis=0)
            interaction_effect = np.random.uniform(1.1, 1.5) * interaction_effect
        
        elif model == 'Cubic':
            interaction_term = np.sum([protein_data[p] for p in proteins], axis=0) ** 3

        elif model == 'Sigmoidal':
            interaction_term = sigmoid(np.sum([protein_data[p] for p in proteins], axis=0))

        else:  # Exponential
            interaction_term = np.exp(np.sum([protein_data[p] for p in proteins], axis=0))
        
        phenotype += interaction_effect * interaction_term
        actual_effects['interaction_effect']['non_linear'] = interaction_effect

    return phenotype, actual_effects

def generate_complex_phenotype(protein_data: pd.DataFrame, 
                               correlated_proteins: Dict[str, Dict[str, float]],
                               correlation_matrix: pd.DataFrame,
                               max_protein_sets: int,
                               max_order: int) -> Tuple[pd.Series, List[Dict]]:
    """Generate a complex phenotype based on multiple sets of interacting proteins."""
    phenotype = pd.Series(0, index=protein_data.index)
    protein_set_info = []

    num_protein_sets = rng.integers(1, max_protein_sets + 1)

    for _ in range(num_protein_sets):
        order = rng.integers(2, max_order + 1)
        proteins = get_uncorrelated_proteins(protein_data, correlated_proteins, order)
        model = rng.choice(['Interaction_only', 'Modifier_protein', 'Dominant', 'No_Interaction', 'Redundant', 'Synergistic',
                            'Cubic', 'Sigmoidal', 'Exponential'])
        base_effect = 0
        marginal_effects = [get_marginal_effect() for _ in range(len(proteins))]

        set_phenotype, actual_effects = generate_phenotype(
            protein_data, model, proteins, correlated_proteins, 
            base_effect, marginal_effects)
        
        phenotype += set_phenotype
        protein_set_info.append(actual_effects)
    
    return phenotype, protein_set_info

def generate_polynomial_phenotype(protein_data: pd.DataFrame, 
                                  correlated_proteins: Dict[str, Dict[str, float]],
                                  max_order: int,
                                  max_terms: int) -> Tuple[pd.Series, Dict]:
    phenotype = pd.Series(0, index=protein_data.index)
    coefficients = {}
    
    proteins = get_uncorrelated_proteins(protein_data, correlated_proteins, max_order)
    
    for order in range(1, max_order + 1):
        combinations = list(itertools.combinations_with_replacement(proteins, order))
        
        if len(combinations) > max_terms:
            combinations = rng.choice(combinations, size=max_terms, replace=False)
        
        for combo in combinations:
            if order == 1:
                coef = get_marginal_effect()
            else:
                coef = get_interaction_effect([get_marginal_effect() for _ in range(order)])
            term = protein_data[list(combo)].prod(axis=1)
            phenotype += coef * term
            coefficients['+'.join(combo)] = coef
    
    actual_effects = {
        'model': 'Polynomial',
        'proteins': proteins,
        'coefficients': coefficients
    }

    return phenotype, actual_effects

def write_csv_header(file) -> None:
    header = [SAMPLE_COLUMN, PHENOTYPE_COLUMN]
    file.write(','.join(header) + '\n')

def write_phenotype_to_csv(file, sample, phenotype) ->  None:
    row = [str(sample), str(phenotype)]
    file.write(','.join(row) + '\n')

def main() -> None:
    parser = argparse.ArgumentParser(description="Simulate phenotypes with configurable parameters.")
    parser.add_argument("--percent_explainability", type=int, default=70, help="Percentage of phenotype variance explained (default: 70)")
    parser.add_argument("--version", type=str, default="V1", help="Version string for output directory (default: V1)")
    args = parser.parse_args()

    PERCENT_EXPLAINABILITY = args.percent_explainability
    VERSION = args.version

    OUTPUT_DIR = BASEDIR + "/results/simulated-phenotypes-her-" + str(PERCENT_EXPLAINABILITY) + "_" + VERSION

    print("Loading and cleaning protein data...")
    clean_protein_data = get_clean_dataframe(CSV)
    numeric_data = clean_protein_data.select_dtypes(include=[np.number])
    protein_data = numeric_data.loc[:, (numeric_data != numeric_data.iloc[0]).any()]
    
    protein_data = protein_data - protein_data.mean()
    
    correlation_matrix_file = os.path.join(BASEDIR, "results", "correlation_matrix.csv")
    if os.path.exists(correlation_matrix_file):
        print(f"Loading existing correlation matrix from {correlation_matrix_file}")
        correlation_matrix = pd.read_csv(correlation_matrix_file, index_col=0)
    else:
        print("Calculating new correlation matrix as it doesn't exist...")
        correlation_matrix = protein_data.corr()
        correlation_matrix.to_csv(correlation_matrix_file)
        print(f"Correlation matrix saved to {correlation_matrix_file}")

    print("Identifying correlated proteins...")
    correlated_proteins = identify_correlated_proteins(correlation_matrix, CORRELATION_THRESHOLD)
    
    interacting_proteins = []
    

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("Generating phenotypes and saving to CSV files...")
    
    models = ['Interaction_only', 'Modifier_protein', 'No_Interaction', 'Dominant', 'Redundant', 'Synergistic',
              'Cubic', 'Sigmoidal', 'Exponential']
    total_iterations = len(models) * (MAX_ORDER - 1)
    current_iteration = 0
    
    for model in models:
        for order in range(2, MAX_ORDER + 1):
            current_iteration += 1
            print(f"Processing model: {model}, order: {order} ({current_iteration}/{total_iterations})")
            
            output_file = os.path.join(OUTPUT_DIR, f"{model}_order_{order}_phenotypes.csv")
            with open(output_file, 'w', encoding="utf-8") as f:
                write_csv_header(f)
                
                proteins = get_uncorrelated_proteins(protein_data, correlated_proteins, order)
                base_effect = 0
                marginal_effects = [get_marginal_effect() for _ in range(len(proteins))]
                
                phenotype, actual_effects = generate_phenotype(
                    protein_data, model, proteins, correlated_proteins, 
                    base_effect, marginal_effects
                )
                
                # Add noise to the complex phenotype based on proportion of pheno variance explained
                noisy_phenotype = add_noise_to_phenotype(phenotype, percent_explainability=PERCENT_EXPLAINABILITY)
                    
                for cell_line, value in zip(protein_data.index, noisy_phenotype):
                    write_phenotype_to_csv(f, cell_line, value)
            
            interacting_proteins.append(actual_effects)
    
    complex_output_file = os.path.join(OUTPUT_DIR, "complex_phenotypes.csv")
    complex_interacting_proteins = []
    with open(complex_output_file, 'w', encoding="utf-8") as f:
        write_csv_header(f)
        phenotype, info = generate_complex_phenotype(
                protein_data, correlated_proteins, 
                correlation_matrix, MAX_TERMS, MAX_ORDER
            )
            
        # Add noise to the complex phenotype based on proportion of pheno variance explained
        noisy_phenotype = add_noise_to_phenotype(phenotype, PERCENT_EXPLAINABILITY)
            
        for cell_line, value in zip(protein_data.index, noisy_phenotype):
            write_phenotype_to_csv(f, cell_line, value)

        complex_interacting_proteins.append(info)
    
    print("Saving interacting proteins info to JSON...")
    interacting_proteins = [convert_ndarrays_to_lists(d) for d in interacting_proteins]
    with open(os.path.join(OUTPUT_DIR, 'interacting_proteins.json'), 'w', encoding="utf-8") as f:
        json.dump(interacting_proteins, f, indent=2)
    
    print("Saving complex interacting proteins info to JSON...")
    complex_interacting_proteins = [convert_ndarrays_to_lists(d) for d in complex_interacting_proteins]
    with open(os.path.join(OUTPUT_DIR, 'complex_interacting_proteins.json'), 'w', encoding="utf-8") as f:
        json.dump(complex_interacting_proteins, f, indent=2)

    # generating phenotypes with polynomial relationships based on order x
    print("Generating polynomial phenotypes...")
    polynomial_output_file = os.path.join(OUTPUT_DIR, "polynomial_phenotypes.csv")
    polynomial_interacting_proteins = []
    with open(polynomial_output_file, 'w', encoding="utf-8") as f:
        write_csv_header(f)

        phenotype, info = generate_polynomial_phenotype(protein_data, correlated_proteins, MAX_ORDER, MAX_TERMS)
        # Add noise to the polynomial phenotype
        noisy_phenotype = add_noise_to_phenotype(phenotype, PERCENT_EXPLAINABILITY)
        for cell_line, value in zip(protein_data.index, noisy_phenotype):
            write_phenotype_to_csv(f, cell_line, value)
        polynomial_interacting_proteins.append(info)
    print("Saving polynomial interacting proteins info to JSON...")
    polynomial_interacting_proteins = convert_ndarrays_to_lists(polynomial_interacting_proteins)
    with open(os.path.join(OUTPUT_DIR, 'polynomial_interacting_proteins.json'), 'w', encoding="utf-8") as f:
        json.dump(polynomial_interacting_proteins, f, indent=2)
    print("Now Go test Syneromics.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        import traceback
        traceback.print_exc()