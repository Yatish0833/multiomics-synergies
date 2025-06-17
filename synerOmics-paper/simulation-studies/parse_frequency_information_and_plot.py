import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def parse_ranger_results():
    base_path = '/Syneromics/results/Ranger/'
    interaction_models = ["Interaction_only", "Modifier_protein", "No_Interaction", "Dominant", "Redundant", "Synergistic", "Cubic", "Sigmoidal", "Exponential"]
    orders = [2]
    heritabilities = [10, 30, 50, 70, 100]
    versions = range(1, 11)  # 1 to 10

    results = []

    for model in interaction_models:
        for order in orders:
            for her in heritabilities:
                for version in versions:
                    folder_name = f"{model}_order_{order}_phenotypes_ranger_results_her{her}_V{version}"
                    folder_path = os.path.join(base_path, folder_name)
                    
                    if not os.path.exists(folder_path):
                        continue

                    tsv_file = f"{model}_order_{order}_occurrences.tsv"
                    tsv_path = os.path.join(folder_path, tsv_file)
                    
                    if not os.path.exists(tsv_path):
                        continue

                    df = pd.read_csv(tsv_path, sep='\t')
                    for _, row in df.iterrows():
                        proteins = row['Protein1'].split(';')[0] + '_' + row['Protein2'].split(';')[0]
                        rank = 0 if pd.isna(row['Rank']) or row['Rank'] == 'N/A' else int(row['Rank'])
                        results.append({
                            'model': model,
                            'order': order,
                            'heritability': her,
                            'version': version,
                            'variable_id': proteins,
                            'count': row['Count'],  # Number of times the interaction combo was found
                            'rank': rank
                        })

    return pd.DataFrame(results)

def save_results(df):
    df.to_csv('parent_child_new_ranger_results.csv', index=False)
    print("Results saved to parent_child_new_ranger_results.csv")

def plot_horizontal_bar_charts(df, order, variable='rank', output_filename_prefix='horizontal_bar'):
    """
    Creates horizontal bar charts for the specified variable, with one subplot per heritability.
    
    Parameters:
    - df: DataFrame with the data
    - order: Order value to filter by
    - variable: Variable to plot (e.g., 'rank' or 'count')
    - output_filename_prefix: Prefix for the output file name
    """
    # Filter data for the specified order
    df_order = df[df['order'] == order]
    
    if df_order.empty:
        print(f"No data available for order {order}")
        return
    

    print(f"\n\n=== Filtering Summary for {variable.capitalize()} ===")
    print("Model | Heritability | Total Points | Non-Zero Points")
    print("-" * 65)
    
    for model in sorted(df_order['model'].unique()):
        for her in sorted(df_order['heritability'].unique()):
            subset = df_order[(df_order['model'] == model) & (df_order['heritability'] == her)][variable]
            total_points = len(subset)
            non_zero_points = len(subset[subset > 0])
            
            print(f"{model:<20} | {her:<12} | {total_points:<12} | {non_zero_points:<15}")
    
    # Define a function to calculate stats based on the variable
    def calc_stats(x):
        if variable == 'rank':
            # For rank, exclude zeros for median calculation
            non_zero = x[x != 0]
            return pd.Series({
                'mean': non_zero.median() if len(non_zero) > 0 else 0
            })
        else:  # variable == 'count'
            # For count, exclude zeros for median calculation
            non_zero = x[x != 0]
            return pd.Series({
                'mean': non_zero.median() if len(non_zero) > 0 else 0
            })
    
    stats = df_order.groupby(['model', 'heritability'])[variable].apply(calc_stats).reset_index()
    stats = stats.set_index(['model', 'heritability', 'level_2']).unstack(level='level_2')
    stats.columns = stats.columns.droplevel(0)  
    stats = stats.reset_index()
    pivoted_data = stats.pivot(index='model', columns='heritability', values='mean')
    heritabilities = sorted(df_order['heritability'].unique())

    fig, axes = plt.subplots(nrows=1, ncols=len(heritabilities), figsize=(20, 8), sharey=True)
    colors = sns.color_palette("colorblind", n_colors=len(pivoted_data.index))
    color_map = dict(zip(pivoted_data.index, colors))
    
    
    # Here we'll sort by the overall mean across all heritabilities for consistency
    overall_means = pivoted_data.mean(axis=1).sort_values(ascending=True)
    consistent_model_order = overall_means.index.tolist()
    
    for i, her in enumerate(heritabilities):
        ax = axes[i]

        her_data = pivoted_data.loc[consistent_model_order, her]
        y_pos = np.arange(len(consistent_model_order))
        values = her_data.values
        bars = ax.barh(y_pos, values, color=[color_map[model] for model in consistent_model_order])

        for j, (model, value) in enumerate(zip(consistent_model_order, values)):
            if value > 0:  
                ax.text(value + 0.1, j, f"{value:.1f}", va='center')

        ax.set_title(f'Heritability: {her}', fontsize=14, fontweight='bold')
    
        ax.set_yticks(y_pos)
        ax.set_yticklabels(consistent_model_order, fontsize=12)
        if i == 0:
            ax.set_ylabel('Interaction Model', fontsize=14, fontweight='bold')
        ax.set_xlabel(f'Median {variable.capitalize()} (excluding zeros)', fontsize=12, fontweight='bold')
        ax.grid(axis='x', linestyle='--', alpha=0.7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    variable_name = "Rank (excluding zeros)" if variable == 'rank' else "Count (excluding zeros)"

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    output_filename = f"{output_filename_prefix}_{variable}_order_{order}.pdf"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")
    
    return fig, axes

def plot_detection_rate_by_heritability(df, order, output_filename_prefix='detection_rate_heritability'):
    """
    Creates a grid of subplots showing detection rate as a function of heritability for each model separately.
    
    Parameters:
    - df: DataFrame with the data
    - order: Order value to filter by
    - output_filename_prefix: Prefix for the output file name
    """
    # Filter data for the specified order
    df_order = df[df['order'] == order]
    
    if df_order.empty:
        print(f"No data available for order {order}")
        return
    
    # Calculate detection rate (percentage of non-zero values) for each model and heritability
    detection_rates = df_order.groupby(['model', 'heritability']).apply(
        lambda x: (x['count'] >=2 ).mean() * 100
    ).reset_index(name='detection_rate')
    

    models = sorted(detection_rates['model'].unique())
    n_models = len(models)
    
    
    n_cols = 3  # 3 columns for 9 models
    n_rows = (n_models + n_cols - 1) // n_cols  # Ceiling division
    
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10), sharey=True, sharex=True)
    axes = axes.flatten()  
    
    
    colors = sns.color_palette("colorblind", n_colors=n_models)
    color_map = dict(zip(models, colors))
    

    for i, model in enumerate(models):
        ax = axes[i]
        model_data = detection_rates[detection_rates['model'] == model]
        
        # Plot detection rate vs heritability for this model
        ax.plot(model_data['heritability'], model_data['detection_rate'], 
                marker='o', linestyle='-', linewidth=2, markersize=8,
                color=color_map[model])
        
        # Calculate weighted arithmetic mean detection rate for this model
        mean_detection_rate = (model_data['detection_rate'] / model_data['heritability']).sum() / (1 / model_data['heritability']).sum()
        
        
        ax.set_title(f"{model} \n(Overall Detection: {mean_detection_rate:.1f}%)", fontsize=12, fontweight='bold')
        

        ax.grid(True, linestyle='--', alpha=0.7)
        

        ax.set_ylim(0, 105)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        

        if i % n_cols == 0:  
            ax.set_ylabel('Avg Detection (%) (Count >= 2)', fontsize=10)
        if i >= n_models - n_cols:  
            ax.set_xlabel('Heritability', fontsize=12, fontweight='bold')
    

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    

    output_filename = f"{output_filename_prefix}_order_{order}.pdf"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")
    

    return fig, axes

def print_data_summary(df):
    print("\nData Summary:")
    print(f"Total number of records: {len(df)}")
    print("\nRecords per order:")
    print(df['order'].value_counts())
    print("\nRecords per model:")
    print(df['model'].value_counts())
    print("\nRank statistics:")
    print(df['rank'].describe())
    print("\nCount statistics:")
    print(df['count'].describe())
    

    print("\nNumber of zero ranks per interaction model:")
    zero_ranks_per_model = df[df['rank'] == 0]['model'].value_counts().sort_index()
    for model, count in zero_ranks_per_model.items():
        print(f"{model}: {count}")
    
    # Total number of zero ranks
    print(f"\nTotal number of zero ranks: {zero_ranks_per_model.sum()}")
    
    print("\nMedian rank (excluding zeros) per interaction model and heritability:")
    median_ranks = df[df['rank'] > 0].groupby(['model', 'heritability'])['rank'].median().unstack(level='heritability')
    print(median_ranks)
    
    print("\nMedian count per interaction model and heritability:")
    median_counts = df.groupby(['model', 'heritability'])['count'].median().unstack(level='heritability')
    print(median_counts)

if __name__ == "__main__":
    # Option 1: Parse the results from original data sources
    results_df = parse_ranger_results()
    
    # Option 2: Load from CSV if already parsed
    # results_df = pd.read_csv('parent_child_new_ranger_results.csv')
    
    if not results_df.empty:
        save_results(results_df)
        print_data_summary(results_df)
        
        
        plot_horizontal_bar_charts(results_df, order=2, variable='rank', 
                                  output_filename_prefix='horizontal_bar')
        
        
        plot_horizontal_bar_charts(results_df, order=2, variable='count', 
                                  output_filename_prefix='horizontal_bar')
        

        plot_detection_rate_by_heritability(results_df, order=2, 
                           output_filename_prefix='detection_rate_heritability')
    else:
        print("No data was parsed. Please check the file paths and content.")