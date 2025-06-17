import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap

def parse_ground_truth_matches():
    """
    Parse all ground_truth_matches.csv files from different directories
    and collate them with version and heritability information.
    """
    base_path = '/Syneromics/scripts/results/Ranger/'
    file_pattern = '**/ground_truth_matches.csv'
    
    all_results = []
    

    for file_path in glob.glob(os.path.join(base_path, file_pattern), recursive=True):
        # Extract folder information to get heritability and version
        folder_path = os.path.dirname(file_path)
        folder_name = os.path.basename(folder_path)
        
        try:
            # Expected format: {model}_order_{order}_phenotypes_ranger_results_her{her}_V{version}
            # or complex_phenotypes_ranger_results_her{her}_V{version}
            parts = folder_name.split('_')
            
            
            her_part = next((part for part in parts if part.startswith('her')), None)
            if her_part:
                her = int(her_part.replace('her', ''))
            else:
                
                her_index = next((i for i, part in enumerate(parts) if 'her' in part), None)
                if her_index is not None and her_index + 1 < len(parts):
                    her_str = parts[her_index + 1].split('_')[0]
                    her = int(her_str)
                else:
                    
                    for part in parts:
                        if 'her' in part:
                            her_str = ''.join(c for c in part if c.isdigit())
                            if her_str:
                                her = int(her_str)
                                break
                    else:
                        print(f"Could not find heritability in folder name: {folder_name}")
                        continue
            
            
            version_part = next((part for part in parts if part.startswith('V') and len(part) > 1), None)
            if version_part:
                version = int(version_part.replace('V', ''))
            else:
               
                v_index = next((i for i, part in enumerate(parts) if part == 'V'), None)
                if v_index is not None and v_index + 1 < len(parts):
                    version_str = parts[v_index + 1].split('_')[0]
                    version = int(version_str)
                else:
                   
                    for part in parts:
                        if part.startswith('V') and len(part) > 1:
                            version_str = ''.join(c for c in part if c.isdigit())
                            if version_str:
                                version = int(version_str)
                                break
                    else:
                        print(f"Could not find version in folder name: {folder_name}")
                        continue
            
           
            if 'complex' in folder_name:
                order = 2  # Default order for complex
            else:
                
                order_index = next((i for i, part in enumerate(parts) if part == 'order'), None)
                if order_index is not None and order_index + 1 < len(parts):
                    order = int(parts[order_index + 1])
                else:
                    order = 2 
            
            
            df = pd.read_csv(file_path)
            
            # The Model column already exists in the file, so we don't need to add it
            # Add heritability, version, and order metadata
            df['heritability'] = her
            df['version'] = version
            df['order'] = order
            
            all_results.append(df)
            
        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            continue

    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        combined_df['Rank'] = pd.to_numeric(combined_df['Rank'], errors='coerce')
        combined_df['Count'] = pd.to_numeric(combined_df['Count'], errors='coerce')
        
        return combined_df
    else:
        print("No valid data files found.")
        return pd.DataFrame()

def save_results(df, filename='/Users/kap037/Desktop/Syneromics/scripts/complex-phenotypes-5-march-2025/Complex_ground_truth_matches_collated.csv'):
    """Save the collated results to a CSV file."""
    df.to_csv(filename, index=False)
    print(f"Results saved to {filename}")

def plot_rank_by_heritability(df):
    """Create bar chart with heritability on x-axis, rank on y-axis, segregated by model."""
    
    if df.empty:
        print("No data available for plotting")
        return

    df_filtered = df[(df['Count'] > 0) & (~df['Rank'].isna())]
    
    if df_filtered.empty:
        print("No valid data for rank plotting after filtering")
        return
    
    # Calculate median rank for each model and heritability
    df_filtered = df[(df['Rank'] > 0) & (~df['Rank'].isna())]
    rank_stats = df_filtered.groupby(['Model', 'heritability'])['Rank'].agg(['median']).reset_index()
    

    pivot_median = rank_stats.pivot(index='Model', columns='heritability', values='median')
    

    plt.figure(figsize=(16, 10))
    

    heritabilities = sorted(df_filtered['heritability'].unique())
    models = pivot_median.index.tolist()
    

    colors = sns.color_palette("colorblind", n_colors=len(models))
    version_cmap = LinearSegmentedColormap.from_list("version_colors", colors)
    

    group_width = 0.8

    bar_width = group_width / len(models)
    

    x = np.arange(len(heritabilities))
    

    for i, model in enumerate(models):
        median_values = pivot_median.loc[model].values

        positions = x + (i - len(models)/2 + 0.5) * bar_width
        

        bars = plt.bar(positions, median_values, width=bar_width, label=model, color=colors[i])
    

    plt.xlabel('Heritability', fontsize=14, labelpad=10, fontweight='bold')
    plt.ylabel('Median Rank (lower is better)', fontsize=14, labelpad=10, fontweight='bold')
    plt.title('Median Rank by Heritability and Model', fontsize=16, fontweight='bold')
    

    plt.xticks(x, heritabilities, fontsize=12)
    plt.yticks(fontsize=12)
    

    plt.legend(title='Model', fontsize=10, title_fontsize=12, loc='upper center', 
               bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=True)
    

    if pivot_median.max().max() / pivot_median.min().min() > 10:
        plt.yscale('log')
        
        plt.ylim(bottom=0.9)
    
   
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    

    plt.tight_layout()
    plt.savefig('/Users/kap037/Desktop/Syneromics/scripts/complex-phenotypes-5-march-2025/rank_by_heritability.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Rank plot saved as rank_by_heritability.png")
    

    print("\nMedian Rank Data:")
    print(pivot_median)

def plot_count_by_heritability(df):
    """Create bar chart with heritability on x-axis, counts on y-axis, segregated by model."""
    if df.empty:
        print("No data available for plotting")
        return
    
    # Calculate median counts for each model and heritability
    df_filtered_counts = df[(df['Count'] > 0) & (~df['Count'].isna())]
    count_stats = df_filtered_counts.groupby(['Model', 'heritability'])['Count'].median().reset_index()
    
    # Create pivot table for median
    pivot_median = count_stats.pivot(index='Model', columns='heritability', values='Count')
    
    
    plt.figure(figsize=(16, 10))
    

    heritabilities = sorted(df['heritability'].unique())
    models = pivot_median.index.tolist()
    

    colors = sns.color_palette("colorblind", n_colors=len(models))
    version_cmap = LinearSegmentedColormap.from_list("version_colors", colors)
    

    group_width = 0.8

    bar_width = group_width / len(models)
    

    x = np.arange(len(heritabilities))
    

    for i, model in enumerate(models):
        median_values = pivot_median.loc[model].values

        positions = x + (i - len(models)/2 + 0.5) * bar_width
        

        bars = plt.bar(positions, median_values, width=bar_width, label=model, color=colors[i])
    

    plt.axhline(y=2, color='red', linestyle='--', linewidth=2, alpha=0.7)
    plt.text(x[-1] + 0.5, 2.1, 'Count = 2', color='red', fontweight='bold')
    

    plt.xlabel('Heritability', fontsize=14, labelpad=10, fontweight='bold')
    plt.ylabel('Median Count', fontsize=14, labelpad=10, fontweight='bold')
    plt.title('Median Count by Heritability and Model', fontsize=16, fontweight='bold')
    
 
    plt.xticks(x, heritabilities, fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.legend(title='Model', fontsize=10, title_fontsize=12, loc='upper center', 
               bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=True)
    
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    

    plt.tight_layout()
    plt.savefig('/Users/kap037/Desktop/Syneromics/scripts/complex-phenotypes-5-march-2025/count_by_heritability.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Count plot saved as count_by_heritability.png")
    

    print("\nMedian Count Data:")
    print(pivot_median)

def plot_protein_set_panels(df):
    """
    Create panels of count vs heritability bar charts segregated by Protein_Set.
    Max 4 panels, one for each unique Protein_Set value.
    """

    if df.empty or 'Protein_Set' not in df.columns:
        print("No data or Protein_Set column available for plotting")
        return
    

    df['Protein_Set'] = pd.to_numeric(df['Protein_Set'], errors='coerce')

    df_filtered = df.dropna(subset=['Protein_Set'])
    
    if df_filtered.empty:
        print("No valid data for protein set panels after filtering")
        return
    

    protein_sets = sorted(df_filtered['Protein_Set'].unique())
    
    # Take at most 4 protein sets
    protein_sets = protein_sets[:min(4, len(protein_sets))]
    

    n_panels = len(protein_sets)
    fig, axs = plt.subplots(1, n_panels, figsize=(5 * n_panels, 6), sharey=True)
    
    if n_panels == 1:
        axs = [axs]

    heritabilities = sorted(df_filtered['heritability'].unique())
    models = sorted(df_filtered['Model'].unique())
    
    colors = sns.color_palette("colorblind", n_colors=len(models))
    color_map = dict(zip(models, colors))
    

    for i, protein_set in enumerate(protein_sets):

        panel_data = df_filtered[df_filtered['Protein_Set'] == protein_set]
        
        if panel_data.empty:
            axs[i].text(0.5, 0.5, f'No data for Protein_Set={protein_set}', 
                       ha='center', va='center', fontsize=12)
            axs[i].set_title(f'Protein_Set = {int(protein_set)}', fontsize=14, fontweight='bold')
            continue
        
        # Calculate median count for each model and heritability
        count_stats = panel_data.groupby(['Model', 'heritability'])['Count'].median().reset_index()
        

        pivot_count = count_stats.pivot(index='Model', columns='heritability', values='Count')
        
        # Check if we have valid data after pivoting
        if pivot_count.empty:
            axs[i].text(0.5, 0.5, f'No valid data for Protein_Set={protein_set}', 
                       ha='center', va='center', fontsize=12)
            axs[i].set_title(f'Protein_Set = {int(protein_set)}', fontsize=14, fontweight='bold')
            continue
        
        group_width = 0.8

        bar_width = group_width / len(pivot_count.index)
        

        x = np.arange(len(pivot_count.columns))
        
        for j, model in enumerate(pivot_count.index):
            values = pivot_count.loc[model].values
            positions = x + (j - len(pivot_count.index)/2 + 0.5) * bar_width
            
            axs[i].bar(positions, values, width=bar_width, label=model if i == 0 else "", 
                      color=color_map[model])
        axs[i].axhline(y=2, color='red', linestyle='--', linewidth=2, alpha=0.7)
        if i == n_panels - 1:  
            axs[i].text(len(pivot_count.columns) - 1 + 0.5, 2.1, 'Count = 2', color='red', fontweight='bold')
        
        axs[i].set_title(f'Protein_Set = {int(protein_set)}', fontsize=14, fontweight='bold')
        axs[i].set_xticks(x)
        axs[i].set_xticklabels(pivot_count.columns, rotation=45, fontsize=10)
        axs[i].grid(axis='y', linestyle='--', alpha=0.7)

        if i == 0:
            axs[i].set_ylabel('Median Count', fontsize=14, fontweight='bold')
        
        axs[i].set_xlabel('Heritability', fontsize=12, fontweight='bold')
    
    handles, labels = [], []
    for model in models:
        handle = plt.Rectangle((0, 0), 1, 1, color=color_map[model])
        handles.append(handle)
        labels.append(model)
    
    fig.legend(handles, labels, 
              loc='lower center', 
              bbox_to_anchor=(0.5, -0.1),
              ncol=min(5, len(models)),
              fontsize=10,
              title='Model',
              title_fontsize=12)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)  

    plt.savefig('/Users/kap037/Desktop/Syneromics/scripts/complex-phenotypes-5-march-2025/protein_set_panels.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Protein set panels saved as protein_set_panels.png")
    
    for protein_set in protein_sets:
        print(f"\nMedian Count Data for Protein_Set = {int(protein_set)}:")
        panel_data = df_filtered[df_filtered['Protein_Set'] == protein_set]
        if not panel_data.empty:
            count_stats = panel_data.groupby(['Model', 'heritability'])['Count'].median().unstack()
            print(count_stats)
        else:
            print("No data available")

def print_data_summary(df):
    """Print summary statistics for the collated data."""
    print("\nData Summary:")
    print(f"Total number of records: {len(df)}")
    
    print("\nRecords per model:")
    print(df['Model'].value_counts())
    
    print("\nRecords per heritability:")
    print(df['heritability'].value_counts().sort_index())
    
    print("\nCount statistics:")
    print(df['Count'].describe())
    
    print("\nRank statistics (excluding zeros and NaN):")
    valid_ranks = df[(df['Count'] > 0) & (~df['Rank'].isna())]['Rank']
    if not valid_ranks.empty:
        print(valid_ranks.describe())
    else:
        print("No valid rank data available")
    
    print("\nZero counts per model:")
    zero_counts = df[df['Count'] == 0]['Model'].value_counts().sort_index()
    for model, count in zero_counts.items():
        print(f"{model}: {count}")
    
    print("\nDetection rate per model (pairs with Count > 0):")
    detection_rates = df.groupby('Model').apply(
        lambda x: (x['Count'] > 0).mean() * 100
    ).sort_values(ascending=False)
    
    for model, rate in detection_rates.items():
        print(f"{model}: {rate:.2f}%")
    
    print("\nDetection rate per heritability (pairs with Count > 0):")
    detection_rates = df.groupby('heritability').apply(
        lambda x: (x['Count'] > 0).mean() * 100
    ).sort_index()
    
    for her, rate in detection_rates.items():
        print(f"{her}: {rate:.2f}%")
        
    # Add detection rate by model and heritability
    print("\nDetection rate by model and heritability (pairs with Count > 0):")
    detection_matrix = df.pivot_table(
        index='Model', 
        columns='heritability',
        values='Count',
        aggfunc=lambda x: (x > 0).mean() * 100
    )
    print(detection_matrix)
    
    # Add median rank by model and heritability
    print("\nMedian rank by model and heritability (only pairs with Count > 0):")
    rank_matrix = df[df['Count'] > 0].pivot_table(
        index='Model', 
        columns='heritability',
        values='Rank',
        aggfunc='median'
    )
    print(rank_matrix)
    
    if 'Protein_Set' in df.columns:
        print("\nProtein_Set distribution:")
        protein_set_counts = df['Protein_Set'].value_counts().sort_index()
        for protein_set, count in protein_set_counts.items():
            print(f"Protein_Set {protein_set}: {count} records")
        
        # Add detection rate by Protein_Set
        print("\nDetection rate by Protein_Set (pairs with Count > 0):")
        protein_set_detection = df.groupby('Protein_Set').apply(
            lambda x: (x['Count'] > 0).mean() * 100
        ).sort_index()
        
        for protein_set, rate in protein_set_detection.items():
            print(f"Protein_Set {protein_set}: {rate:.2f}%")

if __name__ == "__main__":
    results_df = parse_ground_truth_matches()
    
    if not results_df.empty:
        save_results(results_df)
        print_data_summary(results_df)
        plot_rank_by_heritability(results_df)
        plot_count_by_heritability(results_df)

        if 'Protein_Set' in results_df.columns:
            plot_protein_set_panels(results_df)
        else:
            print("Protein_Set column not found in data. Skipping protein set panel visualization.")
    else:
        print("No data was parsed. Please check the file paths and content.")