import json
import sys
import csv
import os

def load_json_data(json_file):
    """Load the JSON data containing model information."""
    with open(json_file, 'r') as f:
        json_data = json.load(f)
    

    if json_data and isinstance(json_data[0], list):
        json_data = json_data[0]
    
    return json_data

def extract_ground_truth_pairs(json_data):
    """Extract all protein pairs that are used in the models."""
    protein_pairs = []
    
    # Count the number of models in the JSON
    # This will be the same for all pairs from this JSON
    model_count = len(json_data)
    
    for model_data in json_data:
        model_type = model_data['model']
        proteins = model_data.get('proteins', [])
        
        # Only process if we have at least 2 proteins
        if len(proteins) >= 2:
            # Create pairs from the proteins list
            for i in range(len(proteins) - 1):
                for j in range(i + 1, len(proteins)):
                    protein_pairs.append({
                        'protein1': proteins[i],
                        'protein2': proteins[j],
                        'model': model_type,
                        'protein_set': model_count  # Add the count of models
                    })
    
    return protein_pairs

def load_ranger_counts(counts_file):
    """Load the protein pair counts from the Ranger results."""
    counts_data = {}
    
    with open(counts_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  
        
        for row in reader:
            if len(row) >= 3:  # var1, var2, count
                var1, var2, count = row[0], row[1], float(row[2])
                
                # Store both orderings of the pair for easier matching
                pair1 = (var1, var2)
                pair2 = (var2, var1)
                
                counts_data[pair1] = count
                counts_data[pair2] = count
    
    # Compute rankings based on counts
    unique_counts = sorted(set(counts_data.values()), reverse=True)
    count_to_rank = {count: i+1 for i, count in enumerate(unique_counts)}
    
    ranks_data = {pair: count_to_rank[count] for pair, count in counts_data.items()}
    
    return counts_data, ranks_data

def match_ground_truth_with_counts(protein_pairs, counts_data, ranks_data):
    """Match ground truth protein pairs with their counts in the Ranger results."""
    matched_results = []
    
    for pair_info in protein_pairs:
        protein1 = pair_info['protein1']
        protein2 = pair_info['protein2']
        model = pair_info['model']
        protein_set = pair_info['protein_set']  
        
        # Check if this pair appears in the counts data
        pair = (protein1, protein2)
        reverse_pair = (protein2, protein1)
        
        count = counts_data.get(pair) or counts_data.get(reverse_pair) or 0
        rank = ranks_data.get(pair) or ranks_data.get(reverse_pair) or 'N/A'
        
        matched_results.append({
            'protein1': protein1,
            'protein2': protein2,
            'model': model,
            'protein_set': protein_set,  # Include in the results
            'count': count,
            'rank': rank
        })
    
    # Sort by count (descending)
    matched_results.sort(key=lambda x: x['count'] if isinstance(x['count'], (int, float)) else 0, reverse=True)
    
    return matched_results

def write_results_to_csv(results, output_file):
    """Write the matched results to a CSV file."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        header = ['Protein1', 'Protein2', 'Model', 'Protein_Set', 'Count', 'Rank']
        writer.writerow(header)
        
        # Write data rows
        for result in results:
            row = [
                result['protein1'],
                result['protein2'],
                result['model'],
                result['protein_set'],
                result['count'],
                result['rank']
            ]
            writer.writerow(row)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <json_file> <ranger_counts_file>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    counts_file = sys.argv[2]
    

    output_dir = os.path.dirname(counts_file)
    output_file = os.path.join(output_dir, "ground_truth_matches.csv")
    

    json_data = load_json_data(json_file)
    protein_pairs = extract_ground_truth_pairs(json_data)
    counts_data, ranks_data = load_ranger_counts(counts_file)
    
    # Match ground truth with counts
    matched_results = match_ground_truth_with_counts(protein_pairs, counts_data, ranks_data)
    

    write_results_to_csv(matched_results, output_file)
    
    print(f"Results written to: {output_file}")
    print(f"Found {len(matched_results)} protein pairs from ground truth.")
    print(f"Number of pairs with non-zero counts: {sum(1 for r in matched_results if r['count'] > 0)}")

if __name__ == "__main__":
    main()