import json
import sys
import csv
import itertools
import os

def extract_proteins(json_data, model, num_proteins):
    for item in json_data:
        if item['model'] == model and len(item['proteins']) == num_proteins:
            proteins = item['proteins']
            correlated = {p: item['correlated_proteins'].get(p, []) for p in proteins}
            return proteins, correlated
    return None, None

def create_combinations(proteins, correlated):
    combinations = []
    for combo in itertools.combinations(proteins, len(proteins)):
        base_combo = list(combo)
        combinations.append((tuple(base_combo), [None] * len(base_combo)))
        for i, protein in enumerate(base_combo):
            for corr in correlated[protein]:
                corr_protein, corr_value = corr.split(':')
                new_combo = base_combo.copy()
                new_combo[i] = corr_protein
                corr_info = [None] * len(base_combo)
                corr_info[i] = (protein, corr_value)
                combinations.append((tuple(new_combo), corr_info))
    return combinations

def load_ranger_data(file_path):
    data = {}
    #rankings = {}
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        #for rank, row in enumerate(sorted(reader, key=lambda x: float(x[-1]), reverse=True), 1):
        #    proteins = tuple(sorted(row[:-1]))
        #    count = float(row[-1])
        #    data[proteins] = count
        #    rankings[proteins] = rank
        for row in reader:
            proteins = tuple(sorted(row[:-1]))
            count = float(row[-1])
            data[proteins] = count

    unique_counts = sorted(set(data.values()), reverse=True)


    count_to_rank = {count: i+1 for i, count in enumerate(unique_counts)}
    
    # Assign ranks to proteins based on their count
    rankings = {proteins: count_to_rank[data[proteins]] for proteins in data}
    return data, rankings

def format_output(combo, corr_info, true_proteins):
    output = []
    for protein, corr in zip(combo, corr_info):
        if protein in true_proteins:
            output.append(f"{protein} (True)")
        elif corr:
            output.append(f"{protein} (Correlated to {corr[0]}, r={corr[1]})")
        else:
            output.append(protein)
    return output

def write_results(combinations, ranger_data, ranger_rankings, true_proteins, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        headers = [f'Protein{i+1}' for i in range(len(next(iter(combinations))[0]))] + ['Count', 'Rank']
        writer.writerow(headers)

        for combo, corr_info in combinations:
            sorted_combo = tuple(sorted(combo))
            count = ranger_data.get(sorted_combo, 0)
            rank = ranger_rankings.get(sorted_combo, 'N/A')
            formatted_combo = format_output(combo, corr_info, true_proteins)
            writer.writerow(formatted_combo + [count, rank])

def main(json_file, ranger_results_file, model_type, num_proteins, pair_file=None):
    with open(json_file, 'r') as f:
        json_data = json.load(f)

    proteins, correlated = extract_proteins(json_data, model_type, num_proteins)
    if not proteins:
        print(f"No matching model '{model_type}' with {num_proteins} proteins found in the JSON file.")
        return

    combinations = create_combinations(proteins, correlated)
    ranger_data, ranger_rankings = load_ranger_data(ranger_results_file)
    
    output_dir = os.path.dirname(ranger_results_file)
    base_output_file = os.path.join(output_dir, f"{model_type}_order_{num_proteins}")

    # Main combinations (pairs or triplets)
    main_output_file = f"{base_output_file}_occurrences.tsv"
    write_results(combinations, ranger_data, ranger_rankings, proteins, main_output_file)
    print(f"Main results written to: {main_output_file}")

    if num_proteins == 3 and pair_file:
        pair_data, pair_rankings = load_ranger_data(pair_file)
        sub_pairs = []
        for combo, corr_info in combinations:
            for i, j in itertools.combinations(range(3), 2):
                sub_pair = (combo[i], combo[j])
                sub_corr = (corr_info[i], corr_info[j])
                sub_pairs.append((sub_pair, sub_corr))
        
        sub_pair_output_file = f"{base_output_file}_subpair_occurrences.tsv"
        write_results(sub_pairs, pair_data, pair_rankings, proteins, sub_pair_output_file)
        print(f"Sub-pair results written to: {sub_pair_output_file}")

if __name__ == "__main__":
    if len(sys.argv) not in [5, 6]:
        print("Usage: python script.py <json_file> <ranger_results_file> <model_type> <num_proteins> [pair_file]")
        sys.exit(1)

    json_file, ranger_results_file, model_type, num_proteins = sys.argv[1:5]
    pair_file = sys.argv[5] if len(sys.argv) == 6 else None
    main(json_file, ranger_results_file, model_type, int(num_proteins), pair_file)