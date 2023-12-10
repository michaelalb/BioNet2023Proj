import os
import json
import collections
from statistics import mean 

# Function to calculate agreement between lists in each postion
def calculate_agreement(lists):
    agreement_vector = []
    for i in range(max([len(lst) for lst in lists])):
        items_at_position = [lst[i] if len(lst) > i else None for lst in lists]
        most_common_item = max(set(items_at_position), key=items_at_position.count)
        agreement = items_at_position.count(most_common_item) / len(lists)
        agreement_vector.append(agreement)
    return agreement_vector

def calculate_postion_distribution(lists):
    distribution_vector = []
    for i in range(max([len(lst) for lst in lists])):
        items_at_position = [lst[i] if len(lst) > i else None for lst in lists]
        distribution_vector.append(collections.Counter(items_at_position))
    return distribution_vector

results_directory = r"ParamOptimizationResults/12_10_2023_06_57"

ranked_gene_lists = {}
for subdir in os.listdir(results_directory):
    if str(subdir).startswith("gene_"):
        with open(os.path.join(results_directory, subdir, "ranked_genes_lists.json"), 'r') as f:
            data = json.load(f)
            for pationt, ranked_gene_list in data.items():
                if pationt in ranked_gene_lists:
                    ranked_gene_lists[pationt].append(ranked_gene_list)
                else:
                    ranked_gene_lists[pationt] =  [ranked_gene_list]

# Calculate agreement per identifier
avg_pationt_agreements = {pationt: calculate_agreement(lists) for pationt, lists in ranked_gene_lists.items()}

pationt_postion_distributions = {pationt: calculate_postion_distribution(lists) for pationt, lists in ranked_gene_lists.items()}

# Calculate agreement per position
avg_position_agreements = []
for i in range(max([len(lst) for lst in avg_pationt_agreements.values()])):
    avg_agreements = mean([lst[i] for lst in avg_pationt_agreements.values() if len(lst) > i])
    avg_position_agreements.append(avg_agreements)

print("Agreement per identifier:")
for identifier, agreement_vector in avg_pationt_agreements.items():
    print(f"{identifier}: {agreement_vector}")

# print("Postion distribution per identifier:")
# for identifier, postion_distribution  in pationt_postion_distributions.items():
#     print(f"{identifier}: {postion_distribution}")


print("\nAgreement per position:")
print(avg_position_agreements)
