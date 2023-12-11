import os
import json
import collections
import pandas as pd
from pathlib import Path
from statistics import mean 

def calc_global_rank(data_dir):
    gene_scores = {}
    for file in Path(data_dir).glob('*.csv'):
        df = pd.read_csv(file, index_col=0)
        for pathway in df.index:
            for snv, weight in df.loc[pathway].items():
                if weight > 0:
                    if snv in gene_scores:
                        gene_scores[snv].add(file.stem)
                    else:
                        gene_scores[snv] = set(file.stem)
    return {snv: len(gene_scores[snv]) for snv in gene_scores.keys()}
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

def get_common_snvs(lists):
    common_snvs = []
    for i in range(max([len(lst) for lst in lists])):
        items_at_position = [lst[i] if len(lst) > i else None for lst in lists]
        most_common_item = max(set(items_at_position), key=items_at_position.count)
        common_snvs.append(most_common_item)
    return common_snvs
results_directory = r"ParamOptimizationResults/12_10_2023_06_57"
# results_directory = r"ParamOptimizationResults/12_09_2023_16_14"

ranked_gene_lists = {}
for subdir in os.listdir(results_directory):
    if str(subdir).startswith("gene_"):
        with open(os.path.join(results_directory, subdir, "ranked_genes_lists.json"), 'r') as f:
            data = json.load(f)
            for patient, ranked_gene_list in data.items():
                if patient in ranked_gene_lists:
                    ranked_gene_lists[patient].append(ranked_gene_list)
                else:
                    ranked_gene_lists[patient] =  [ranked_gene_list]

# Calculate agreement per identifier
avg_patient_agreements = {patient: calculate_agreement(lists) for patient, lists in ranked_gene_lists.items()}

patient_postion_distributions = {patient: calculate_postion_distribution(lists) for patient, lists in ranked_gene_lists.items()}

# Calculate agreement per position
avg_position_agreements = []
for i in range(max([len(lst) for lst in avg_patient_agreements.values()])):
    avg_agreements = mean([lst[i] for lst in avg_patient_agreements.values() if len(lst) > i])
    avg_position_agreements.append(avg_agreements)

data_dir=r'Data\DriverMaxSetApproximation\BaseData'
global_rank = calc_global_rank(data_dir)

global_rank = sorted(global_rank, key=lambda x: global_rank[x], reverse=True)
patient_ranks = {}
for file in Path(data_dir).glob('*.csv'):
    df = pd.read_csv(file, index_col=0)
    snvs = df.columns.values.tolist()
    snvs = sorted(snvs, key=lambda x: global_rank.index(x) if x in global_rank else len(global_rank))
    patient_ranks[file.stem] = snvs

agreements = {}
for patient, rank in patient_ranks.items():
    common_rank = get_common_snvs(ranked_gene_lists[patient])
    agreements[patient] = calculate_agreement([rank, common_rank])    

avg_position_agreements = []
for i in range(max([len(lst) for lst in agreements.values()])):
    avg_agreements = mean([lst[i] for lst in agreements.values() if len(lst) > i])
    avg_position_agreements.append(avg_agreements)
print(get_common_snvs(ranked_gene_lists['TCGA.DM.A28H.01']))
print(patient_ranks['TCGA.DM.A28H.01'])
y = calc_global_rank(data_dir)
print(dict(sorted(y.items(), key=lambda x:x[1],reverse=True)))
# for patient, rank in patient_ranks.items():
#     print(f"{patient}:")
#     print(rank)
#     print(get_common_snvs(ranked_gene_lists[patient]))    
# print("Agreement per identifier:")
# for identifier, agreement_vector in agreements.items():
#     print(f"{identifier}: {agreement_vector}")

# # # print("Postion distribution per identifier:")
# # # for identifier, postion_distribution  in patient_postion_distributions.items():
# # #     print(f"{identifier}: {postion_distribution}")


# print("\nAgreement per position:")
# print(avg_position_agreements)
