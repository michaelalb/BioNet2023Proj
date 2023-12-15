import os
import re
import json
import collections
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
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


def agreement_analysis():
    results_directory = r"ParamOptimizationResults/12_09_2023_16_14"  #  12_10_2023_06_57

    ranked_gene_lists = {}
    for subdir in os.listdir(results_directory):
        if str(subdir).startswith("gene_"):
            file_path = os.path.join(results_directory, subdir, "SingleRunResults/ranked_genes_lists.json")
            if not os.path.isfile(file_path):
                continue
            with open(file_path, 'r') as f:
                data = json.load(f)
                for pationt, ranked_gene_list in data.items():
                    if pationt in ranked_gene_lists:
                        ranked_gene_lists[pationt].append(ranked_gene_list)
                    else:
                        ranked_gene_lists[pationt] = [ranked_gene_list]

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

    print("Postion distribution per identifier:")
    for identifier, postion_distribution  in pationt_postion_distributions.items():
        print(f"{identifier}: {postion_distribution}")


    print("\nAgreement per position:")
    print(avg_position_agreements)

def get_points(results_directory):
    patients_per_gene = {}
    with open(os.path.join(results_directory, 'SingleRunResults/ranked_genes_lists.json')) as f:
        ranked_genes_lists = json.load(f)
    for patient, ranked_genes_list in ranked_genes_lists.items():
        for gene in ranked_genes_list:
            patients_per_gene[gene] = patients_per_gene.get(gene, []) + [patient]
    patients_per_gene = {gene: len(patients) for gene, patients in patients_per_gene.items()}
    petient_count_per_rank = {}
    with open(os.path.join(results_directory, 'SingleRunResults/ranked_genes_lists.json')) as f:
        ranked_genes_lists = json.load(f)
    for ranked_genes_list in ranked_genes_lists.values():
        for i, gene in enumerate(ranked_genes_list):
            petient_count_per_rank[i] = petient_count_per_rank.get(i, []) + [patients_per_gene[gene]]
    x = []
    y = []
    for i, petient_count in petient_count_per_rank.items():
         for j in petient_count:
             x.append(i)
             y.append(j)
    return(x,y)
    
def plot_gene_list_length_distribution(results_directory,file_regex ,parm_regex, param_name):
    gene_list_lengths_by_param = {}
    for subdir in Path(results_directory).glob(file_regex):
        with open(os.path.join(subdir, 'SingleRunResults/ranked_genes_lists.json')) as f:
            ranked_genes_lists = json.load(f)
        gene_list_lengths = [len(gene_list) for gene_list in ranked_genes_lists.values()]
        param = float(re.search(parm_regex, str(subdir)).group(1))
        gene_list_lengths_by_param[param] = gene_list_lengths
    gene_list_lengths_by_param = dict(sorted(gene_list_lengths_by_param.items(), key=lambda item: item[0], reverse=True))
    plt.boxplot(gene_list_lengths_by_param.values(), labels=gene_list_lengths_by_param.keys())
    plt.xlabel(param_name)
    plt.ylabel('Gene List Length')
    plt.title('Gene List Length Distribution')
    plt.show()

def plot_unique_genes_count_distribution(results_directory,file_regex ,parm_regex, param_name):
    unqiue_genes_count = {}
    for subdir in Path(results_directory).glob(file_regex):
        with open(os.path.join(subdir, 'SingleRunResults/ranked_genes_lists.json')) as f:
            ranked_genes_lists = json.load(f)
        unqiue_genes = set([gene for list in ranked_genes_lists.values() for gene in list])
        param = float(re.search(parm_regex, str(subdir)).group(1))
        unqiue_genes_count[param] = len(unqiue_genes)
    unqiue_genes_count = dict(sorted(unqiue_genes_count.items(), key=lambda item: item[0], reverse=True))
    plt.plot(unqiue_genes_count.keys(), unqiue_genes_count.values(), '--o', markersize=5)
    plt.xlabel(param_name)
    plt.ylabel('Unique Genes Count')
    plt.title('Unique Genes Count Distribution')
    plt.show()

def plot_genes_occurrence_distribution(results_directory,file_regex ,parm_regex, param_name):
    gene_occurrences_by_param = {}
    for subdir in Path(results_directory).glob(file_regex):
        with open(os.path.join(subdir, 'SingleRunResults/ranked_genes_lists.json')) as f:
            ranked_genes_lists = json.load(f)
        gene_occurrences = collections.Counter([gene for list in ranked_genes_lists.values() for gene in list])
        param = float(re.search(parm_regex, str(subdir)).group(1))
        gene_occurrences_by_param[param] = [occurrences for gene, occurrences in gene_occurrences.items()]
    gene_occurrences_by_param = dict(sorted(gene_occurrences_by_param.items(), key=lambda item: item[0], reverse=True))
    plt.boxplot(gene_occurrences_by_param.values(), labels=gene_occurrences_by_param.keys())
    plt.xlabel(param_name)
    plt.ylabel('Gene occurrence')
    plt.title('Gene Occurrences')
    plt.show()

if __name__ == "__main__":
    plot_gene_list_length_distribution(r"ParamOptimizationResults/12_13_2023_05_32",'alpha_param=*_beta_param=0',r'alpha_param=([^_]*)_.*','Alpha')
    plot_gene_list_length_distribution(r"ParamOptimizationResults/12_13_2023_05_32",'alpha_param=0_beta_param=*',r'.*beta_param=(.*)','Beta')
    plot_unique_genes_count_distribution(r"ParamOptimizationResults/12_13_2023_05_32",'alpha_param=*_beta_param=0',r'alpha_param=([^_]*)_.*','Alpha')
    plot_unique_genes_count_distribution(r"ParamOptimizationResults/12_13_2023_05_32",'alpha_param=0_beta_param=*',r'.*beta_param=(.*)','Beta')
    plot_genes_occurrence_distribution(r"ParamOptimizationResults/12_13_2023_05_32",'alpha_param=*_beta_param=0',r'alpha_param=([^_]*)_.*','Alpha')
    plot_genes_occurrence_distribution(r"ParamOptimizationResults/12_13_2023_05_32",'alpha_param=0_beta_param=*',r'.*beta_param=(.*)','Beta')
    # i = 0
    # results_directory = "ParamOptimizationResults/12_11_2023_20_02"
    # fig, axs = plt.subplots(len(os.listdir(results_directory)) - 1, figsize=(16,9),sharex=True)
    # for subdir in os.listdir(results_directory):
    #     if str(subdir).startswith("alpha_"):
    #         x,y = get_points(os.path.join(results_directory, subdir))
    #         axs[i].plot(x, y, 'o', markersize=1)
    #         axs[i].set_title(subdir)
    #         axs[i].set_ylim([0, 100])
    #         i += 1
    
    # plt.show()