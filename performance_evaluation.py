import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns

TOP_X = 20


def calculate_precision(ranked_list, gold_standard):
    precision_vector = []
    for i in range(1, len(ranked_list) + 1):
        intersection = len(set(ranked_list[:i]).intersection(gold_standard))
        precision_vector.append(intersection / i)
    precision_vector.extend([precision_vector[-1]] * max(0, 100 - len(ranked_list)))
    return precision_vector


def calculate_recall(ranked_list, gold_standard):
    recall_vector = []
    for i in range(1, len(ranked_list) + 1):
        intersection = len(set(ranked_list[:i]).intersection(gold_standard))
        recall_vector.append(intersection / len(gold_standard))
    recall_vector.extend([recall_vector[-1]] * max(0, 100 - len(ranked_list)))
    return recall_vector


def load_patient_snps():
    patient_snps = {}
    for file in Path('Data/DriverMaxSetApproximation/BaseData').glob('*.csv'):
        patient_name = file.stem
        df = pd.read_csv(file, index_col=0)
        patient_snps[patient_name] = df.columns.tolist()
    return patient_snps


def check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers):
    precision_matrices = {}
    recall_matrices = {}
    f1_matrices = {}
    
    for patient in ranked_genes_lists.keys():
        ranked_genes_list = ranked_genes_lists[patient]
        patient_snp = patient_snps[patient]
        
        if not set(gold_standard_drivers).intersection(patient_snp):
            print(f"No known drivers with SNV mutations for patient {patient}")
            continue
        
        if set(gold_standard_drivers).intersection(ranked_genes_list):
            precision_matrices[patient] = calculate_precision(
                ranked_genes_list[:min(TOP_X, len(ranked_genes_list))], gold_standard_drivers)[:TOP_X]
            
            intersected_genes = set(gold_standard_drivers).intersection(patient_snp)
            recall_matrices[patient] = calculate_recall(
                ranked_genes_list[:min(TOP_X, len(ranked_genes_list))], intersected_genes)[:TOP_X]
            
            curr_f1 = (2 * np.array(precision_matrices[patient]) * np.array(recall_matrices[patient])) / (
                        np.array(precision_matrices[patient]) + np.array(recall_matrices[patient]))
            curr_f1 = np.nan_to_num(curr_f1, nan=0)
            f1_matrices[patient] = curr_f1
    
    precision_matrix = np.array(list(precision_matrices.values()))
    recall_matrix = np.array(list(recall_matrices.values()))
    f1_matrix = np.array(list(f1_matrices.values()))
    
    precision_means = np.mean(precision_matrix, axis=0)
    recall_means = np.mean(recall_matrix, axis=0)
    f1_means = np.mean(f1_matrix, axis=0)
    
    return {
        "precision": precision_means,
        "recall": recall_means,
        "f1": f1_means
    }


def plot_performances(performances, save_path=None):
    sns.set()
    figure, axis = plt.subplots(1, 3, figsize=(16,9))
    plt1 = axis[0]
    plt2 = axis[1]
    plt3 = axis[2]

    for name, performance in performances.items():
        plt1.plot(range(len(performance['precision'])),performance['precision'], label=name, marker='o',markersize=3)
    plt1.set_xticks(range(0,20,2))
    plt1.set_xlabel('Top N Genes')
    plt1.set_ylabel('Average Precision')
    plt1.set_title('Precision')
    plt1.set_box_aspect(9/16)

    for name, performance in performances.items():
        plt2.plot(range(len(performance['recall'])),performance['recall'], label=name, marker='o',markersize=3)
    plt2.set_xticks(range(0,20,2))
    plt2.set_xlabel('Top N Genes')
    plt2.set_ylabel('Average Recall')
    plt2.set_title('Recall')
    plt2.set_box_aspect(9/16)

    for name, performance in performances.items():
        plt3.plot(range(len(performance['f1'])),performance['f1'], label=name, marker='o',markersize=3)
    plt3.set_xticks(range(0,20,2))
    plt3.set_xlabel('Top N Genes')
    plt3.set_ylabel('Average F1')
    plt3.set_title('F1')
    plt3.set_box_aspect(9/16)

    handles, labels = plt1.get_legend_handles_labels()
    figure.legend(handles, labels, loc='lower center')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()


if __name__ == '__main__':
    print("loading data")
    patient_snps = load_patient_snps()
    # ranked_genes_lists = json.load(open('./ranked_genes_lists.json'))
    # gene_wights = json.load(open('./gene_weights.json'))
    # wight_per_rank = {}
    # for patient in ranked_genes_lists.keys():
    #     for i, gene in enumerate(ranked_genes_lists[patient]):
    #         if i >= TOP_X:
    #             break
    #         wight_per_rank[i+1] = wight_per_rank.get(i+1, []) + [gene_wights.get(gene,0)]


    # #box plot of weights per rank
    # keys = list(wight_per_rank.keys())
    # values = list(wight_per_rank.values())

    # #Create a box plot
    # plt.boxplot(values, labels=keys)
    # plt.xlabel('Rank')
    # plt.ylabel('Weight')
    # plt.title('Gene Weights per Rank')
    # plt.grid(True)
    # plt.show()


    PRODIGY_results = json.load(open('./Data/PRODIGY_results.json'))
    # global_ranked_genes_lists = json.load(open('./sorted_gene_names_by_weight.json'))
    # global_ranked_genes_lists
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    print("calculating performances")
    all_performances = {}
    for d in os.listdir('./ParamOptimizationResults/12_11_2023_22_05/'):
        if d.startswith('alpha'):
            with open(os.path.join('./ParamOptimizationResults/12_11_2023_22_05/', d, 'ranked_genes_lists.json')) as f:
                ranked_genes_lists = json.load(f)
            all_performances[d.split('=')[1]] =check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers)
    PRODIGY_performances = check_performances(PRODIGY_results, patient_snps, gold_standard_drivers)
    all_performances['PRODIGY'] = PRODIGY_performances
    # global_performances = check_performances({'TCGA.A6.2671.01':global_ranked_genes_lists}, patient_snps, gold_standard_drivers)
    #all_performances = {k: all_performances[k] for k in ('-50', '-1','-0.01', '-0.1','0.01','0.1','1','50','PRODIGY')}
    plot_performances(all_performances)#, 'global': global_performances})

    # df = pd.DataFrame()
    # for patient in set(our_performances['recall'].keys()).intersection(PRODIGY_performances['recall'].keys()) :
    #     res = {}
    #     for i in range(min(min(TOP_X, len(our_performances['recall'][patient])), len(PRODIGY_performances['recall'][patient]))):
    #         res[i] = our_performances['recall'][patient][i] - PRODIGY_performances['recall'][patient][i]
    #     df[patient] = pd.Series(res)
    # df.T.to_csv('diff.csv')
