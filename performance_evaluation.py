import seaborn as sns
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from pathlib import Path

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
                ranked_genes_list[:min(20, len(ranked_genes_list))], gold_standard_drivers)[:20]
            
            intersected_genes = set(gold_standard_drivers).intersection(patient_snp)
            recall_matrices[patient] = calculate_recall(
                ranked_genes_list[:min(20, len(ranked_genes_list))], intersected_genes)[:20]
            
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

def plot_performances(performances):
    for name, performance in performances.items():
        plt.plot(performance['precision'], label=name)    
    plt.title('Precision')
    plt.legend()
    plt.show()
    for name, performance in performances.items():
        plt.plot(performance['recall'], label=name)    
    plt.title('Recall')
    plt.legend()
    plt.show()
    for name, performance in performances.items():
        plt.plot(performance['f1'], label=name)    
    plt.title('f1')
    plt.legend()
    plt.show()
    

print("loading data")
patient_snps = load_patient_snps()
ranked_genes_lists = json.load(open('./patients_with_ranked_genes_by_weight.json'))
PRODIGY_results = json.load(open('./Data/PRODIGY_results.json'))
global_ranked_genes_lists = json.load(open('./sorted_gene_names_by_weight.json'))
global_ranked_genes_lists
gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
print("calculating performances")
our_performances = check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers)
PRODIGY_performances = check_performances(PRODIGY_results, patient_snps, gold_standard_drivers)

global_performances = check_performances({'TCGA.A6.2671.01':global_ranked_genes_lists}, patient_snps, gold_standard_drivers)
plot_performances({'our algotithem': our_performances, 'PRODIGY': PRODIGY_performances})#, 'global': global_performances})