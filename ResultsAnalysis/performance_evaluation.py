import os
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


from LinearProgrammingSolution.GraphHandlers import load_patient_snps
import Utils

TOP_X = 20


def calculate_precision(ranked_list, gold_standard):
    """
    calculate the fraction of the first i genes that are known drivers for every length i
    :param ranked_list: a ranked list of potential drivers
    :param gold_standard: a list of known drivers
    :return: a vector of precision values for every length i
    """
    precision_vector = []
    for i in range(1, len(ranked_list) + 1):
        intersection = len(set(ranked_list[:i]).intersection(gold_standard))
        precision_vector.append(intersection / i)
    precision_vector.extend([precision_vector[-1]] * max(0, 100 - len(ranked_list)))
    return precision_vector


def calculate_recall(ranked_list, gold_standard):
    """
    caculate the fraction of known drivers that are in the first i genes for every length i
    :param ranked_list: a ranked list of potential drivers
    :param gold_standard: a list of known drivers
    :return: a vector of recall values for every length i
    """
    recall_vector = []
    for i in range(1, len(ranked_list) + 1):
        intersection = len(set(ranked_list[:i]).intersection(gold_standard))
        recall_vector.append(intersection / len(gold_standard))
    recall_vector.extend([recall_vector[-1]] * max(0, 100 - len(ranked_list)))
    return recall_vector


def check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers):
    """
    calculate the mean precision, recall and f1 for every length i upto TOP_X
    :param ranked_genes_lists: a dictionary of ranked genes lists for every patient
    :param patient_snps: a dictionary of snps for every patient
    :param gold_standard_drivers: a list of known drivers
    :return: vectors of mean precision, recall and f1 for every length i
    """
    precision_matrices = {}
    recall_matrices = {}
    f1_matrices = {}
    
    for patient in ranked_genes_lists.keys():
        ranked_genes_list = ranked_genes_lists[patient]
        patient_snp = patient_snps[patient]
        
        if not set(gold_standard_drivers).intersection(patient_snp):
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


def plot_performances(performances, is_sorted, save_path=None):
    """
    plot the performances of the different parameters
    :param performances: a dictionary of performances for every parameter combination and PRODIGY
    :param is_sorted: whether the dictionary is sorted by a parameter values and should be colored accordingly
    :param save_path: path to save the plot to or None to show the plot
    """
    sns.set()
    figure, axis = plt.subplots(1, 3, figsize=(16, 9))
    if is_sorted:
        colormap = plt.cm.get_cmap('RdBu', len(performances) -1)
    for i, metric in enumerate(['precision', 'recall', 'f1']):
        plot = axis[i]
        for j, (name, performance) in enumerate(performances.items()):
            lines = plot.plot(range(len(performance[metric])), performance[metric], label=name, marker='o',markersize=3)
            if is_sorted:
                lines[0].set_color(colormap(j) if name != 'PRODIGY' else 'black')
        plot.set_xticks(range(0, TOP_X, 2))
        plot.set_xlabel('Top N Genes')
        plot.set_ylabel('Average {}'.format(metric.capitalize()))
        plot.set_title(metric.capitalize())
        plot.set_box_aspect(9/16)
    handles, labels = axis[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='upper right')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()


def performance_evaluation_main(results_dir):
    patient_snps = load_patient_snps()
    PRODIGY_results = json.load(open('./Data/PRODIGY_results.json'))
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    res = Utils.load_results(results_dir)
    all_performances = {}
    for (alpha, beta, gamma), ranked_genes_lists in res.items():
        all_performances[(alpha, beta, gamma)] = check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers)
    PRODIGY_performances = check_performances(PRODIGY_results, patient_snps, gold_standard_drivers)
    all_performances = {alpha: all_performances[(alpha, beta, gamma)] for alpha, beta, gamma in all_performances.keys() if (beta == 0) and (gamma == 1)}
    all_performances = dict(sorted(all_performances.items(), key=lambda item: item[0]))
    all_performances['PRODIGY'] = PRODIGY_performances
    plot_performances(all_performances, is_sorted=True)