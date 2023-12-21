import collections
import json
import os
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import Utils
from LinearProgrammingSolution.GraphHandlers import load_patient_snps
from ResultsAnalysis.performance_evaluation import check_performances, plot_performances


def plot_gene_list_length_distribution(results_by_param):
    """
    This function plots the distribution of gene list lengths for every parameter value for each param.
    :param results_by_param: dictionary of param name to dictionary of results for parameter values.
    """
    sns.set()
    figure, axis = plt.subplots(1, len(results_by_param), figsize=(16, 9))

    for i, (param_name, results) in enumerate(results_by_param.items()):
        gene_list_lengths_by_param = {}
        plot = axis[i]
        for param_value, ranked_genes_lists in results.items():
            gene_list_lengths = [len(gene_list) for gene_list in ranked_genes_lists.values()]
            gene_list_lengths_by_param[param_value] = gene_list_lengths
        gene_list_lengths_by_param = dict(sorted(gene_list_lengths_by_param.items(), key=lambda item: item[0]))
        plot.boxplot(gene_list_lengths_by_param.values(), labels=gene_list_lengths_by_param.keys())
        plot.set_xlabel(param_name)
        plot.set_ylabel('Gene List Length')
        plot.set_title('Gene List Length Distribution for {}'.format(param_name))
    plt.show()

def plot_unique_genes_count_distribution(results, param_name):
    """
    This function plots the distribution of unique genes count for every parameter value.
    :param results: dictionary of results for parameter values.
    :param param_name: name of the parameter to plot the distribution for.
    """
    unqiue_genes_count = {}
    for param, ranked_genes_lists in results.items():
        unqiue_genes = set([gene for list in ranked_genes_lists.values() for gene in list])
        unqiue_genes_count[param] = len(unqiue_genes)
    unqiue_genes_count = dict(sorted(unqiue_genes_count.items(), key=lambda item: item[0], reverse=True))
    plt.plot(unqiue_genes_count.keys(), unqiue_genes_count.values(), '--o', markersize=5)
    plt.xlabel(param_name)
    plt.ylabel('Unique Genes Count')
    plt.title('Unique Genes Count Distribution')
    plt.show()

def plot_genes_occurrence_distribution(results, param_name):
    """
    This function plots the distribution of gene occurrences for every parameter value.
    :param results: dictionary of results for parameter values.
    :param param_name: name of the parameter to plot the distribution for.
    """
    gene_occurrences_by_param = {}
    for param, ranked_genes_lists in results.items():
        gene_occurrences = collections.Counter([gene for list in ranked_genes_lists.values() for gene in list])
        gene_occurrences_by_param[param] = [occurrences for gene, occurrences in gene_occurrences.items()]
    gene_occurrences_by_param = dict(sorted(gene_occurrences_by_param.items(), key=lambda item: item[0], reverse=True))
    plt.boxplot(gene_occurrences_by_param.values(), labels=gene_occurrences_by_param.keys())
    plt.xlabel(param_name)
    plt.ylabel('Gene occurrence')
    plt.title('Gene Occurrences')
    plt.show()

def plot_top_genes_occurrences_table(results_by_param,TOP_X = 6):
    """
    This function plots a table of occurrences across pations of the top genes for every parameter valuevalue for each param.
    :param results_by_param: dictionary of param name to dictionary of results for parameter values.
    :param TOP_X: number of top genes to include.
    """
    sns.set()
    figure, axis = plt.subplots(1, len(results_by_param), figsize=(16, 9))

    for i, (param_name, results) in enumerate(results_by_param.items()):
        plot = axis[i]
        top_genes_occurrences = {}
        for param, ranked_genes_lists in results.items():
            gene_occurrences = collections.Counter([gene for list in ranked_genes_lists.values() for gene in list])
            top_genes = {gene: occurrences for gene, occurrences in gene_occurrences.most_common(TOP_X)}
            top_genes_occurrences[param] = top_genes
        top_genes_occurrences = dict(sorted(top_genes_occurrences.items(), key=lambda item: item[0], reverse=True))
        df = pd.DataFrame(top_genes_occurrences)
        df = df.fillna(0)
        sns.heatmap(data=df,annot=True, cmap='Blues', fmt='g', ax=plot)
        plot.set_title(f'Top Genes Occurrences by {param_name}')
        plot.set_xlabel(param_name)
    plt.show()
    # top_genes_occurrences = {}
    # for param, ranked_genes_lists in results.items():
    #     gene_occurrences = collections.Counter([gene for list in ranked_genes_lists.values() for gene in list])
    #     top_genes = {gene: occurrences for gene, occurrences in gene_occurrences.most_common(TOP_X)}
    #     top_genes_occurrences[param] = top_genes
    # top_genes_occurrences = dict(sorted(top_genes_occurrences.items(), key=lambda item: item[0], reverse=True))
    # df = pd.DataFrame(top_genes_occurrences)
    # df = df.fillna(0)
    # sns.set(font_scale=0.8)
    # sns.heatmap(data=df,annot=True, cmap='Blues', fmt='g')
    # plt.show()


def parameters_analysis_main(results_dir):
    res = Utils.load_results(results_dir)
    res_by_alpha = {alpha: res[(alpha, beta, gamma)] for alpha, beta, gamma in res.keys() if (beta == 0) and (gamma == 1)}
    res_by_beta = {beta: res[(alpha, beta, gamma)] for alpha, beta, gamma in res.keys() if (alpha == 0) and (gamma == 1)}
    res_by_gamma = {gamma: res[(alpha, beta, gamma)] for alpha, beta, gamma in res.keys() if (alpha == 0) and (beta == 0)}
    res_by_param = {'alpha': res_by_alpha, 'beta': res_by_beta, 'gamma': res_by_gamma}
    # plot_gene_list_length_distribution(res_by_param)
    # plot_unique_genes_count_distribution(res_by_alpha, 'alpha')
    # plot_unique_genes_count_distribution(res_by_beta, 'beta')
    # plot_unique_genes_count_distribution(res_by_gamma, 'gamma')
    # plot_genes_occurrence_distribution(res_by_alpha, 'alpha')
    # plot_genes_occurrence_distribution(res_by_beta, 'beta')
    # plot_genes_occurrence_distribution(res_by_gamma, 'gamma')
    plot_top_genes_occurrences_table({'alpha': res_by_alpha, 'beta': res_by_beta}, TOP_X=10)
    # plot_top_genes_occurrences_table(res_by_beta, 'beta')
    # plot_top_genes_occurrences_table(res_by_gamma, 'gamma')


def get_best_perf(results_file: str, target_value_index: int = 0,
                  required_metrics=None) -> dict:
    required_metrics = ['precision', 'recall', 'f1'] if required_metrics is None else required_metrics
    # required_metrics_str = '_'.join(required_metrics)
    with(open(results_file)) as f:
        all_res_dict = json.load(f)['search_results']
    # get prodigy results
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    PRODIGY_results = json.load(open('./Data/PRODIGY_results.json'))
    patient_snps = load_patient_snps()
    PRODIGY_performances = check_performances(PRODIGY_results, patient_snps, gold_standard_drivers)
    # get best results
    if target_value_index <= 0 or target_value_index > len(gold_standard_drivers):
        target_value_index = 19
    better_counter = 0
    best_perf = {}
    for specific_res in all_res_dict.values():
        better = True
        for per_metric in required_metrics:
            if PRODIGY_performances[per_metric][target_value_index] > specific_res[per_metric][target_value_index]:
                better = False
        if better:
            better_counter += 1
            print(
                f'fount better perf than prodigy: {specific_res["alpha_param"]=}, {specific_res["beta_param"]=}, {specific_res["gamma_param"]=}')
            print(specific_res['recall'][target_value_index])
            if best_perf == {}:
                best_perf = specific_res
            else:
                best = True
                for per_metric in required_metrics:
                    if best_perf[per_metric][target_value_index] > specific_res[per_metric][target_value_index]:
                        best = False
                if best:
                    best_perf = specific_res
    print(f'found {better_counter} better results than prodigy')
    print(f'best perf was: {best_perf["alpha_param"]=}, {best_perf["beta_param"]=}, {best_perf["gamma_param"]=}')
    return best_perf


def plot_gamma_tradeoff(results_file, alpha_val, beta_val):
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    PRODIGY_results = json.load(open('./Data/PRODIGY_results.json'))
    patient_snps = load_patient_snps()
    PRODIGY_performances = check_performances(PRODIGY_results, patient_snps, gold_standard_drivers)
    with(open(results_file)) as f:
        all_res_dict = json.load(f)['search_results']
    all_res_dict = {key: value for key, value in all_res_dict.items() if
                    value['alpha_param'] == alpha_val and value['beta_param'] == beta_val}
    perf_dict = {k: {
        "precision": v["precision"],
        "recall": v["recall"],
        "f1": v["f1"]} for k, v in all_res_dict.items()}
    perf_dict['PRODIGY'] = PRODIGY_performances
    results_directory = os.path.dirname(results_file)
    plot_performances(perf_dict, is_sorted=False,
                      save_path=str(Path(results_directory) / f'gamma_comp_{alpha_val=}_{beta_val=}.png'))
