import collections
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import Utils
    
def plot_gene_list_length_distribution(results, param_name):
    """
    This function plots the distribution of gene list lengths for every parameter value.
    :param results: dictionary of results for parameter values.
    :param param_name: name of the parameter to plot the distribution for.
    """
    gene_list_lengths_by_param = {}
    for param, ranked_genes_lists in results.items():
        gene_list_lengths = [len(gene_list) for gene_list in ranked_genes_lists.values()]
        gene_list_lengths_by_param[param] = gene_list_lengths
    gene_list_lengths_by_param = dict(sorted(gene_list_lengths_by_param.items(), key=lambda item: item[0]))
    plt.boxplot(gene_list_lengths_by_param.values(), labels=gene_list_lengths_by_param.keys())
    plt.xlabel(param_name)
    plt.ylabel('Gene List Length')
    plt.title('Gene List Length Distribution')
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

def plot_top_genes_occurrences_table(results, param_name,TOP_X = 6):
    """
    This function plots a table of occurrences across pations of the top genes for every parameter value.
    :param results: dictionary of results for parameter values.
    :param param_name: name of the parameter to plot the table for.
    :param TOP_X: number of top genes to include.
    """
    top_genes_occurrences = {}
    for param, ranked_genes_lists in results.items():
        gene_occurrences = collections.Counter([gene for list in ranked_genes_lists.values() for gene in list])
        top_genes = {gene: occurrences for gene, occurrences in gene_occurrences.most_common(TOP_X)}
        top_genes_occurrences[param] = top_genes
    top_genes_occurrences = dict(sorted(top_genes_occurrences.items(), key=lambda item: item[0], reverse=True))
    df = pd.DataFrame(top_genes_occurrences)
    df = df.fillna(0)
    sns.set(font_scale=0.8)
    sns.heatmap(data=df,annot=True, cmap='Blues', fmt='g')
    plt.show()

def parameters_analysis_main():
    res = Utils.load_results(r"ParamOptimizationResults/12_13_2023_05_32")
    res_by_alpha = {alpha: res[(alpha, beta)] for alpha, beta in res.keys() if (beta == 0)}
    res_by_beta = {beta: res[(alpha, beta)] for alpha, beta in res.keys() if (alpha == 0)}
    res = Utils.load_results(r"ParamOptimizationResults/12_14_2023_22_20")
    res_by_gamma = {gamma: res[(alpha, beta, gamma)] for alpha, beta, gamma in res.keys() if (alpha == 0) and (beta == 0)}
    
    plot_gene_list_length_distribution(res_by_alpha, 'alpha')
    plot_gene_list_length_distribution(res_by_beta, 'beta')
    plot_gene_list_length_distribution(res_by_gamma, 'gamma')
    plot_unique_genes_count_distribution(res_by_alpha, 'alpha')
    plot_unique_genes_count_distribution(res_by_beta, 'beta')
    plot_unique_genes_count_distribution(res_by_gamma, 'gamma')
    plot_genes_occurrence_distribution(res_by_alpha, 'alpha')
    plot_genes_occurrence_distribution(res_by_beta, 'beta')
    plot_genes_occurrence_distribution(res_by_gamma, 'gamma')
    plot_top_genes_occurrences_table(res_by_alpha, 'alpha')
    plot_top_genes_occurrences_table(res_by_beta, 'beta')
    plot_top_genes_occurrences_table(res_by_gamma, 'gamma')
    
    
    
    
    
    
    