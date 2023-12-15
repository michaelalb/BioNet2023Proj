import json
from datetime import datetime
from typing import Union, Tuple

import numpy as np

from LinearProgrammingSolution.GraphHandlers import *
from LinearProgrammingSolution.MatchingDataHandler import MatchingDataHandler
from LinearProgrammingSolution.MatchingSolver import MatchingSolver
from LinearProgrammingSolution.MatchingVisualizer import *


def run_ilp_analysis(path_to_data: str,
                     should_draw_graph: bool,
                     alpha: float = 0.2,
                     beta: float = 0.2,
                     gamma: float = 1,
                     should_save_files: bool = True,
                     base_path: str = './'):
    """
    This function runs the entire ILP analysis.
    :param path_to_data: Path to load patient data into graph.
    :param should_draw_graph: Should draw the graph before and after optimization.
    :param alpha: alpha parameter for ILP - explained in optimization function.
    :param beta: beta parameter for ILP - explained in optimization function.
    :param gamma: gamma parameter for ILP - explained in optimization function.
    :param should_save_files: Should save the graph before and after optimization.
    :param base_path: Base path to save the files.
    :return:
    cover set - set of nodes that are in the min cover set.
    not cover set - set of nodes that are not in the min cover set.
    bottom cover set - set of nodes that are in the min cover set and are on the "bottom" of the graph - genes.
    top cover set - set of nodes that are in the min cover set and are on the "top" of the graph -
    patients_pathway pairs.
    new_graph - the graph after optimization.
    orig_graph - the original graph.
    adjusted_gene_weights - the weights of the edges after optimization.
    """
    matching_data_handler = MatchingDataHandler(path_to_data)
    print(f'loading data - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    matching_data_handler.load_data()
    print(f'creating graph - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    orig_graph = matching_data_handler.get_graph()
    matching_solver = MatchingSolver()
    print(f"saving graph picture - {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    if should_draw_graph:
        draw_graph(orig_graph, save=True, name='before.png')
    print(f"finding min cover set - {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, adjusted_gene_weights = \
        matching_solver.find_min_cover_set(orig_graph, alpha, beta, gamma,
                                           should_save_files=should_save_files,
                                           base_path=base_path)
    print(f"saving optimized graph picture - {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    if should_draw_graph:
        draw_graph(new_graph, save=True, name='after.png')
    return cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights


def get_sorted_genes_by_wight_from_graph(optimized_graph: nx.Graph = None, should_save_files: bool = True,
                                         base_path: str = './'):
    """
    This function returns a list of gene names sorted by their weight - extracted from the graph.
    :param optimized_graph:
    :param should_save_files:
    :param base_path:
    :return:
    """
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    print(f'loading graph - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    new_graph = load_graph_from_file('SingleRunResults/new_graph.pkl') if optimized_graph is None else optimized_graph
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    gene_weights = calculate_gene_edge_weight_sum(new_graph)
    sorted_gene_names_by_weight = sorted(gene_weights, key=lambda x: gene_weights.get(x), reverse=True)
    if should_save_files:
        with open(str(Path(base_path) / 'sorted_gene_names_by_weight.json'), 'w+') as f:
            json.dump(sorted_gene_names_by_weight, f, indent=4)
        with open(str(Path(base_path) / 'gene_weights.json'), 'w+') as f:
            json.dump(gene_weights, f, indent=4)
    return sorted_gene_names_by_weight, gene_weights


def get_sorted_genes_by_wight_from_dict(optimized_gene_weights: dict, should_save_files: bool = True,
                                        base_path: str = './') -> list:
    """
    This function returns a list of gene names sorted by their weight.
    :param optimized_gene_weights: dictionary with the gene names as keys and their weight as values.
    :param should_save_files: whether to save the results to a file.
    :param base_path: path to save the results to.
    :return: a list of gene names sorted by their weight.
    """
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    sorted_gene_names_by_weight = sorted(optimized_gene_weights,
                                         key=lambda x: optimized_gene_weights.get(x), reverse=True)
    if should_save_files:
        with open(str(Path(base_path) / 'sorted_gene_names_by_weight.json'), 'w+') as f:
            json.dump(sorted_gene_names_by_weight, f, indent=4)
    return sorted_gene_names_by_weight


def get_param_range(param_limits: dict, total_number_of_steps: int,
                    param_name: str = 'gene_param') -> Union[np.ndarray, list]:
    """
    This function returns the range of values to check for a given parameter.
    :param param_limits: dictionary with the limits for the parameters - usage explained in main.
    :param total_number_of_steps: maximum number of steps to take for each parameter.
    :param param_name: name of the parameter to get the range for - must exist in param_limits.
    :return: an iterable of values to check for the given parameter.
    """
    if param_limits[param_name].get('strict_vals') is not None:
        return param_limits[param_name]['strict_vals']
    param_range = np.arange(param_limits[param_name]['left_bound'],
                            param_limits[param_name]['right_bound'],
                            param_limits[param_name]['step_size'])
    if len(param_range) > total_number_of_steps:
        param_range = np.linspace(param_limits[param_name]['left_bound'],
                                  param_limits[param_name]['right_bound'],
                                  total_number_of_steps)
    return param_range


def set_up_param_ranges(param_limits: dict, total_number_of_steps: int) -> Tuple[Union[np.ndarray, list],
                                                                                 Union[np.ndarray, list],
                                                                                 Union[np.ndarray, list]]:
    """
    This function sets up the parameter ranges for the parameter search.
    :param param_limits: dictionary with the limits for the parameters - usage explained in main.
    :param total_number_of_steps: maximum number of steps to take for each parameter.
    :return:
    """
    alpha_param_range = get_param_range(param_limits, total_number_of_steps, 'alpha')
    beta_param_range = get_param_range(param_limits, total_number_of_steps, 'beta')
    gamma_param_range = get_param_range(param_limits, total_number_of_steps, 'gamma')
    return alpha_param_range, beta_param_range, gamma_param_range
