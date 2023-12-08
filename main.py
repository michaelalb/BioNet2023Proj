import json
from datetime import datetime

import numpy as np

from WeightedBiPartideGraphMatching.GraphHandlers import *
from WeightedBiPartideGraphMatching.MatchingDataHandler import MatchingDataHandler
from WeightedBiPartideGraphMatching.MatchingSolver import MatchingSolver
from WeightedBiPartideGraphMatching.MatchingVisualizer import *
from performance_evaluation import check_performances


def run_ilp_analysis(path_to_data: str,
                     should_draw_graph: bool,
                     new_gene_penalty: float = 0.2,
                     should_save_files: bool = True):
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
    cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph = matching_solver.find_min_cover_set(
        orig_graph, new_gene_penalty, should_save_files=should_save_files)
    print(f"saving optimized graph picture - {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    if should_draw_graph:
        draw_graph(new_graph, save=True, name='after.png')
    return cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph


def get_sorted_genes_by_wight(optimized_graph: nx.Graph = None, should_save_files: bool = True):
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    print(f'loading graph - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    new_graph = load_graph_from_file('new_graph.pkl') if optimized_graph is None else optimized_graph
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    gene_weights = calculate_gene_edge_weight_sum(new_graph)
    sorted_gene_names_by_weight = sorted(gene_weights, key=lambda x: gene_weights.get(x), reverse=True)
    if should_save_files:
        with open('./sorted_gene_names_by_weight.json', 'w+') as f:
            json.dump(sorted_gene_names_by_weight, f)
        with open('./gene_weights.json', 'w+') as f:
            json.dump(gene_weights, f)
    return sorted_gene_names_by_weight, gene_weights


# def get_graph_with_top_k_edges_from_graph_by_weight_sum(optimized_graph: nx.Graph, k: int,
#                                                         gene_weights: dict = None,
#                                                         should_draw_graph: bool = False,
#                                                         should_save_files: bool = True):
#     patients_with_ranked_genes_by_weight = get_rank_per_patient_from_base_data(optimized_graph=optimized_graph,
#                                                                                should_save_files=should_save_files)
#
#     # print(f'getting top {k} sub graph - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
#     # top_k_sub_graph = get_graph_with_top_k_edges_from_graph_by_weight_sum(optimized_graph, gene_weights, k)
#     # print(f'drawing graph - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
#     # if should_draw_graph:
#     #     draw_graph(top_k_sub_graph, save=True, name=f'top_{k}.png')
#     # if should_save_files:
#     #     with open('./patients_with_ranked_genes_by_weight.json', 'w+') as f:
#     #         json.dump(patients_with_ranked_genes_by_weight, f)
#     #     with open(f'./top_{k}_genes_graph.pkl', 'wb+') as f:
#     #         pickle.dump(top_k_sub_graph, f)
#     # print(f'done - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
#     # return top_k_sub_graph, patients_with_ranked_genes_by_weight
#     return patients_with_ranked_genes_by_weight


def param_search(gene_number_to_optimize: int = 5, total_number_of_steps: int = 50):
    steps_dict = {'gene_param': {},
                  "gene_number_to_optimize": gene_number_to_optimize}
    right_bound = 0.5
    left_bound = 0.1
    best_performance_param = 0.2
    best_performance = 0
    step_size = 0.05
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    patient_snps = load_patient_snps()

    for gene_param in np.arange(left_bound, right_bound, step_size):
        print(f'Optimizing gene_param: {gene_param}')
        cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph = \
            run_ilp_analysis(path_to_data='Data/DriverMaxSetApproximation/BaseData',
                             should_draw_graph=False,
                             new_gene_penalty=gene_param,
                             should_save_files=False)
        sorted_gene_names_by_weight, gene_weights = get_sorted_genes_by_wight(new_graph, should_save_files=False)
        ranked_genes_lists = get_rank_per_patient_from_base_data(sorted_gene_names_by_weight, patient_snps)
        our_performances = check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers)
        steps_dict['gene_param'][gene_param] = {
            'precision': our_performances['precision'],
            'recall': our_performances['recall'],
            'f1': our_performances['f1']
        }
        target_performance = our_performances['precision'][gene_number_to_optimize]
        print(f'performance for gene_param {gene_param} is {target_performance}')
        if target_performance > best_performance:
            best_performance = target_performance
            best_performance_param = gene_param
    with open('./Data/DriverMaxSetApproximation/param_search.json', 'w+') as f:
        json.dump(steps_dict, f)
    print(f'best param is {best_performance_param} with performance of {best_performance}')
    return steps_dict


def main(should_calc_optimized_graph: bool = False, path_to_base_data: str = 'Data/DriverMaxSetApproximation/BaseData',
         should_draw_graph: bool = False, should_calc_gene_weights: bool = False, should_calc_sub_graph: bool = False,
         gene_weights: dict = None, new_graph: nx.Graph = None,
         k: int = 20, should_perform_param_search: bool = False):
    if should_calc_optimized_graph:
        cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph = \
            run_ilp_analysis(path_to_data=path_to_base_data,
                             should_draw_graph=should_draw_graph)
    if should_calc_gene_weights:
        sorted_gene_names_by_weight, gene_weights = get_rank_per_patient_from_base_data(new_graph)

    if should_calc_sub_graph:
        # for sub graphing
        top_k_sub_graph, patients_with_ranked_genes_by_weight = get_graph_with_top_k_edges_from_graph_by_weight_sum(
            new_graph, k=k, gene_weights=gene_weights, should_draw_graph=False, should_save_files=False)

    if should_perform_param_search:
        # for param search
        param_search()


if __name__ == '__main__':
    main(should_perform_param_search=True)
