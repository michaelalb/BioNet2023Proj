import json
from datetime import datetime
from pathlib import Path

import numpy as np

from WeightedBiPartideGraphMatching.GraphHandlers import *
from WeightedBiPartideGraphMatching.MatchingDataHandler import MatchingDataHandler
from WeightedBiPartideGraphMatching.MatchingSolver import MatchingSolver
from WeightedBiPartideGraphMatching.MatchingVisualizer import *
from performance_evaluation import check_performances


def run_ilp_analysis(path_to_data: str,
                     should_draw_graph: bool,
                     alpha: float = 0.2,
                     should_save_files: bool = True,
                     base_path: str = './'):
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
    cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, adjusted_gene_weights =\
        matching_solver.find_min_cover_set(orig_graph, alpha,
                                           should_save_files=should_save_files,
                                           base_path=base_path)
    print(f"saving optimized graph picture - {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    if should_draw_graph:
        draw_graph(new_graph, save=True, name='after.png')
    return cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights


def get_sorted_genes_by_wight_from_graph(optimized_graph: nx.Graph = None, should_save_files: bool = True,
                                         base_path: str = './'):
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    print(f'loading graph - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    new_graph = load_graph_from_file('new_graph.pkl') if optimized_graph is None else optimized_graph
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
                                        base_path: str = './'):
    print(f'calculating gene weights - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
    sorted_gene_names_by_weight = sorted(optimized_gene_weights,
                                         key=lambda x: optimized_gene_weights.get(x), reverse=True)
    if should_save_files:
        with open(str(Path(base_path) / 'sorted_gene_names_by_weight.json'), 'w+') as f:
            json.dump(sorted_gene_names_by_weight, f, indent=4)
    return sorted_gene_names_by_weight


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

def get_param_range(param_limits: dict, total_number_of_steps: int,
                    param_name: str = 'gene_param'):
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


def set_up_param_ranges(param_limits: dict, total_number_of_steps: int):
    gene_param_range = get_param_range(param_limits, total_number_of_steps, 'gene_param')
    gene_penalty_patient_discount_range = get_param_range(param_limits, total_number_of_steps,
                                                          'gene_penalty_patient_discount')
    return gene_param_range, gene_penalty_patient_discount_range


def param_search(param_limits: dict,
                 gene_number_to_optimize: int = 5, total_number_of_steps: int = 50):
    steps_dict = {'search_results': {},
                  "gene_number_to_optimize": gene_number_to_optimize}
    gene_param_range, gene_penalty_patient_discount_range = set_up_param_ranges(param_limits, total_number_of_steps)
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    patient_snps = load_patient_snps()
    best_performance, best_performance_gene_param, best_performance_gene_penalty_patient_discount_param = 0, 0, 0
    # path set up
    base_run_path = Path(f'./ParamOptimizationResults/{datetime.now().strftime("%m_%d_%Y_%H_%M")}')

    for gene_param in gene_param_range:
        for gene_penalty_patient_discount_param in gene_penalty_patient_discount_range:
            current_run_path = Path(base_run_path / f'{gene_param=}_{gene_penalty_patient_discount_param=}')
            current_run_path.mkdir(parents=True, exist_ok=True)
            print(f'Optimizing {gene_param=} - {gene_penalty_patient_discount_param=}'
                  f' - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
            cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights = \
                run_ilp_analysis(path_to_data='Data/DriverMaxSetApproximation/BaseData',
                                 should_draw_graph=False,
                                 new_gene_penalty=gene_param,
                                 should_save_files=True,
                                 gene_penalty_patient_discount=gene_penalty_patient_discount_param,
                                 base_path=str(current_run_path))

            # sorted_gene_names_by_weight = get_sorted_genes_by_wight_from_dict(adjusted_gene_weights,
            #                                                                   should_save_files=True,
            #                                                                   base_path=str(current_run_path))
            # optimized_patient_genes = get_patient_genes_from_graph(new_graph)
            sorted_gene_names_by_weight, gene_weights = \
                get_sorted_genes_by_wight_from_graph(new_graph, should_save_files=True, base_path=str(current_run_path))

            ranked_genes_lists = get_rank_per_patient_from_base_data(sorted_gene_names_by_weight,
                                                                     patient_snps)
            with open(str(current_run_path / 'ranked_genes_lists.json'), 'w+') as f:
                json.dump(ranked_genes_lists, f, indent=4)
            our_performances = check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers)
            target_performance = our_performances['precision'][gene_number_to_optimize - 1]
            steps_dict['search_results'][str((gene_param, gene_penalty_patient_discount_param))] = {
                "gene_param": gene_param,
                "gene_penalty_patient_discount_param": gene_penalty_patient_discount_param,
                'precision': list(our_performances['precision']),
                'recall': list(our_performances['recall']),
                'f1': list(our_performances['f1']),
                'target_performance': target_performance,
            }
            target_performance = our_performances['precision'][gene_number_to_optimize - 1]
            print(f'performance for gene_param {gene_param} is {target_performance}')
            if target_performance > best_performance:
                best_performance = target_performance
                best_performance_gene_param = gene_param
                best_performance_gene_penalty_patient_discount_param = gene_penalty_patient_discount_param
            with open(str(base_run_path / 'param_search.json'), 'w+') as f:
                json.dump(steps_dict, f, indent=4)
    with open(str(base_run_path / 'param_search.json'), 'w+') as f:
        json.dump(steps_dict, f, indent=4)
    print(f'best performance of {best_performance} - with params {best_performance_gene_param=}'
          f' - {best_performance_gene_penalty_patient_discount_param=}')
    return steps_dict

def run_single_ilp_analysis(alpha: float, current_run_path: Path):
    patient_snps = load_patient_snps()
    cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights = \
        run_ilp_analysis(path_to_data='Data/DriverMaxSetApproximation/BaseData',
                            should_draw_graph=False,
                            alpha=alpha,
                            should_save_files=True,
                            base_path=str(current_run_path))

    # sorted_gene_names_by_weight = get_sorted_genes_by_wight_from_dict(adjusted_gene_weights,
    #                                                                   should_save_files=True,
    #                                                                   base_path=str(current_run_path))
    # optimized_patient_genes = get_patient_genes_from_graph(new_graph)
    sorted_gene_names_by_weight, gene_weights = \
        get_sorted_genes_by_wight_from_graph(new_graph, should_save_files=True, base_path=str(current_run_path))

    ranked_genes_lists = get_rank_per_patient_from_base_data(sorted_gene_names_by_weight,
                                                                patient_snps)
    with open(str(current_run_path / 'ranked_genes_lists.json'), 'w+') as f:
        json.dump(ranked_genes_lists, f, indent=4)

def main(should_calc_optimized_graph: bool = False, path_to_base_data: str = 'Data/DriverMaxSetApproximation/BaseData',
         should_draw_graph: bool = False, should_calc_gene_weights: bool = False, should_calc_sub_graph: bool = False,
         gene_weights: dict = None, new_graph: nx.Graph = None,
         k: int = 20, should_perform_param_search: bool = False):
    if should_calc_optimized_graph:
        cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights = \
            run_ilp_analysis(path_to_data=path_to_base_data,
                             should_draw_graph=should_draw_graph)
    if should_calc_gene_weights:
        sorted_gene_names_by_weight, gene_weights = get_rank_per_patient_from_base_data(new_graph)

    # if should_calc_sub_graph:
    #     # for sub graphing
    #     top_k_sub_graph, patients_with_ranked_genes_by_weight = get_graph_with_top_k_edges_from_graph_by_weight_sum(
    #         new_graph, k=k, gene_weights=gene_weights, should_draw_graph=False, should_save_files=False)

    if should_perform_param_search:
        # for param search
        param_limits = {
            'gene_param':
                {
                    'strict_vals': [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 10, 50],
                    'left_bound': 0.1,
                    'right_bound': 0.15,
                    'step_size': 0.05
                },
            "gene_penalty_patient_discount":
                {
                    'strict_vals': [1],
                    'left_bound': 0.1,
                    'right_bound': 0.15,
                    'step_size': 0.05
                }
        }
        param_search(param_limits=param_limits, gene_number_to_optimize=5, total_number_of_steps=50)


if __name__ == '__main__':
    # main(should_perform_param_search=True)
    run_single_ilp_analysis(alpha=0.2, current_run_path=Path('.'))
