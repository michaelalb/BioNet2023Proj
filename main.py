import json
from datetime import datetime
from pathlib import Path

import networkx as nx

from LinearProgrammingSolution.GraphHandlers import load_patient_snps, get_patient_genes_from_graph, \
    get_rank_per_patient_from_base_data
from Utils import set_up_param_ranges, run_ilp_analysis, get_sorted_genes_by_wight_from_dict
from ResultsAnalysis.performance_evaluation import check_performances, plot_performances, performance_evaluation_main
from ResultsAnalysis.parameters_analysis import parameters_analysis_main


def param_search(param_limits: dict,
                 gene_number_to_optimize: int = 5, total_number_of_steps: int = 50) -> dict:
    """
    This function performs a parameter search for the best alpha, beta and gamma parameters for the algorithm.
    :param param_limits: dictionary with the limits for the parameters - usage explained in main.
    :param gene_number_to_optimize: number of genes to optimize for - looking at precision of the n-th gene.
    :param total_number_of_steps: maximum number of steps to perform for each parameter.
    :return: dictionary with results of the parameter search.
    """
    steps_dict = {'search_results': {},
                  "gene_number_to_optimize": gene_number_to_optimize}
    alpha_param_range, beta_param_range, gamma_param_range = set_up_param_ranges(param_limits, total_number_of_steps)
    gold_standard_drivers = json.load(open('./Data/gold_standard_drivers.json'))
    PRODIGY_results = json.load(open('./Data/PRODIGY_results.json'))
    patient_snps = load_patient_snps()
    PRODIGY_performances = check_performances(PRODIGY_results, patient_snps, gold_standard_drivers, )
    best_performance, best_performance_alpha_param, best_performance_beta_param, best_performance_gamma_param = 0, 0, 0, 0
    # path set up
    base_run_path = Path(f'./ParamOptimizationResults/{datetime.now().strftime("%m_%d_%Y_%H_%M")}')

    for alpha_param in alpha_param_range:
        for beta_param in beta_param_range:
            for gamma_param in gamma_param_range:
                current_run_path = Path(base_run_path / f'{alpha_param=}_{beta_param=}_{gamma_param=}')
                current_run_path.mkdir(parents=True, exist_ok=True)
                print(
                    f'Optimizing {alpha_param=} {beta_param=} {gamma_param=} - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
                cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights = \
                    run_ilp_analysis(path_to_data='Data/DriverMaxSetApproximation/BaseData',
                                     should_draw_graph=False,
                                     should_save_files=True,
                                     alpha=alpha_param,
                                     beta=beta_param,
                                     gamma=gamma_param,
                                     base_path=str(current_run_path))

                sorted_gene_names_by_weight = get_sorted_genes_by_wight_from_dict(adjusted_gene_weights,
                                                                                  should_save_files=True,
                                                                                  base_path=str(current_run_path))
                optimized_patient_genes = get_patient_genes_from_graph(new_graph)
                ranked_genes_lists = get_rank_per_patient_from_base_data(sorted_gene_names_by_weight,
                                                                         optimized_patient_genes)
                with open(str(current_run_path / 'ranked_genes_lists.json'), 'w+') as f:
                    json.dump(ranked_genes_lists, f, indent=4)
                our_performances = check_performances(ranked_genes_lists, patient_snps, gold_standard_drivers)
                target_performance = our_performances['precision'][gene_number_to_optimize - 1]
                steps_dict['search_results'][str((alpha_param, beta_param, gamma_param))] = {
                    "alpha_param": alpha_param,
                    "beta_param": beta_param,
                    "gamma_param": gamma_param,
                    'precision': list(our_performances['precision']),
                    'recall': list(our_performances['recall']),
                    'f1': list(our_performances['f1']),
                    'target_performance': target_performance,
                }
                target_performance = our_performances['precision'][gene_number_to_optimize - 1]
                print(f'performance for {alpha_param=}, {beta_param=} {gamma_param=} is {target_performance}')
                if target_performance > best_performance:
                    best_performance = target_performance
                    best_performance_alpha_param = alpha_param
                    best_performance_beta_param = beta_param
                    best_performance_gamma_param = gamma_param
                with open(str(base_run_path / 'param_search.json'), 'w+') as f:
                    json.dump(steps_dict, f, indent=4)
                print(f'plotting perf for {alpha_param=} - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}')
                plot_performances(
                    {'our algotithem': our_performances, 'PRODIGY': PRODIGY_performances},
                    save_path=str(current_run_path / 'performances.png'))
        perf_dict = {k: {
            "precision": v["precision"],
            "recall": v["recall"],
            "f1": v["f1"]} for k, v in steps_dict['search_results'].items()
            if str(k).find(str(alpha_param)) != -1 and str(k).find(str(alpha_param)) < str(k).find(',')}
        perf_dict['PRODIGY'] = PRODIGY_performances
        plot_performances(perf_dict, save_path=str(base_run_path / f'all_performances_{alpha_param=}.png'))
    with open(str(base_run_path / 'param_search.json'), 'w+') as f:
        json.dump(steps_dict, f, indent=4)
    perf_dict = {k: {
        "precision": v["precision"],
        "recall": v["recall"],
        "f1": v["f1"]} for k, v in steps_dict['search_results'].items()}
    perf_dict['PRODIGY'] = PRODIGY_performances
    plot_performances(perf_dict, save_path=str(base_run_path / 'all_performances.png'))
    print(
        f'best performance of {best_performance} - with params {best_performance_alpha_param=} {best_performance_beta_param=} {best_performance_gamma_param=}')
    return steps_dict


def run_single_ilp_analysis(alpha: float, beta: float, gamma: float, current_run_path: Path):
    """
    runs a single ILP analysis with the given parameters.
    :param alpha: controls the weight multiplier of edges by gene globality.
    :param beta: controls the weight multiplier of edges by gene locality.
    :param gamma: controls the number of genes allowed to be assigned to each patient_pathway.
    :param current_run_path: the path to save the results to.
    :return: None the results are saved to the given path.
    """
    cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights = \
        run_ilp_analysis(path_to_data='Data/DriverMaxSetApproximation/BaseData',
                         should_draw_graph=False,
                         alpha=alpha,
                         beta=beta,
                         gamma=gamma,
                         should_save_files=True,
                         base_path=str(current_run_path))

    sorted_gene_names_by_weight = get_sorted_genes_by_wight_from_dict(adjusted_gene_weights,
                                                                      should_save_files=True,
                                                                      base_path=str(current_run_path))
    optimized_patient_genes = get_patient_genes_from_graph(new_graph)
    ranked_genes_lists = get_rank_per_patient_from_base_data(sorted_gene_names_by_weight,
                                                             optimized_patient_genes)
    with open(str(current_run_path / 'ranked_genes_lists.json'), 'w+') as f:
        json.dump(ranked_genes_lists, f, indent=4)


def optimization_main(should_calc_optimized_graph: bool = False,
                      path_to_base_data: str = 'Data/DriverMaxSetApproximation/BaseData',
                      should_draw_graph: bool = False, should_calc_gene_weights: bool = False,
                      new_graph: nx.Graph = None, should_perform_param_search: bool = False,
                      param_limits: dict = None) -> None:
    """
    This is the main orchestrator of the project. It is responsible for running the different parts of the project.
    :param should_calc_optimized_graph: This flag indicates whether the optimized graph should be calculated from
    scratch - This parameter should be used if you'd like to simply run the ILP analysis only.
    :param should_calc_gene_weights: This flag indicates whether we should use the new graph created / provided by the
    ILP to get ranks per patient. Can either be used with should_calc_optimized_graph or with a provided graph using the
    new_graph param.
    :param should_perform_param_search: This flag indicates whether we should perform a parameter search for the ILP.
    Must be used with param_limits.
    :param should_draw_graph: This flag indicates whether we should draw the graphs, original and optimized.
    :param path_to_base_data: This is the path to the base data folder, it is used to save the files from all steps.
    :param new_graph: This is a graph that is provided to the function, it is used to get ranks per patient - can be
    used with should_calc_gene_weights.
    :param param_limits: This is a dictionary that contains the limits for the parameter search. Must be used with
    should_perform_param_search.
    :return: None - the communication with the function output is done through the files that are saved.
    """
    if should_calc_optimized_graph:
        cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, orig_graph, adjusted_gene_weights = \
            run_ilp_analysis(path_to_data=path_to_base_data,
                             should_draw_graph=should_draw_graph)
    if should_calc_gene_weights:
        sorted_gene_names_by_weight = get_rank_per_patient_from_base_data(new_graph)

    if should_perform_param_search:
        assert param_limits is not None, 'param_limits must be provided for param search'
        param_search(param_limits=param_limits, gene_number_to_optimize=5, total_number_of_steps=50)


if __name__ == '__main__':
    # for param search - alpha, beta, gamma
    # strict vals - values to check for each param - overrides left_bound, right_bound, step_size
    # left_bound - left bound for param search
    # right_bound - right bound for param search
    # step_size - step size for param search
    param_limits = {
        'alpha':
            {
                'strict_vals': [25, 0],
                'left_bound': 0.1,
                'right_bound': 0.15,
                'step_size': 0.05
            },
        'beta':
            {
                'strict_vals': [0.01, 0.05, 0.1, 0.2, 0.5, 2, 5, 10, 25, 0],
                'left_bound': 0.1,
                'right_bound': 0.15,
                'step_size': 0.05
            },
        'gamma':
            {
                'strict_vals': [1]
            }
    }
    performance_evaluation_main()
    # parameters_analysis_main()
    # optimization_main(should_perform_param_search=True)
    # run_single_ilp_analysis(alpha=0.2, beta=0.2, current_run_path=Path('./SingleRunResults/'))
    optimization_main(should_perform_param_search=True, param_limits=param_limits)
    # run_single_ilp_analysis(alpha=0.01, beta=0.01,gamma=20, current_run_path=Path('./SingleRunResults/'))
