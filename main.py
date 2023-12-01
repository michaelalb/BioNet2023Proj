import json

from WeightedBiPartideGraphMatching.GraphHandlers import *
from WeightedBiPartideGraphMatching.MatchingDataHandler import MatchingDataHandler
from WeightedBiPartideGraphMatching.MatchingSolver import MatchingSolver
from WeightedBiPartideGraphMatching.MatchingVisualizer import *

SHOULD_CALC_NEW_GRAPH = False

if __name__ == '__main__':
    if SHOULD_CALC_NEW_GRAPH:
        path_to_data = 'Data/DriverMaxSetApproximation/BaseData'
        matching_data_handler = MatchingDataHandler(path_to_data)
        print('loading data')
        matching_data_handler.load_data()
        print('creating graph')
        graph = matching_data_handler.get_graph()
        matching_solver = MatchingSolver()
        print("saving graph picture")
        draw_graph(graph, save=True, name='before.png')
        print("finding min cover set")
        cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph = matching_solver.find_min_cover_set(graph)
        print("saving optimized graph picture")
        draw_graph(new_graph, save=True, name='after.png')
    else:
        print('loading graph')
        new_graph = load_graph_from_file('new_graph.pkl')
        print('calculating gene weights')
        gene_weights = calculate_gene_edge_weight_sum(new_graph)
        sorted_gene_names_by_weight = sorted(gene_weights, key=lambda x: gene_weights.get(x), reverse=True)

        with open('./sorted_gene_names_by_weight.json', 'w+') as f:
            json.dump(sorted_gene_names_by_weight, f)
        with open('./gene_weights.json', 'w+') as f:
            json.dump(gene_weights, f)

        # for sub graphing
        patients_with_ranked_genes_by_weight = get_rank_per_patient(new_graph, sorted_gene_names_by_weight)
        with open('./patients_with_ranked_genes_by_weight.json', 'w+') as f:
            json.dump(patients_with_ranked_genes_by_weight, f)
        # k = 20
        # print('getting top k sub graph')
        # top_k_sub_graph = get_graph_with_top_k_edges_from_graph_by_weight_sum(new_graph, gene_weights, k)
        # print('drawing graph')
        # draw_graph(top_k_sub_graph, save=True, name=f'top_{k}.png')
        # with open(f'./top_{k}_genes_graph.pkl', 'wb+') as f:
        #     pickle.dump(top_k_sub_graph, f)
        # print('a')