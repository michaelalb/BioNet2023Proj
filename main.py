from WeightedBiPartideGraphMatching.MatchingDataHandler import MatchingDataHandler
from WeightedBiPartideGraphMatching.MatchingSolver import MatchingSolver
from WeightedBiPartideGraphMatching.MatchingVisualizer import *

if __name__ == '__main__':
    path_to_data = 'Data/DriverMaxSetApproximation/BaseData'
    matching_data_handler = MatchingDataHandler(path_to_data)
    matching_data_handler.load_data()
    graph = matching_data_handler.get_graph()
    matching_solver = MatchingSolver()
    # draw_graph(graph, patients=['TCGA.3L.AA1B.01'])

    cover_set = matching_solver.find_min_cover_set(graph)

    print('nodes')
    print(graph.nodes(data=True))
    print('edges')
    print(graph.edges(data=True))