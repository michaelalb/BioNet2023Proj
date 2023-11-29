from WeightedBiPartideGraphMatching.MatchingDataHandler import MatchingDataHandler
from WeightedBiPartideGraphMatching.MatchingSolver import MatchingSolver
from WeightedBiPartideGraphMatching.MatchingVisualizer import *

if __name__ == '__main__':
    path_to_data = 'Data/DriverMaxSetApproximation/BaseData'
    matching_data_handler = MatchingDataHandler(path_to_data)
    matching_data_handler.load_data()
    graph = matching_data_handler.get_graph()
    matching_solver = MatchingSolver()
    draw_graph(graph, save=True, name='before.png')

    cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph = matching_solver.find_min_cover_set(graph)
    draw_graph(new_graph, save=True, name='after.png')