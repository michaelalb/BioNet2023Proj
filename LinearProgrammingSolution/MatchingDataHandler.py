import pandas as pd
import networkx as nx
from pathlib import Path


class MatchingDataHandler:
    """
    This class is responsible for loading the data from the csv files and creating the graph.
    each patient-pathway is a top-node(0) in the graph, and each snv is a bottom-node(1) in the graph.
    The edges are weighted by the cell value in the csv file.
    zero values are not added to the graph.
    """
    def __init__(self, directory):
        self.directory = Path(directory)
        self._graph = nx.Graph()

    def load_data(self):
        for file in self.directory.glob('*.csv'):
            patient_name = file.stem
            df = pd.read_csv(file, index_col=0)
            self._add_to_graph(df, patient_name)

    def _add_to_graph(self, df, patient_name):
        for pathway in df.index:
            patient_pathway_node = (patient_name, pathway)
            self._graph.add_node(patient_pathway_node, bipartite=0)  # Add the patient-pathway node to one side

            for snv, weight in df.loc[pathway].items():
                if weight != 0:
                    self._graph.add_node(snv, bipartite=1)  # Add the snv to the other side
                    self._graph.add_edge(patient_pathway_node, snv,
                                         weight=weight)  # Add the edge weighted by the cell value

    def get_graph(self):
        return self._graph
