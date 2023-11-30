import pickle

import networkx as nx


def create_new_bipartite_graph(top_nodes, bottom_nodes, edges):
    """
    Creates a new bipartide graph with the given number of nodes
    :param top_nodes: The number of nodes in the top set
    :param bottom_nodes: The number of nodes in the bottom setq
    :param edges: The edges of the graph
    :return: The created bipartide graph
    """
    # Create a new bipartite graph
    new_graph = nx.Graph()

    # Add nodes with the bipartite attribute
    new_graph.add_nodes_from(top_nodes, bipartite=0)
    new_graph.add_nodes_from(bottom_nodes, bipartite=1)

    # Add edges
    for i, edge in enumerate(edges):
        # no weight
        if len(edge) == 2:
            u, v = edge
            new_graph.add_edge(u, v)
        # weight
        else:
            u, v, data = edge
            new_graph.add_edge(u, v, weight=data['weight'])

    return new_graph


def load_graph_from_file(path_to_graph):
    """
    Loads a graph from a file
    :param path_to_graph: The path to the graph
    :return: The loaded graph
    """
    with open(path_to_graph, 'rb') as f:
        graph = pickle.load(f)
    return graph


def calculate_gene_edge_weight_sum(graph):
    # Get the nodes for the specified side
    nodes = [name for name, data in graph.nodes(data=True) if isinstance(name, str)]

    weight_sums = {}
    for node in nodes:
        # Get the edges for the current node
        edges = graph.edges(node, data=True)
        weight_sum = 0
        for _, _, data in edges:
            weight_sum += abs(data['weight'])
        weight_sums[node] = weight_sum

    return weight_sums


def get_graph_with_top_k_edges_from_graph_by_weight_sum(graph, node_weight_dict, k):
    # Sort the dictionary by weight in descending order and get the top k nodes
    top_k_nodes = sorted(node_weight_dict, key=node_weight_dict.get, reverse=True)[:k]

    # Create a subgraph with only the top k nodes and their associated edges
    subgraph = graph.subgraph(top_k_nodes)

    # Add the neighbors of the top k nodes to the subgraph
    for node in top_k_nodes:
        neighbors = graph.neighbors(node)
        for neighbor in neighbors:
            if neighbor not in subgraph:
                edges_to_add = graph.edges(neighbor)
                subgraph.add_edges_from(edges_to_add)

    return subgraph
