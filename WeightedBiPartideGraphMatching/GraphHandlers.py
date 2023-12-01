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
    nei = [i for node in top_k_nodes for i in graph.neighbors(node)]

    # Create a subgraph with only the top k nodes and their associated edges
    subgraph = graph.subgraph(top_k_nodes + nei).copy()
    print(f"number of nodes in original graph: {len(graph.nodes())} reduced to {len(subgraph.nodes())}")
    print(f"number of edges in original graph: {len(graph.edges())} reduced to {len(subgraph.edges())}")
    return subgraph


def get_rank_per_patient(graph, sorted_genes):
    patient_node_names = [name for name, data in graph.nodes(data=True) if isinstance(name, tuple)]
    patient_names = list(set([name[0] if name[0].find('.') != -1 else name[1] for name in patient_node_names]))

    patient_genes = {}
    # get all patient genes
    for patient_name in patient_names:
        patient_genes[patient_name] = []
        for node in graph.nodes():
            if isinstance(node, tuple) and (node[0] == patient_name or node[1] == patient_name):
                gene = [i for i in graph.neighbors(node)]
                assert len(gene) == 1
                patient_genes[patient_name].append(gene[0])
        patient_genes[patient_name] = list(set(patient_genes[patient_name]))
    sorted_patient_genes = {patient: sort_list_by_reference(sorted_genes, genes) for
                            patient, genes in patient_genes.items()}
    return sorted_patient_genes


def sort_list_by_reference(reference_list, unordered_list):
    return sorted(unordered_list, key=reference_list.index)
