import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import bipartite


def draw_graph(graph, patients=None, save=False, name=None):
    graph = graph.copy()
    plt.figure(figsize=(62, 62))
    patient_pathway = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 0]
    SNVs = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 1]
    pathways = list(set([name[1] for name in list(patient_pathway)]))
    if patients is None:
        subgraph = graph
    else:
        # Get all nodes connected to the given patients
        subgraph = nx.Graph()
        subgraph.add_nodes_from(SNVs, bipartite=1)
        for patient in patients:
            for pathway in pathways:
                if (patient, pathway) in graph.nodes():
                    subgraph.add_node((patient, pathway), bipartite=0)
                    subgraph.add_weighted_edges_from(graph.edges((patient, pathway),data=True))

        if not nx.is_connected(subgraph):
            # Get the largest connected component
            largest_cc = max(nx.connected_components(subgraph), key=len)

            # Remove all nodes that are not in the largest connected component
            subgraph = subgraph.subgraph(largest_cc).copy()

        # Separate by group
    l = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 0]
    r = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 1]
    pos = {}

    # Update position for node from each group
    pos.update((node, (1, float(index) / len(l))) for index, node in enumerate(l))
    pos.update((node, (2, float(index) / len(r))) for index, node in enumerate(r))

    nx.draw(subgraph, pos=pos, with_labels=True)
    if save:
        plt.savefig(name)
    else:
        plt.show()