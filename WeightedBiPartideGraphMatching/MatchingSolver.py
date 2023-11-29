from pulp import LpProblem, LpVariable, lpSum, LpInteger, LpMaximize
import json
import networkx as nx
import pickle


class MatchingSolver:

    def _custom_sort(self, t):
        "sort so gene  nodes are first, then pathway nodes"
        # todo: assert this actually does this
        if isinstance(t[0], str):
            return t
        else:
            if len(t) == 2:
                return t[1], t[0]
            else:
                return t[1], t[0], t[2]


    def find_min_cover_set(self, graph: nx.Graph):
        prob = LpProblem("Maximum_Weight_Cover_Set", LpMaximize)

        # Create a binary variable to state that a node is included in the cover
        sorted_edges = [self._custom_sort(i) for i in graph.edges()]
        sorted_edges_data = [self._custom_sort(i) for i in graph.edges(data=True)]

        top_nodes = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 0]
        bottom_nodes = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 1]
        patient_pathway_nodes = LpVariable.dicts("x", top_nodes, 0, 1, LpInteger)
        gene_nodes = LpVariable.dicts("y", bottom_nodes, 0, 1, LpInteger)
        edges = LpVariable.dicts("e", sorted_edges, 0, 1, LpInteger)

        # Objective function
        prob += lpSum([edges[(n, n1)] * abs(d['weight']) for n, n1, d in sorted_edges_data])
        # Constraints
        # Constraint - 1 : Edge-Vertex relationship
        for (gene_node, pathway_node) in sorted_edges:
            # goes over edges so that pathway nodes are first, then gene nodes
            # if we choose an edge both bottom and top nodes must be chosen
            prob += edges[(gene_node, pathway_node)] <= patient_pathway_nodes[pathway_node]
            prob += edges[(gene_node, pathway_node)] <= gene_nodes[gene_node]

        # Constraint - 2 : for each pathway node only one gene node can be chosen
        for pathway_node in top_nodes:
            prob += lpSum([edges[(gene_node, pathway_node)] for gene_node in graph.neighbors(pathway_node)]) <= 1

        # Constraint - 3 : dont allow disconnected nodes
        for node in top_nodes:
            node_edges = [edge for edge in sorted_edges if node in edge]
            prob += lpSum([edges[edge] for edge in node_edges]) >= patient_pathway_nodes[node]
        for node in bottom_nodes:
            node_edges = [edge for edge in sorted_edges if node in edge]
            prob += lpSum([edges[edge] for edge in node_edges]) >= gene_nodes[node]
        prob.solve()

        # Return the nodes included in the cover
        bottom_cover_set = [node for node in bottom_nodes if gene_nodes[node].varValue == 1]
        top_cover_set = [node for node in top_nodes if patient_pathway_nodes[node].varValue == 1]
        cover_set = bottom_cover_set + top_cover_set
        not_cover_set = [node for node in graph.nodes() if node not in cover_set]
        chosen_edges = [edge for edge in sorted_edges if edges[edge].varValue == 1]

        # print optimization results
        print(f"choose {len(bottom_cover_set)} genes out of {len(bottom_nodes)} and {len(top_cover_set)} pathways out of {len(top_nodes)}")
        print(f"total weight: {prob.objective.value()}")
        print(f'chosen edges: {len(chosen_edges)} out of {len(sorted_edges)}')
        print(f'total sum of weights: {sum([abs(d["weight"]) for _, _, d in sorted_edges_data])}')

        # save results to json file
        data = {
            'bottom_cover_set': bottom_cover_set,
            'top_cover_set': top_cover_set,
            'chosen_edges': chosen_edges,
            'not_cover_set': not_cover_set,
            'total_weight': prob.objective.value(),
            'total_sum_of_weights': sum([abs(d["weight"]) for _, _, d  in sorted_edges_data]),
        }
        with open('./cover_set.json', 'w+') as f:
            json.dump(data, f)

        # Create a new bipartite graph
        new_graph = nx.Graph()

        # Add nodes with the bipartite attribute
        new_graph.add_nodes_from(top_cover_set, bipartite=0)
        new_graph.add_nodes_from(bottom_cover_set, bipartite=1)

        # Add edges
        for edge in chosen_edges:
            u, v = edge
            new_graph.add_edge(u, v)#, weight=weight)

        with open('./new_graph.pickle.pkl', 'wb+') as f:
            pickle.dump(new_graph, f)

        c = {}
        for gene_node in bottom_cover_set:
            deg = new_graph.degree(gene_node)
            c[deg] = c.get(deg, 0) + 1
        print(c)
        # todo:
        # 1. add weights to new graph
        # 2. build out to functions
        # 3. filter out weakly connected genes

        return cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph