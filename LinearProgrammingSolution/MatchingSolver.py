from datetime import datetime
from pathlib import Path

from pulp import LpProblem, LpVariable, lpSum, LpInteger, LpMaximize
import json
import networkx as nx
import pickle

from LinearProgrammingSolution.GraphHandlers import create_new_bipartite_graph


class MatchingSolver:
    """
    This class solves the ILP problem of finding the maximum weighted matching in a bipartite graph.
    """

    def _custom_sort(self, t):
        "sort so gene  nodes are first, then pathway nodes"
        if isinstance(t[0], str):
            return t
        else:
            if len(t) == 2:
                return t[1], t[0]
            else:
                return t[1], t[0], t[2]

    def get_nodes_edges_from_ILP_res(self, graph: nx.Graph,
                                     sorted_edges_data,
                                     edges,
                                     patient_pathway_nodes,
                                     gene_nodes,
                                     bottom_nodes,
                                     top_nodes):
        """
        Get the nodes and edges from the ILP result
        :param graph: The original graph
        :param sorted_edges_data: sorted edges data with weight
        :param edges: ILP dict of edges with results
        :param patient_pathway_nodes: ILP dict of  all original top nodes with results
        :param gene_nodes: ILP dict of all original bottom nodes with results
        :param bottom_nodes:list of all original bottom nodes with results - genes
        :param top_nodes: list of all original top nodes with results - patient pathways
        :return:
        """

        bottom_cover_set = [node for node in bottom_nodes if gene_nodes[node].varValue == 1]
        top_cover_set = [node for node in top_nodes if patient_pathway_nodes[node].varValue == 1]
        cover_set = bottom_cover_set + top_cover_set
        not_cover_set = [node for node in graph.nodes() if node not in cover_set]
        chosen_edges = [edge for edge in sorted_edges_data if edges[edge[:2]].varValue == 1]
        return bottom_cover_set, top_cover_set, cover_set, not_cover_set, chosen_edges

    def print_and_save_ILP_res(self, bottom_cover_set, bottom_nodes, top_cover_set, top_nodes, prob,
                               sorted_edges, sorted_edges_data, chosen_edges, not_cover_set,
                               new_graph, gene_weights, gene_adjustments, gene_wight_multipliers,
                               should_save_files=True, base_path='./'):
        """
        a wrapper function for printing and saving the ILP results.
        """
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
            'adjusted_gene_weights': gene_weights,
            'gene_adjustments': gene_adjustments,
            'gene_wight_multipliers': gene_wight_multipliers,
            'total_sum_of_original_weights': sum([abs(d["weight"]) for _, _, d in sorted_edges_data]),
        }
        if should_save_files:
            with open(str(Path(base_path) / f'cover_set.json'), 'w+') as f:
                json.dump(data, f, indent=4)

            with open(str(Path(base_path) / 'new_graph.pkl'), 'wb+') as f:
                pickle.dump(new_graph, f)

    def find_min_cover_set(self, graph: nx.Graph, alpha, beta, gamma, should_save_files=True, base_path='./'):
        """
        This function solves the ILP problem of finding the maximum weighted matching in a bipartite graph.
        :param graph: The original graph.
        :param alpha: a parameter to control the globality of driver genes you wish to optimize for.
        It does so by adjusting the weights of the edges by the fraction of distinct patients affected by the gene
        across the entire cohort.
        :param beta: a parameter to control the locality of driver genes you wish to optimize for.
        It does so by adjusting the weights of the edges by the fraction of distinct pathways within a sepcific patient
        affected by the gene.
        :param gamma: a parameter to control the biological assumption you with to optimize.
        It sets the final amount of genes allowed to affect each pathway in the optimized graph.
        :param should_save_files: a boolean parameter to control whether to save the results to a file.
        :param base_path: the path to save the results to.
        :return:
        cover set - set of nodes that are in the min cover set.
        not cover set - set of nodes that are not in the min cover set.
        bottom cover set - set of nodes that are in the min cover set and are on the "bottom" of the graph - genes.
        top cover set - set of nodes that are in the min cover set and are on the "top" of the graph -
        patients_pathway pairs.
        new_graph - the graph after optimization.
        orig_graph - the original graph.
        adjusted_gene_weights - the weights of the edges after optimization.
        """
        print(f"Starting ILP with {alpha=}, {beta=}, {gamma=}")
        prob = LpProblem("Maximum_Weight_Cover_Set", LpMaximize)
        
        # Create a binary variable to state that a node is included in the cover
        sorted_edges = [self._custom_sort(i) for i in graph.edges()]
        sorted_edges_data = [self._custom_sort(i) for i in graph.edges(data=True)]

        top_nodes = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 0]
        bottom_nodes = [name for name, data in graph.nodes(data=True) if data['bipartite'] == 1]
        print("gene_patient_pairs generation")
        gene_patient_pairs = list(set([(gene, patient_name) for (gene, (patient_name, pathway)) in sorted_edges]))
        print("patient_pathway_nodes generation")
        patient_pathway_nodes = LpVariable.dicts("x", top_nodes, 0, 1, LpInteger)
        gene_nodes = LpVariable.dicts("y", bottom_nodes, 0, 1, LpInteger)
        # gene_patient_nodes = LpVariable.dicts("z", gene_patient_pairs, 0, 1, LpInteger)
        edges = LpVariable.dicts("e", sorted_edges, 0, 1, LpInteger)
        number_of_patients = len(set([patient_name for (patient_name, pathway) in top_nodes]))
        print("Done creating ILP vars")
        # Objective function
        print("Creating objective function")
        print("edge_wights generation")
        gene_patient_pathway_weight_counts, gene_patient_pathway_weight_multipliers = {}, {}
        patient_pathway_count = {}
        # Constraints
        # Constraint - 1 : Adjust weights of edges according to the number of patients covered
        print("Constraint - 1 : Adjust weights of edges according to the number of patients covered")
        for gene, (patient_name, pathway) in sorted_edges:
            gene_patient = (gene, patient_name)
            gene_patient_pathway_weight_counts[gene_patient] = gene_patient_pathway_weight_counts.get(
                gene_patient, 0) + 1
            if patient_name not in patient_pathway_count:
                patient_pathway_count[patient_name] = []
            patient_pathway_count[patient_name].append(pathway)
        for gene, patient_name in gene_patient_pathway_weight_counts:
            patient_pathways = len(set(patient_pathway_count[patient_name]))
            gene_patient_pathway_weight_multipliers[(gene, patient_name)] = 1 + (
                    beta * gene_patient_pathway_weight_counts[(gene, patient_name)] / patient_pathways)

        gene_wight_multipliers = {}
        # Constraints
        # Constraint - 2 : Adjust weights of edges according to the number of patients covered
        print("Constraint - 2 : Adjust weights of edges according to the number of pathways in patient covered")
        for gene in gene_nodes:
            covered_patients = len([patient_name for (gene1, patient_name) in
                                    gene_patient_pairs if gene == gene1])
            gene_wight_multipliers[gene] = 1 + alpha * covered_patients / number_of_patients
        prob += lpSum([edges[(gene, pathway)] * abs(d['weight']) * gene_wight_multipliers[gene] * gene_patient_pathway_weight_multipliers[(gene, pathway[0])]
                       for gene, pathway, d in sorted_edges_data])

        # Constraint - 3 : Edge-Vertex relationship
        print("Constraint - 3 : Edge-Vertex relationship")
        for (gene_node, pathway_node) in sorted_edges:
            # goes over edges so that pathway nodes are first, then gene nodes
            # if we choose an edge both bottom and top nodes must be chosen
            prob += edges[(gene_node, pathway_node)] <= patient_pathway_nodes[pathway_node]
            prob += edges[(gene_node, pathway_node)] <= gene_nodes[gene_node]

        # Constraint - 4 : for each pathway node only 'gamma' gene node can be chosen
        print("Constraint - 4 : for each pathway node only 'gamma' gene node can be chosen")
        for pathway_node in top_nodes:
            prob += lpSum([edges[(gene_node, pathway_node)] for gene_node in graph.neighbors(pathway_node)]) <= gamma

        # Constraint - 5 : dont allow disconnected nodes
        print("Constraint - 5 : dont allow disconnected nodes")
        print("Creating constraint for top nodes")
        for node in top_nodes:
            node_edges = [edge for edge in sorted_edges if node in edge]
            prob += lpSum([edges[edge] for edge in node_edges]) >= patient_pathway_nodes[node]
        print("Creating constraint for bottom nodes")
        for node in bottom_nodes:
            node_edges = [edge for edge in sorted_edges if node in edge]
            prob += lpSum([edges[edge] for edge in node_edges]) >= gene_nodes[node]
        print(f"Done creating ILP constraints - {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        print("solving")
        prob.solve()

        # Return the nodes included in the cover
        bottom_cover_set, top_cover_set, cover_set, not_cover_set, chosen_edges = self.get_nodes_edges_from_ILP_res(
            graph, sorted_edges_data, edges, patient_pathway_nodes, gene_nodes, bottom_nodes, top_nodes)

        # get adjusted gene weights
        gene_weights, gene_adjustments = self._get_gene_weights_by_penalty(
            chosen_edges, alpha, number_of_patients)

        # create new graph with only the chosen edges
        new_graph = create_new_bipartite_graph(bottom_cover_set, top_cover_set, chosen_edges)

        self.print_and_save_ILP_res(bottom_cover_set, bottom_nodes, top_cover_set, top_nodes, prob,
                                    sorted_edges, sorted_edges_data, chosen_edges, not_cover_set, new_graph,
                                    gene_weights, gene_adjustments, gene_wight_multipliers,
                                    should_save_files=should_save_files, base_path=base_path)

        return cover_set, not_cover_set, bottom_cover_set, top_cover_set, new_graph, gene_weights

    def _get_gene_weights_by_penalty(self, chosen_edges: list, alpha: float, number_of_patients):
        gene_weights = {}
        gene_patient_counts = {}
        adjusted_gene_weights = {}
        gene_adjustments = {}
        for edge in chosen_edges:
            if isinstance(edge[0], str):
                gene = edge[0]
                patient = edge[1][0]
            else:
                gene = edge[1]
                patient = edge[0][0]
            if gene not in gene_weights:
                gene_weights[gene] = 0
                gene_patient_counts[gene] = []
            gene_weights[gene] += edge[2]['weight']
            gene_patient_counts[gene].append(patient)

        for gene in gene_weights:
            number_of_covered_patients = len(set(gene_patient_counts[gene]))
            current_gene_wight_multiplier = 1 + alpha * number_of_covered_patients / number_of_patients
            adjusted_gene_weights[gene] = gene_weights.get(gene) * current_gene_wight_multiplier
            gene_adjustments[gene] = current_gene_wight_multiplier
        return adjusted_gene_weights, gene_adjustments
