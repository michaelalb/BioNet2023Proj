from pulp import LpProblem, LpMinimize, LpVariable, lpSum
import networkx as nx


class MatchingSolver:

    def find_min_cover_set(self, graph):
        prob = LpProblem("Minimum_Weight_Cover_Set", LpMinimize)

        # Create a binary variable to state that a node is included in the cover
        x = LpVariable.dicts("x", graph.nodes(), 0, 1, cat='Integer')

        # Objective function
        prob += lpSum([x[node] * weight for node, weight in nx.get_node_attributes(graph, 'weight').items()])

        # Constraints
        for edge in graph.edges():
            prob += x[edge[0]] + x[edge[1]] >= 1

        prob.solve()

        # Return the nodes included in the cover
        cover_set = [node for node in x if x[node].varValue == 1]
        return cover_set
