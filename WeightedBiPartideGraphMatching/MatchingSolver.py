from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpInteger


class MatchingSolver:

    def find_min_cover_set(self, graph, gamma=0.8, alpha=1, beta=0.4):
        prob = LpProblem("Minimum_Weight_Cover_Set", LpMinimize)

        # Create a binary variable to state that a node is included in the cover
        x = LpVariable.dicts("x", graph.nodes(), 0, 1, LpInteger)
        y = LpVariable.dicts("y", graph.nodes(), 0, 1, LpInteger)
        e = LpVariable.dicts("e", (tuple(sorted(edge)) for edge in graph.edges()), 0, 1, LpInteger)

        # Objective function
        prob += lpSum([x[node] for node in graph.nodes()])

        # Constraints
        for i in graph.nodes():
            for j in graph.neighbors(i):
                # Constraint 1: xi=eij
                prob += x[i] - e[(i, j)] == 0

                # Constraint 2: Σi eij wij≥yj γ λj Σi wij
                wij = graph[i][j]['weight']
                lambda_j = 1.0 / wij  # λ values are the reciprocal of the weights
                prob += lpSum([e[tuple(sorted((i, j)))] * wij for i in graph.neighbors(j)]) >= y[
                    j] * gamma * lambda_j * wij

                # Constraint 3: Σj yj≥α |Y|
        prob += lpSum([y[j] for j in graph.nodes()]) >= alpha * len(graph.nodes())

        # Constraint 4: For each patient, the top β fraction of the expression altered genes with highest weights (λj) are always covered.
        for p in graph.nodes():
            sorted_neighbors = sorted(graph.neighbors(p), key=lambda x: 1.0 / graph[p][x]['weight'], reverse=True)
            top_beta = sorted_neighbors[:int(beta * len(sorted_neighbors))]
            for j in top_beta:
                prob += y[j] == 1

        prob.solve()

        # Return the nodes included in the cover
        cover_set = [node for node in x if x[node].varValue == 1]
        return cover_set