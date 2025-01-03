import qubovert as qv
from qubovert.sim import anneal_quso, anneal_qubo
from qubovert import boolean_var


def create_Q_matrix(qubo, num_variables):
    Q = [[0 for _ in range(num_variables)] for _ in range(num_variables)]
    offset = 0
    NVar = len(qubo.variables)
    for key in qubo.keys():
        if len(key) == 1:
            Q[key[0]][key[0]] = int(qubo[key])
        elif len(key) == 2:
            Q[key[0]][key[1]] = int(qubo[key])
        elif len(key) == 0:
            offset = qubo[key]

    return NVar, Q, offset


def create_qv_model(Q, varnames):
    x = {i: boolean_var('x(%d)' % varnames[i]) for i in range(len(varnames))}
    model = 0
    for i in range(len(Q)):
        for j in range(i+1, len(Q)):
            model += Q[i][j]*x[i]*x[j]
    for i in range(len(Q)):
        model += Q[i][i]*x[i]
    return model


def eliminateNodes(Q, nodes_list, nodes_assignments):
    newQ = [[0 for _ in range(len(nodes_list) - len(nodes_assignments))] for _ in range(len(nodes_list) - len(nodes_assignments))]
    new_nodes_list = [i for i in nodes_list if i not in nodes_assignments.keys()]
    new_i = 0
    for i in range(len(nodes_list)):
        new_j = new_i
        for j in range(i, len(nodes_list)):
            if (nodes_list[i] not in nodes_assignments) and (nodes_list[j] not in nodes_assignments):
                newQ[new_i][new_j] += Q[i][j]
            if (nodes_list[i] in nodes_assignments) and (nodes_list[j] not in nodes_assignments):
                if nodes_assignments[nodes_list[i]]:
                    newQ[new_j][new_j] += Q[i][j]
            if (nodes_list[i] not in nodes_assignments) and (nodes_list[j] in nodes_assignments):
                if nodes_assignments[nodes_list[j]]:
                    newQ[new_i][new_i] += Q[i][j]
            if nodes_list[j] not in nodes_assignments:
                new_j += 1
        if nodes_list[i] not in nodes_assignments:
            new_i += 1
    return newQ, new_nodes_list
