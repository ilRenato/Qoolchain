import math
import qubovert as qv
from create_Q_matrix import create_qv_model, eliminateNodes
from CompareResult import evaluateResult
from qubovert.sim import anneal_quso, anneal_qubo
from qubovert import boolean_var
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../Toolchain_implementation")))
import QuboPreprocessing


def simulated_annealing(resfile, model):
    qubo = model.to_qubo()
    anneal_res = qv.sim.anneal_qubo(qubo, num_anneals=30, anneal_duration=10000)
    opt_value = anneal_res.best.value
    resfile.write("Energy " + format(anneal_res.best.value) + "\n")
    resfile.write("Solution " + format(anneal_res.best.state) + "\n")
    return opt_value


def toolchain_simulation(resfile, statfile, Q, offset, opt_value):
    subQ, varnames, persistencies, lower_bound = QuboPreprocessing.PreQuboWrap(Q, offset, bytes(statfile, encoding='utf8'))
    total_solution = {}
    # QUBO subfunctions are scanned in the opposite way since the toolchain returns last the complete
    # strongly connected components that have to be solved first
    for i in range(len(varnames) - 1, -1, -1):
        # skip subqubo if it has no elements, thus every persistency has been fixed
        if len(varnames[i]) != 0:
            # check for already solved variables
            solved_vars = {}
            solved_vars_found = False
            for var in varnames[i]:
                if var in total_solution:
                    solved_vars[var] = total_solution[var]
                    solved_vars_found = True
            # If there are solved variables substitute them in the function
            if solved_vars_found:
                currentQ, current_varnames = eliminateNodes(subQ[i], varnames[i], solved_vars)
            else:
                currentQ = subQ[i]
                current_varnames = varnames[i]
            # set variables value to True in case of zero matrix
            mat0 = [[0 for _ in range(len(current_varnames))] for _ in range(len(current_varnames))]
            if currentQ == mat0:
                for var in current_varnames:
                    total_solution[var] = True
            else:
                model = create_qv_model(currentQ, current_varnames)
                # Check if there are nodes with no connections that can be assigned to any binary value
                for element in current_varnames:
                    # If an element has no connections it does not appear in the qubovert model
                    if ('x(' + format(element) + ')') not in model.variables:
                        total_solution[element] = True
                subqubo = model.to_qubo()
                anneal_res = qv.sim.anneal_qubo(subqubo, num_anneals=30, anneal_duration=10000)
                for key, val in model.convert_solution(anneal_res.best.state).items():
                    new_key = int(key.split('(')[1].split(')')[0])
                    if val == 1:
                        new_value = True
                    else:
                        new_value = False
                    total_solution[new_key] = new_value
                # print(model.convert_solution(anneal_res.best.state))
    total_solution.update(persistencies)
    cpp_value = evaluateResult(Q, total_solution, offset)
    resfile.write("CppEnergy " + format(cpp_value) + "\n")
    resfile.write("CppSolution " + format(total_solution) + "\n")
    resfile.write("Lower bound " + format(lower_bound) + "\n")
    if lower_bound > opt_value:
        resfile.write("Lower_bound_Error\n")
    else:
        resfile.write("Lower_bound_OK\n")
    if cpp_value == opt_value:
        resfile.write("OK\n")
    else:
        resfile.write("Error " + format(cpp_value - opt_value) + "\n")


def toolchain_shannon_simulation(resfile, statfile, Q, offset, min_dimension, opt_value):
    sh_iterations = 5
    largest_size_qubo = 0
    subQ, varnames, persistencies, lower_bound = QuboPreprocessing.PreQuboWrap(Q, offset, bytes(statfile, encoding='utf8'))
    total_solution = {}
    # QUBO subfunctions are scanned in the opposite way since the toolchain returns last the complete
    # strongly connected components that have to be solved first
    for i in range(len(varnames) - 1, -1, -1):
        # skip subqubo if it has no elements, thus every persistency has been fixed
        if len(varnames[i]) != 0:
            # check for already solved variables
            solved_vars = {}
            solved_vars_found = False
            for var in varnames[i]:
                if var in total_solution:
                    solved_vars[var] = total_solution[var]
                    solved_vars_found = True
            # If there are solved variables substitute them in the function
            if solved_vars_found:
                currentQ, current_varnames = eliminateNodes(subQ[i], varnames[i], solved_vars)
            else:
                currentQ = subQ[i]
                current_varnames = varnames[i]
            # set variables value to True in case of zero matrix
            mat0 = [[0 for _ in range(len(current_varnames))] for _ in range(len(current_varnames))]
            if currentQ == mat0:
                for var in current_varnames:
                    total_solution[var] = True
            else:
                if len(current_varnames) < min_dimension:
                    if len(current_varnames) > largest_size_qubo:
                        largest_size_qubo = len(current_varnames)
                    model = create_qv_model(currentQ, current_varnames)
                    # Check if there are nodes with no connections that can be assigned to any binary value
                    for element in current_varnames:
                        # If an element has no connections it does not appear in the qubovert model
                        if ('x(' + format(element) + ')') not in model.variables:
                            total_solution[element] = True
                    subqubo = model.to_qubo()
                    anneal_res = qv.sim.anneal_qubo(subqubo, num_anneals=30, anneal_duration=10000)
                    for key, val in model.convert_solution(anneal_res.best.state).items():
                        new_key = int(key.split('(')[1].split(')')[0])
                        if val == 1:
                            new_value = True
                        else:
                            new_value = False
                        total_solution[new_key] = new_value
                    # print(model.convert_solution(anneal_res.best.state))
                else:
                    shanQubos_list = QuboPreprocessing.ShannonDecompositionWrap(currentQ, sh_iterations, 0)
                    solutions = []
                    results = []
                    for item, shanPersist in shanQubos_list:
                        tmp_solution = {}
                        for subf, shanVarnames in item:
                            if len(shanVarnames) != 0:
                                # check for already solved variables
                                solved_vars = {}
                                solved_vars_found = False
                                for var in shanVarnames:
                                    if var in tmp_solution:
                                        solved_vars[var] = tmp_solution[var]
                                        solved_vars_found = True
                                # If there are solved variables substitute them in the function
                                if solved_vars_found:
                                    new_subf, new_shanVarnames = eliminateNodes(subf, shanVarnames, solved_vars)
                                else:
                                    new_subf = subf
                                    new_shanVarnames = shanVarnames
                                # set variables value to True in case of zero matrix
                                mat0 = [[0 for _ in range(len(new_shanVarnames))] for _ in range(len(new_shanVarnames))]
                                if new_subf == mat0:
                                    for var in new_shanVarnames:
                                        tmp_solution[var] = True
                                else:
                                    model = create_qv_model(new_subf, new_shanVarnames)
                                    # Check if there are nodes with no connections that can be assigned to any binary value
                                    for element in new_shanVarnames:
                                        # If an element has no connections it does not appear in the qubovert model
                                        if ('x(' + format(element) + ')') not in model.variables:
                                            tmp_solution[element] = True
                                    if len(model.variables) > largest_size_qubo:
                                        largest_size_qubo = len(model.variables)
                                    subqubo = model.to_qubo()
                                    anneal_res = qv.sim.anneal_qubo(subqubo, num_anneals=30, anneal_duration=10000)
                                    for key, val in (model.convert_solution(anneal_res.best.state)).items():
                                        new_key = int(key.split('(')[1].split(')')[0])
                                        if val == 1:
                                            new_value = True
                                        else:
                                            new_value = False
                                        tmp_solution[new_key] = new_value
                        tmp_solution.update(shanPersist)
                        cpp_value = evaluateResult(currentQ, tmp_solution, 0)
                        solutions.append(tmp_solution)
                        results.append(cpp_value)
                    best_res = min(results)
                    best_sol_index = results.index(best_res)
                    best_sol = solutions[best_sol_index]
                    sol_with_names = {}
                    for k, v in best_sol.items():
                        sol_with_names[current_varnames[k]] = v
                    total_solution.update(sol_with_names)
                    # print(sol_with_names)

    total_solution.update(persistencies)
    cpp_value = evaluateResult(Q, total_solution, offset)
    resfile.write("ShannonEnergy " + format(cpp_value) + "\n")
    resfile.write("ShannonSolution " + format(total_solution) + "\n")
    resfile.write("Largest_Shannon_QUBO " + format(largest_size_qubo) + "\n")
    if cpp_value <= opt_value:
        resfile.write("Shannon_OK\n")
    else:
        resfile.write("Shannon_Error " + format(cpp_value - opt_value) + "\n")
