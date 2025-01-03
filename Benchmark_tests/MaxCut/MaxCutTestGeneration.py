import numpy as np
import random
import math
from qubovert import boolean_var
import qubovert as qv

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from create_Q_matrix import create_Q_matrix
import ToolchainTestScripts as Toolchain

MaxProblemDimention = 60
MinProblemDimention = 10
ProblemForEachDimention = 10
MinShannonDimension = 5

# Probability percentage of having an edge between each couple of nodes
Density = [2, 4, 8, 16]

for Pone in Density:
    # File for Qoolchain statistics, hence file path starts from the Qoolchain directory
    statfile = "./../Benchmark_tests/MaxCut/MaxCutStats_P{}.txt".format(Pone)
    statfile_sh = "./../Benchmark_tests/MaxCut/MaxCutStats_P{}_Shannon.txt".format(Pone)
    if os.path.exists(statfile):
        os.remove(statfile)

    for numberOfNodes in range(MinProblemDimention, MaxProblemDimention + 1, 2):
        print(numberOfNodes)
        for Problem in range(ProblemForEachDimention):
            print("Problem N" + format(Problem))
            fileName = "Results/MaxCut_" + format(numberOfNodes) + "_nodes_" + "P" + format(Pone) + "_" +\
                                   format(Problem) + "_n.txt"

            try:
                f = open(fileName, "w")
            except:
                print("Error in file creation\n")
                exit(1)

            graph_ok = False
            while not graph_ok:

                Weigth = np.zeros((numberOfNodes, numberOfNodes))
                for i in range(numberOfNodes):
                    for j in range(i+1, numberOfNodes):
                        prob = random.choices([0, 1], weights=(100-Pone, Pone), k=1)[0]
                        if prob == 1:
                            Weigth.itemset((i, j), random.choice(range(1, 7)))
                        else:
                            Weigth.itemset((i, j), 0)

                x = {i: boolean_var('x(%d)' % i) for i in range(numberOfNodes)}

                model = 0

                for i in range(numberOfNodes):
                    for j in range(i+1, numberOfNodes):
                        model -= Weigth.item((i,j))*(x[j]+x[i]-2*x[i]*x[j])

                # Add the variables that are not already present in the model
                for i in range(numberOfNodes):
                    if ('x(' + format(i) + ')') not in model.variables:
                        model[('x(' + format(i) + ')'), ] = 0

                qubo = model.to_qubo()

                if len(qubo.variables) == numberOfNodes:
                    graph_ok = True

            d, Q, off = create_Q_matrix(qubo, numberOfNodes)

            for i in range(numberOfNodes):
                for j in range(numberOfNodes):
                    f.write(format(Q[i][j]) + " ")
                f.write("\n")
            f.write("\n")
            f.close()

            fileNameSol = "Results/MaxCut_" + format(numberOfNodes) + "_nodes_" + "P" + format(Pone) + "_" + \
                          format(Problem) + "_n_sol.txt"

            try:
                fSol = open(fileNameSol, "w")
            except:
                print("Error in file creation\n")
                exit(1)

            print("Problem " + format(numberOfNodes) + " " + format(Pone) + " " + format(Problem) + " Simulated annealing")
            opt_value = Toolchain.simulated_annealing(fSol, model)

            print("Problem " + format(numberOfNodes) + " " + format(Pone) + " " + format(Problem) + " Toolchain")
            Toolchain.toolchain_simulation(fSol, statfile, Q, off, opt_value)

            print("Problem " + format(numberOfNodes) + " " + format(Pone) + " " + format(Problem) + " Shannon")
            Toolchain.toolchain_shannon_simulation(fSol, statfile_sh, Q, 0, 5, opt_value)

            # Calculate maximum value and upper bound for number of qubits estimation
            print("Maximum and upper bound calculation")
            model = -model
            qubo = model.to_qubo()
            d, Q, off = create_Q_matrix(qubo, numberOfNodes)
            opt_value = Toolchain.simulated_annealing(fSol, model)
            Toolchain.toolchain_simulation(fSol, "", Q, off, opt_value)

            fSol.close()