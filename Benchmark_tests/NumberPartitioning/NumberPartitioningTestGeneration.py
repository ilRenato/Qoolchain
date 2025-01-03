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

import QuboPreprocessing

MaxProblemDimention = 40
MinProblemDimention = 10
ProblemForEachDimention = 10
MinShannonDimension = 5

# Maximum number present in the set
MaxNumList = [3, 5, 10, 15]

for MaxNum in MaxNumList:
    # File for Qoolchain statistics, hence file path starts from the Qoolchain directory
    statfile = "./../Benchmark_tests/NumberPartitioning/NumPartStats_P{}.txt".format(MaxNum)
    statfile_sh = "./../Benchmark_tests/NumberPartitioning/NumPartStats_P{}_Shannon.txt".format(MaxNum)
    if os.path.exists(statfile):
        os.remove(statfile)

    for numberOfNumbers in range(MinProblemDimention, MaxProblemDimention + 1, 2):
        print(numberOfNumbers)
        for Problem in range(ProblemForEachDimention):
            fileName = "Results/NumberPartitioning_" + format(numberOfNumbers) + "_nodes_" + "M" + format(MaxNum) + "_" + \
                       format(Problem) + "_n.txt"

            try:
                f = open(fileName, "w")
            except:
                print("Error in file creation\n")
                exit(1)

            S = random.choices(range(1, MaxNum), k=numberOfNumbers)

            c = sum(S)

            x = {i: boolean_var('x(%d)' % i) for i in range(numberOfNumbers)}

            model = 0
            temp = 0
            for i in range(numberOfNumbers):
                temp += 2 * S[i] * x[i]

            model = (c - temp) ** 2

            qubo = model.to_qubo()

            d, Q, off = create_Q_matrix(qubo, numberOfNumbers)

            for i in range(numberOfNumbers):
                for j in range(numberOfNumbers):
                    f.write(format(Q[i][j]) + " ")
                f.write("\n")
            f.write("\n")
            f.close()

            fileNameSol = "Results/NumberPartitioning_" + format(numberOfNumbers) + "_nodes_" + "M" + format(MaxNum) +\
                          "_" + format(Problem) + "_n_sol.txt"

            try:
                fSol = open(fileNameSol, "w")
            except:
                print("Error in file creation\n")
                exit(1)

            print("Problem " + format(numberOfNumbers) + " " + format(MaxNum) + " " + format(Problem) + " Simulated annealing")
            opt_value = Toolchain.simulated_annealing(fSol, model)

            print("Problem " + format(numberOfNumbers) + " " + format(MaxNum) + " " + format(Problem) + " Toolchain")
            Toolchain.toolchain_simulation(fSol, statfile, Q, off, opt_value)

            print("Problem " + format(numberOfNumbers) + " " + format(MaxNum) + " " + format(Problem) + " Shannon")
            Toolchain.toolchain_shannon_simulation(fSol, statfile_sh, Q, 0, 6, opt_value)

            # Calculate maximum value and upper bound for number of qubits estimation
            print("Maximum and upper bound calculation")
            model = -model
            qubo = model.to_qubo()
            d, Q, off = create_Q_matrix(qubo, numberOfNumbers)
            opt_value = Toolchain.simulated_annealing(fSol, model)
            Toolchain.toolchain_simulation(fSol, "", Q, off, opt_value)

            fSol.close()
