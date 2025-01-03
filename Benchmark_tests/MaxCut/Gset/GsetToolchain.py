from Converter import create_QMatrix
import numpy as np
from qubovert import boolean_var

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from create_Q_matrix import create_qv_model, create_Q_matrix
import ToolchainTestScripts as Toolchain

MinShannonDimension = 5

# File for Qoolchain statistics, hence file path starts from the Qoolchain directory
statfile = "./../Benchmark_tests/MaxCut/Gset/GsetStats.txt"
if os.path.exists(statfile):
    os.remove(statfile)
# File for Qoolchain statistics, hence file path starts from the Qoolchain directory
statfile_sh = "./../Benchmark_tests/MaxCut/Gset/GsetStats_shannon.txt"
if os.path.exists(statfile_sh):
    os.remove(statfile_sh)

result_file = "BestKnownValue.txt"

i = 0
for gfile in os.listdir("GSets_u800"):
    gfile_path = "GSets_u800/" + gfile

    QMatrix, nodes = create_QMatrix(gfile_path)
    model = create_qv_model(QMatrix, list(range(nodes)))

    fileNameSol = "Results/" + gfile + ".txt"

    try:
        fSol = open(fileNameSol, "w")
    except:
        print("Error in file creation\n")
        exit(1)

    try:
        rfile = open(result_file, "r")
    except:
        print("Error in file creation\n")
        exit(1)

    for line in rfile:
        splitted_line = line.split(':')
        if splitted_line[0] == gfile.split('.')[0]:
            best_known_value = int(splitted_line[1])

    print("Problem {} ".format(i) + gfile + " Simulated annealing")
    Toolchain.simulated_annealing(fSol, model)

    print("Problem " + gfile + " Toolchain")
    Toolchain.toolchain_simulation(fSol, statfile, QMatrix, 0, best_known_value)

    print("Problem " + gfile + " Shannon")
    Toolchain.toolchain_shannon_simulation(fSol, statfile_sh, QMatrix, 0, MinShannonDimension, best_known_value)

    # Calculate maximum value and upper bound for number of qubits estimation
    print("Maximum and upper bound calculation")
    model = -model
    qubo = model.to_qubo()
    d, Q, off = create_Q_matrix(qubo, len(QMatrix))
    opt_value = Toolchain.simulated_annealing(fSol, model)
    Toolchain.toolchain_simulation(fSol, "", Q, off, opt_value)

    fSol.close()
    rfile.close()
    i += 1
