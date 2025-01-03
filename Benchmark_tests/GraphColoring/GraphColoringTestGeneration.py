import sys
import os
from qubovert import boolean_var
from DIMACS_converter import DIMACStoList

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import ToolchainTestScripts as Toolchain
from create_Q_matrix import create_Q_matrix

MinShannonDimension = 5

# File for Qoolchain statistics, hence file path starts from the Qoolchain directory
statfile = "./../Benchmark_tests/GraphColoring/GraphColoringStats.txt"
if os.path.exists(statfile):
    os.remove(statfile)
# File for Qoolchain statistics, hence file path starts from the Qoolchain directory
statfile_sh = "./../Benchmark_tests/GraphColoring/GraphColoringStats_Shannon.txt"
if os.path.exists(statfile_sh):
    os.remove(statfile_sh)

chrn_file = "ChromaticNumbers.txt"

for filename in os.listdir("Benchmarks"):
    file_path = "Benchmarks/" + filename

    # Open chromatic number file
    try:
        chrn_f = open(chrn_file, "r")
    except:
        print(f'File {chrn_file} not found')
        exit(1)

    # Check chromatic number for the current graph
    for line in chrn_f:
        splitted_line = line.split()
        if splitted_line[0] == filename:
            chr_n = int(splitted_line[1])
    chrn_f.close()

    n_nodes, n_edges, graph = DIMACStoList(file_path)

    x = {(i, c): boolean_var('x(%d, %d)' % (i, c)) for i in range(n_nodes) for c in range(chr_n)}

    model = 0
    # Each node must have exactly one color
    for i in range(n_nodes):
        model += (sum(x[i, c] for c in range(chr_n)) - 1) ** 2

    # Adjacent nodes must have different colors
    for edge in graph:
        for c in range(chr_n):
            model += x[edge[0]-1, c] * x[edge[1]-1, c]

    qubo = model.to_qubo()
    d, Q, off = create_Q_matrix(qubo, n_nodes * chr_n)

    try:
        f = open("Results/" + filename.split('.')[0] + "_n.txt", "w")
    except:
        print(f'Error in {filename} QUBO matrix file creation')
        exit(1)

    nz_terms = 0
    for i in range(n_nodes * chr_n):
        for j in range(n_nodes * chr_n):
            f.write(format(Q[i][j]) + " ")
            if j > i and Q[i][j] != 0:
                nz_terms += 1  # Count nonzero terms to calculate the density
        f.write("\n")
    f.write("\n")
    f.close()

    try:
        fSol = open("Results/" + filename.split('.')[0] + "_n_sol.txt", "w")
    except:
        print(f'Error in {filename} solution file creation')
        exit(1)

    print("Problem " + filename + " Simulated annealing")
    opt_value = Toolchain.simulated_annealing(fSol, model)

    print("Problem " + filename + " Toolchain")
    Toolchain.toolchain_simulation(fSol, statfile, Q, off, opt_value)

    print("Problem " + filename + " Shannon")
    Toolchain.toolchain_shannon_simulation(fSol, statfile_sh, Q, 0, MinShannonDimension, opt_value)

    # Calculate maximum value and upper bound for number of qubits estimation
    print("Maximum and upper bound calculation")
    model = -model
    qubo = model.to_qubo()
    d, Q, off = create_Q_matrix(qubo, n_nodes * chr_n)
    opt_value = Toolchain.simulated_annealing(fSol, model)
    Toolchain.toolchain_simulation(fSol, "", Q, off, opt_value)

    # Add the density of the QUBO matrix to the solution file
    fSol.write("Density: " + format(float(nz_terms) / ((((n_nodes * chr_n) ** 2) / 2) - (n_nodes * chr_n))) + "\n")

    fSol.close()
