import sys

sys.path.insert(0, r'./../Dwave_Toolchain')
import dwave.preprocessing
from dwave.preprocessing.composites import FixVariablesComposite
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.samplers import SimulatedAnnealingSampler
from create_Q_matrix import create_qv_model
import dimod

MaxProblemDimension = 100
MinProblemDimension = 10
ProblemForEachDimension = 10

# Probability percentage of having an edge between each couple of nodes
Density = [80, 90, 95, 98]

sampler = SimulatedAnnealingSampler()

try:
    resfile = open("Dwave_toolchain_res.txt", "w")
except:
    print("Error in file creation\n")
    exit(1)

try:
    bound_file = open("Dwave_bounds.txt", "w")
except:
    print("Error in file creation\n")
    exit(1)

for Pone in Density:
    resfile.write("Density: " + format(Pone) + "\n")
    bound_file.write("Density: " + format(Pone) + "\n")
    for numberOfNodes in range(MinProblemDimension, MaxProblemDimension + 1, 2):
        resfile.write(format(numberOfNodes) + '\n')
        bound_file.write(format(numberOfNodes) + '\n')
        for Problem in range(ProblemForEachDimension):
            print("P" + format(Pone) + ", Nodes: " + format(numberOfNodes) + ", N" + format(Problem) + "\n")
            fileName = "Results/MaxClique_" + format(numberOfNodes) + "_nodes_" + "P" + format(Pone) + "_" + \
                       format(Problem) + "_n.txt"

            try:
                f = open(fileName, "r")
            except:
                print(fileName)
                print("Error in file opening\n")
                exit(1)

            Q = []
            for line in f:
                splitted_line = line.split(' ')
                if len(splitted_line) == 1:
                    continue
                Q.append([int(num) for num in splitted_line if num != '\n'])

            model = create_qv_model(Q, [i for i in range(numberOfNodes)])
            qubo = model.to_qubo()

            # Fix variables
            sampler_persistency = FixVariablesComposite(sampler, algorithm='roof_duality')
            _, fixed_persistency = sampler_persistency.sample_qubo(qubo.Q, strict=False)

            # Lower bound
            low_bound = roof_duality(dimod.BinaryQuadraticModel.from_qubo(qubo.Q))[0]

            # Upper bound
            for i in range(len(Q)):
                for j in range(len(Q[i])):
                    Q[i][j] = -Q[i][j]
            model = create_qv_model(Q, [i for i in range(numberOfNodes)])
            qubo = model.to_qubo()
            up_bound = -roof_duality(dimod.BinaryQuadraticModel.from_qubo(qubo.Q))[0]

            # Calculate range
            f_range = up_bound - low_bound
            bound_file.write(format(f_range) + '\n')

            resfile.write(format(len(fixed_persistency)) + ' ' + format(fixed_persistency) + '\n')

resfile.close()
bound_file.close()
