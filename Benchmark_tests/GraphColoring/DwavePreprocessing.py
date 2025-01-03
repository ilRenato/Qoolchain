import sys

sys.path.insert(0, r'./../Dwave_Toolchain')
#import dwave.preprocessing
from dwave.preprocessing.composites import FixVariablesComposite
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.samplers import SimulatedAnnealingSampler
from create_Q_matrix import create_qv_model
import dimod
import os

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

for filename in os.listdir("data"):
    try:
        f = open("data/" + filename, "r")
    except:
        print(f'Error in {filename} file opening')
        exit(1)

    print(filename)
    resfile.write(filename + '\n')
    bound_file.write(filename + '\n')

    Q = []
    for line in f:
        splitted_line = line.split()
        if len(splitted_line) <= 1:
            continue
        Q.append([int(num) for num in splitted_line if num != '\n'])
    f.close()

    model = create_qv_model(Q, [i for i in range(len(Q))])
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
    model = create_qv_model(Q, [i for i in range(len(Q))])
    qubo = model.to_qubo()
    up_bound = -roof_duality(dimod.BinaryQuadraticModel.from_qubo(qubo.Q))[0]

    # Calculate range
    f_range = up_bound - low_bound
    bound_file.write(format(f_range) + '\n')

    resfile.write(format(len(fixed_persistency)) + ' ' + format(fixed_persistency) + '\n')

resfile.close()
bound_file.close()
