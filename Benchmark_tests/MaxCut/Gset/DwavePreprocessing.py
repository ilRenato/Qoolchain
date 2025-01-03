import sys

sys.path.insert(0, r'./../../Dwave_Toolchain')
import dwave.preprocessing
from dwave.preprocessing.composites import FixVariablesComposite
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.samplers import SimulatedAnnealingSampler
from create_Q_matrix import create_qv_model
import dimod
import os
from Converter import create_QMatrix

MAX_gset = 21

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

for i in range(1, MAX_gset + 1):
    filename = "G" + format(i) + ".txt"
    try:
        f = open("GSets_u800/" + filename, "r")
    except:
        print(f'Error in {filename} file opening')
        exit(1)

    print(filename)
    resfile.write(filename + '\n')
    bound_file.write(filename + '\n')

    Q, nodes = create_QMatrix("GSets_u800/" + filename)

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
