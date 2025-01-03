import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

import os
import sys
sys.path.insert(0, r'./..')
from bounds_utils import *
from create_Q_matrix import create_qv_model
from qubovert import boolean_var

M = [3, 5, 10, 15]
MaxProblemDimension = 40
MinProblemDimension = 10
ProblemForEachDimension = 10
N_problems = int(((MaxProblemDimension - MinProblemDimension) / 2) + 1)

try:
    errfile = open("error.log", "w")
except:
    print("Error in file creation\n")
    exit(1)

try:
    shannon_errfile = open("sh_error.log", "w")
except:
    print("Error in Shannon file creation\n")
    exit(1)

# Dwave data
try:
    dwave_file = open("Dwave_bounds.txt", "r")
except:
    print("Error in Dwave file opening\n")
    exit(1)

dwave_list = []
i = -1
for line in dwave_file:
    splitted_line = line.split()
    if splitted_line[0] == "Density:":
        i += 1
        dwave_list.append({})
    elif splitted_line[0].isdigit():
        n_nodes = int(splitted_line[0])
        dwave_list[i][n_nodes] = []
    else:
        dwave_list[i][n_nodes].append(float(splitted_line[0]))
dwave_file.close()

err_dict = {}
err_bet_dict = {}
err_wor_dict = {}
value_err_list = []
i_dwave = 0
for maxnum in M:
    err_dict[maxnum] = []
    err_bet_dict[maxnum] = []
    err_wor_dict[maxnum] = []
    for dim in range(MinProblemDimension, MaxProblemDimension + 1, 2):
        lower_bound_list = []
        posneg_list = []
        naive_list = []
        opt_value_list = []
        for num_prob in range(ProblemForEachDimension):

            try:
                resfile = open("Results/NumberPartitioning_" + format(dim) + "_nodes_M" + format(maxnum) + "_" + format(num_prob) + "_n_sol.txt", "r")
            except:
                print("Error in file opening\n")
                exit(1)

            flag = False
            for line in resfile:
                splitted_line = line.split(' ')
                if splitted_line[0] == "Error":
                    errfile.write("ERROR: NumberPartitioning_" + format(dim) + "_nodes_M" + format(maxnum) + "_" + format(num_prob) + "_n_sol.txt ")
                    errfile.write(format(int(float(splitted_line[1]))) + "\n")
                    value_err_list.append(int(splitted_line[1]))
                if splitted_line[0] == "OK\n":
                    value_err_list.append(0)
                if splitted_line[0] == "Shannon_Error":
                    shannon_errfile.write("ERROR: NumberPartitioning_" + format(dim) + "_nodes_M" + format(maxnum) + "_" + format(num_prob) + "_n_sol.txt ")
                    shannon_errfile.write(format(int(float(splitted_line[1]))) + "\n")
                if splitted_line[0] == "Energy":
                    if not flag:
                        true_range = -float(splitted_line[1])
                        flag = True
                    else:
                        true_range -= float(splitted_line[1])
                        flag = False
                        opt_value_list.append(true_range)
                if splitted_line[0] == "Lower":
                    if flag:
                        est_range = -float(splitted_line[2])
                    else:
                        est_range -= float(splitted_line[2])
                        lower_bound_list.append(est_range)
            resfile.close()

            # Bound evaluation with different methods
            try:
                resfile = open("Results/NumberPartitioning_" + format(dim) + "_nodes_M" + format(maxnum) + "_" + format(num_prob) + "_n.txt", "r")
            except:
                print("Error in file opening\n")
                exit(1)

            Q = []
            for i in range(dim):
                Q.append([int(x) for x in resfile.readline().split()])
            resfile.close()
            model = create_qv_model(Q, list(range(dim)))
            qubo = model.to_qubo()
            lower_bet, upper_bet = bounds_pos_neg(qubo)
            range_bet = upper_bet - lower_bet
            posneg_list.append(range_bet)
            lower_wor, upper_wor = bounds_naive(qubo)
            range_wor = upper_wor - lower_wor
            naive_list.append(range_wor)

        err_avg = np.average([(lower_bound_list[i] - opt_value_list[i]) * 100 / opt_value_list[i] for i in range(len(opt_value_list))])
        err_dict[maxnum].append(err_avg)
        posneg_avg = np.average([(posneg_list[i] - opt_value_list[i]) * 100 / opt_value_list[i] for i in range(len(opt_value_list))])
        err_bet_dict[maxnum].append(posneg_avg)
        naive_avg = np.average([(naive_list[i] - opt_value_list[i]) * 100 / opt_value_list[i] for i in range(len(opt_value_list))])
        err_wor_dict[maxnum].append(naive_avg)
        dwave_list[i_dwave][dim] = np.average([(dwave_list[i_dwave][dim][i] - opt_value_list[i]) * 100 / opt_value_list[i] for i in range(len(opt_value_list))])
    i_dwave += 1
value_err = np.average(value_err_list)
max_err = max(value_err_list)
print("Average error: " + format(value_err))
print("Maximum error: " + format(max_err))
errfile.close()
shannon_errfile.close()

cmap = plt.get_cmap('RdYlGn')
colors = cmap(np.linspace(1, 0, 4))
for i in range(len(M)):
    plt.plot(np.linspace(MinProblemDimension, MaxProblemDimension, N_problems), err_dict[M[i]], label=format(M[i]) + "%", linewidth=2, color=colors[i])
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Error over estimated range [%]', fontsize=15)
plt.legend(title="Density")
plt.savefig("Lower_bound.eps", format='eps', bbox_inches='tight')
plt.savefig("Lower_bound.png", format='png', bbox_inches='tight')
plt.savefig("Lower_bound.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Other methods plot
patches = []
for i in range(len(M)):
    plt.plot(np.linspace(MinProblemDimension, MaxProblemDimension, N_problems), err_dict[M[i]], linewidth=2, linestyle='dashed', color=colors[i])
    plt.plot(np.linspace(MinProblemDimension, MaxProblemDimension, N_problems), err_bet_dict[M[i]], linewidth=2, linestyle='dotted', color=colors[i])
    patches.append(mpatches.Patch(color=colors[i], label=format(M[i]) + "%"))
patches.append(mlines.Line2D([], [], color='0', linestyle='dashed', label="posneg"))
patches.append(mlines.Line2D([], [], color='0', linestyle='dotted', label="naive"))
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Error over estimated range [%]', fontsize=15)
plt.legend(title="Density", handles=patches)
plt.savefig("Lower_bound_compare.eps", format='eps', bbox_inches='tight')
plt.savefig("Lower_bound_compare.png", format='png', bbox_inches='tight')
plt.savefig("Lower_bound_compare.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Dwave comparison
patches = []
for i in range(len(M)):
    plt.plot(np.linspace(MinProblemDimension, MaxProblemDimension, N_problems), err_dict[M[i]], linewidth=2, color=colors[i])
    plt.plot(np.linspace(MinProblemDimension, MaxProblemDimension, N_problems), dwave_list[i].values(), linewidth=2, linestyle='dashed', color=colors[i])
    patches.append(mpatches.Patch(color=colors[i], label=format(M[i]) + "%"))
patches.append(mlines.Line2D([], [], color='0', linestyle='solid', label="Qoolchain"))
patches.append(mlines.Line2D([], [], color='0', linestyle='dashed', label="Dwave"))
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Error over estimated range [%]', fontsize=15)
plt.legend(title="Density", handles=patches)
plt.savefig("Range_dwave.eps", format='eps', bbox_inches='tight')
plt.savefig("Range_dwave.png", format='png', bbox_inches='tight')
plt.savefig("Range_dwave.pdf", format='pdf', bbox_inches='tight')
plt.close()
