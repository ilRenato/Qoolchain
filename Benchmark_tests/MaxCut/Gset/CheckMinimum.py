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
from Converter import create_QMatrix

MAX_gset = 21

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

dwave_dict = {}
for line in dwave_file:
    splitted_line = line.split()
    if not splitted_line[0].replace('.', '').isdigit():
        benchmark = splitted_line[0][:-4]
    else:
        dwave_dict[benchmark] = float(splitted_line[0])
dwave_file.close()

err_dict = {}
err_bet_dict = {}
err_wor_dict = {}
opt_value_dict = {}
value_err_list = []
for i in range(1, MAX_gset + 1):
    filename = "G" + format(i) + ".txt"
    try:
        resfile = open("Results/" + filename, "r")
    except:
        print(f'Error in {filename} file opening\n')
        exit(1)

    flag = False
    for line in resfile:
        splitted_line = line.split(' ')
        if splitted_line[0] == "Error":
            errfile.write("ERROR: " + filename)
            errfile.write(format(int(float(splitted_line[1]))) + "\n")
            value_err_list.append(int(splitted_line[1]))
        if splitted_line[0] == "OK\n":
            value_err_list.append(0)
        if splitted_line[0] == "Shannon_Error":
            shannon_errfile.write("ERROR: " + filename)
            shannon_errfile.write(format(int(float(splitted_line[1]))) + "\n")
        if splitted_line[0] == "Energy":
            if not flag:
                true_range = -float(splitted_line[1])
                flag = True
            else:
                true_range -= float(splitted_line[1])
                flag = False
                opt_value = true_range
                opt_value_dict[filename[:-4]] = opt_value
        if splitted_line[0] == "Lower":
            if flag:
                est_range = -float(splitted_line[2])
            else:
                est_range -= float(splitted_line[2])
                lower_bound = est_range
    resfile.close()
    err_dict[filename[:-4]] = (lower_bound - opt_value) * 100 / opt_value
    dwave_dict[filename[:-4]] = (dwave_dict[filename[:-4]] - opt_value) * 100 / opt_value

    # Bound evaluation with different methods
    Q, _ = create_QMatrix("GSets/" + filename)
    model = create_qv_model(Q, list(range(len(Q))))
    qubo = model.to_qubo()
    lower_bet, upper_bet = bounds_pos_neg(qubo)
    range_bet = upper_bet - lower_bet
    err_bet_dict[filename[:-4]] = (range_bet - opt_value_dict[filename[:-4]]) * 100 / opt_value_dict[filename[:-4]]
    lower_wor, upper_wor = bounds_naive(qubo)
    range_wor = upper_wor - lower_wor
    err_wor_dict[filename[:-4]] = (range_wor - opt_value_dict[filename[:-4]]) * 100 / opt_value_dict[filename[:-4]]

value_err = np.average(value_err_list)
max_err = max(value_err_list)
print("Average error: " + format(value_err))
print("Maximum error: " + format(max_err))
errfile.close()
shannon_errfile.close()

x = np.arange(len(err_dict))
plt.bar(x, err_dict.values())
plt.xlabel('Benchmark', fontsize=15)
plt.xticks(x, err_dict.keys(), rotation=30)
plt.ylabel('Error over estimated range [%]', fontsize=15)
plt.savefig("Lower_bound.eps", format='eps', bbox_inches='tight')
plt.savefig("Lower_bound.png", format='png', bbox_inches='tight')
plt.savefig("Lower_bound.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Other methods plot
cmap = plt.get_cmap('plasma')
colors = cmap(np.linspace(1, 0, 3))
width = 0.25
plt.bar(x, err_dict.values(), label='Qoolchain', color=colors[0], width=width)
plt.bar(x + width, err_bet_dict.values(), label='posneg', color=colors[1], width=width)
plt.bar(x + 2*width, err_wor_dict.values(), label='naive', color=colors[2], width=width)
plt.xlabel('Benchmark', fontsize=15)
plt.xticks(x + width, err_dict.keys(), rotation=30)
plt.ylabel('Error over estimated range [%]', fontsize=15)
plt.legend(loc='upper left')
plt.savefig("Lower_bound_compare.eps", format='eps', bbox_inches='tight')
plt.savefig("Lower_bound_compare.png", format='png', bbox_inches='tight')
plt.savefig("Lower_bound_compare.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Dwave comparison
cmap = plt.get_cmap('plasma')
colors = cmap(np.linspace(1, 0, 2))
plt.bar(x, err_dict.values(), label='Qoolchain', color=colors[0], width=width)
plt.bar(x + width, dwave_dict.values(), label='dwave', color=colors[1], width=width)
plt.xlabel('Benchmark', fontsize=15)
plt.xticks(x + width/2, err_dict.keys(), rotation=30)
plt.ylabel('Error over estimated range [%]', fontsize=15)
plt.legend(loc='upper left')
plt.savefig("Range_dwave.eps", format='eps', bbox_inches='tight')
plt.savefig("Range_dwave.png", format='png', bbox_inches='tight')
plt.savefig("Range_dwave.pdf", format='pdf', bbox_inches='tight')
plt.close()
