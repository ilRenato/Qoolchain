import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import rc
import os

nodes_bench = {'1-FullIns_3': 120, '1-FullIns_4': 465, '2-FullIns_3': 260, '3-FullIns_3': 480, '4-FullIns_3': 798,
               '5-FullIns_3': 1232, 'myciel3': 44, 'myciel4': 115, 'myciel5': 282}

try:
    resfile = open("GraphColoringStats.txt", "r")
except:
    print("Error in file opening\n")
    exit(1)

persistency_dict = {}
all_pers_dict = {}
trivial_dec_dict = {}
var_reduction_dict = {}
scc_dec_dict = {}
probing_dict = {}
largest_q_dict = {}
tot_var_dict = {}
num_of_vars = 0
for line in resfile:
    if line[0] != '-':
        splitted_line = line.split()
        param = splitted_line[0]
        value = float(splitted_line[1].split('%')[0])
        if param == "Num_of_variables":
            num_of_vars = value
            if value not in persistency_dict:
                persistency_dict[value] = []
            if value not in all_pers_dict:
                all_pers_dict[value] = []
            if value not in trivial_dec_dict:
                trivial_dec_dict[value] = []
            if value not in scc_dec_dict:
                scc_dec_dict[value] = []
            if value not in probing_dict:
                probing_dict[value] = []
            if value not in var_reduction_dict:
                var_reduction_dict[value] = []
            if value not in largest_q_dict:
                largest_q_dict[value] = []
            if value not in tot_var_dict:
                tot_var_dict[value] = []
        if param == "Persistencies_percentage":
            persistency_dict[num_of_vars] = value
        if param == "Num_trivial_Decomposition":
            trivial_dec_dict[num_of_vars] = value
        if param == "Num_SCC_Decomposition":
            scc_dec_dict[num_of_vars] = value
        if param == "All_persistencies":
            all_pers_dict[num_of_vars] = (value * 100 / num_of_vars)
        if param == "Probing_persistencies":
            probing_dict[num_of_vars] = (value * 100 / num_of_vars)
        if param == "Variable_reduction_percentage":
            var_reduction_dict[num_of_vars] = value
        if param == "Largest_subfunction":
            largest_q_dict[num_of_vars] = value
for key in largest_q_dict:
    tot_var_dict[key] = ((key - largest_q_dict[key]) * 100 / key)

resfile.close()

# Collect results on Shannon decompositions
shannon_dict = {}
shannon_red_dict = {}
for filename in os.listdir("Results"):
    if not filename[-6:] == "_n.txt":
        try:
            sh_file = open("Results/" + filename, "r")
        except:
            print("Error in Shannon file opening\n")
            exit(1)

        dim = nodes_bench[filename[:-10]]

        for line in sh_file:
            splitted_line = line.split(' ')
            if splitted_line[0] == "Largest_Shannon_QUBO":
                shannon_dict[filename[:-10]] = (dim - int(splitted_line[1])) * 100 / dim
                shannon_red_dict[filename[:-10]] = (largest_q_dict[dim] - int(splitted_line[1])) * 100 / dim
        sh_file.close()


# Persistency plot
plt.bar(persistency_dict.keys(), persistency_dict.values())
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Persistencies [%]', fontsize=15)
plt.title('Persistencies fixed', fontsize=20)
#plt.legend()
plt.savefig("Persistencies.eps", format='eps', bbox_inches='tight')
plt.savefig("Persistencies.png", format='png', bbox_inches='tight')
plt.savefig("Persistencies.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Trivial decomposition plot
plt.bar(trivial_dec_dict.keys(), trivial_dec_dict.values())
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Number of decompositions', fontsize=15)
plt.title('Trivial decompositions', fontsize=20)
#plt.legend()
plt.savefig("Trivial_decomposition.eps", format='eps', bbox_inches='tight')
plt.savefig("Trivial_decomposition.png", format='png', bbox_inches='tight')
plt.savefig("Trivial_decomposition.pdf", format='pdf', bbox_inches='tight')
plt.close()

# SCC decomposition plot
plt.bar(scc_dec_dict.keys(), scc_dec_dict.values())
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Number of decompositions', fontsize=15)
plt.title('SCC decompositions', fontsize=20)
#plt.legend()
plt.savefig("SCC_decompositions.eps", format='eps', bbox_inches='tight')
plt.savefig("SCC_decompositions.png", format='png', bbox_inches='tight')
plt.savefig("SCC_decompositions.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Total variables reduction plot
plt.bar(var_reduction_dict.keys(), var_reduction_dict.values())
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Variable reduction [%]', fontsize=15)
plt.title('Total variable reduction', fontsize=20)
#plt.legend()
plt.savefig("Var_reduction.eps", format='eps', bbox_inches='tight')
plt.savefig("Var_reduction.png", format='png', bbox_inches='tight')
plt.savefig("Var_reduction.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Shannon reduction
plt.bar(shannon_dict.keys(), shannon_dict.values())
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Variable reduction [%]', fontsize=15)
#plt.legend()
plt.savefig("dwave-shannon.eps", format='eps', bbox_inches='tight')
plt.savefig("dwave-shannon.png", format='png', bbox_inches='tight')
plt.savefig("dwave-shannon.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Plot of all techniques together
pers = {}
prob = {}
tot = {}
for key in shannon_dict.keys():
    pers[key] = all_pers_dict[nodes_bench[key]]
    prob[key] = probing_dict[nodes_bench[key]]
    tot[key] = tot_var_dict[nodes_bench[key]]
cmap = plt.get_cmap('plasma')
colors = cmap(np.linspace(1, 0, 4))
x = list(range(len(shannon_dict)))
bottom = [0] * len(persistency_dict)
plt.bar(x, pers.values(), color=colors[0], label='Persistencies', bottom=bottom)
bottom = [bottom[i] + list(pers.values())[i] for i in range(len(bottom))]
plt.bar(x, prob.values(), color=colors[1], label='Probing', bottom=bottom)
bottom = [bottom[i] + list(pers.values())[i] + list(prob.values())[i] for i in range(len(bottom))]
plt.bar(x, tot.values(), color=colors[2], label='Decompositions', bottom=bottom)
bottom = [bottom[i] + list(pers.values())[i] + list(prob.values())[i] + list(tot.values())[i] for i in range(len(bottom))]
plt.bar(x, shannon_dict.values(), color=colors[3], label='Shannon', bottom=bottom)
plt.xlabel('Number of nodes', fontsize=15)
plt.xticks(x, labels=shannon_dict.keys(), rotation=30)
plt.ylabel('Variable reduction [%]', fontsize=15)
plt.legend()
plt.savefig("Reduction_comparison.eps", format='eps', bbox_inches='tight')
plt.savefig("Reduction_comparison.png", format='png', bbox_inches='tight')
plt.savefig("Reduction_comparison.pdf", format='pdf', bbox_inches='tight')
plt.close()
