import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import rc

P = [6, 8, 12, 16]
MaxProblemDimension = 60
MinProblemDimension = 10
ProblemForEachDimension = 10
persistency_list = []
all_pers_list = []
trivial_dec_list = []
var_reduction_list = []
scc_dec_list = []
probing_list = []
shannon_sizes = []
shannon_reduction = []
largest_qubo = []
tot_var = []
for prob in P:
    try:
        resfile = open("MaxCutStats_P{}.txt".format(prob), "r")
    except:
        print("Error in file creation\n")
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
                persistency_dict[num_of_vars].append(value)
            if param == "Num_trivial_Decomposition":
                trivial_dec_dict[num_of_vars].append(value)
            if param == "Num_SCC_Decomposition":
                scc_dec_dict[num_of_vars].append(value)
            if param == "All_persistencies":
                all_pers_dict[num_of_vars].append(value)
            if param == "Probing_persistencies":
                probing_dict[num_of_vars].append(value)
            if param == "Variable_reduction_percentage":
                var_reduction_dict[num_of_vars].append(value)
            if param == "Largest_subfunction":
                largest_q_dict[num_of_vars].append(value)
    for key in persistency_dict:
        pers_avg = np.average(persistency_dict[key])
        persistency_dict[key] = pers_avg
    for key in all_pers_dict:
        all_avg = np.average(all_pers_dict[key])
        all_pers_dict[key] = (all_avg * 100 / key)
    for key in trivial_dec_dict:
        triv_avg = np.average(trivial_dec_dict[key])
        trivial_dec_dict[key] = triv_avg
    for key in scc_dec_dict:
        scc_avg = np.average(scc_dec_dict[key])
        scc_dec_dict[key] = scc_avg
    for key in probing_dict:
        prob_avg = np.average(probing_dict[key])
        probing_dict[key] = (prob_avg * 100 / key)
    for key in var_reduction_dict:
        red_avg = np.average(var_reduction_dict[key])
        var_reduction_dict[key] = red_avg
    for key in largest_q_dict:
        q_avg = np.average(largest_q_dict[key])
        tot_var_dict[key] = ((key - q_avg) * 100 / key)
    persistency_list.append(persistency_dict)
    all_pers_list.append(all_pers_dict)
    trivial_dec_list.append(trivial_dec_dict)
    var_reduction_list.append(var_reduction_dict)
    scc_dec_list.append(scc_dec_dict)
    probing_list.append(probing_dict)
    largest_qubo.append(largest_q_dict)
    tot_var.append(tot_var_dict)

    resfile.close()

    # Collect results on Shannon decompositions
    shannon_dict = {}
    shannon_red_dict = {}
    for dim in range(MinProblemDimension, MaxProblemDimension + 1, 2):
        shannon_sum = 0
        reduction_sum = 0
        for n_prob in range(ProblemForEachDimension):
            try:
                sh_file = open("Results/MaxCut_" + format(dim) + "_nodes_P" + format(prob) + "_" + format(n_prob) + "_n_sol.txt", "r")
            except:
                print("Error in Shannon file opening\n")
                exit(1)

            for line in sh_file:
                splitted_line = line.split(' ')
                if splitted_line[0] == "Largest_Shannon_QUBO":
                    shannon_sum += (dim - int(splitted_line[1])) * 100 / dim
                    reduction_sum += (largest_q_dict[dim][n_prob] - int(splitted_line[1])) * 100 / dim
        shannon_dict[dim] = shannon_sum / ProblemForEachDimension
        shannon_red_dict[dim] = reduction_sum / ProblemForEachDimension
    shannon_sizes.append(shannon_dict)
    shannon_reduction.append(shannon_red_dict)

# Collect results on Dwave toolchain
try:
    dwave_file = open("Dwave_toolchain_res.txt", "r")
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
    elif len(splitted_line) == 1:
        n_nodes = int(splitted_line[0])
        dwave_list[i][n_nodes] = []
    else:
        dwave_list[i][n_nodes].append(int(splitted_line[0]))
for i in range(len(dwave_list)):
    for key in dwave_list[i]:
        dwave_list[i][key] = (np.average(dwave_list[i][key]) / key) * 100


# Persistency plot
for i in range(len(P)):
    plt.plot(persistency_list[i].keys(), persistency_list[i].values(), label="Density "+format(P[i])+"%", linewidth=2)
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Persistencies [%]', fontsize=15)
plt.title('Persistencies fixed', fontsize=20)
plt.legend()
plt.savefig("Persistencies.eps", format='eps', bbox_inches='tight')
plt.savefig("Persistencies.png", format='png', bbox_inches='tight')
plt.savefig("Persistencies.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Trivial decomposition plot
for i in range(len(P)):
    plt.plot(trivial_dec_list[i].keys(), trivial_dec_list[i].values(), label="Density "+format(P[i])+"%", linewidth=2)
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Number of decompositions', fontsize=15)
plt.title('Trivial decompositions', fontsize=20)
plt.legend()
plt.savefig("Trivial_decomposition.eps", format='eps', bbox_inches='tight')
plt.savefig("Trivial_decomposition.png", format='png', bbox_inches='tight')
plt.savefig("Trivial_decomposition.pdf", format='pdf', bbox_inches='tight')
plt.close()

# SCC decomposition plot
for i in range(len(P)):
    plt.plot(scc_dec_list[i].keys(), scc_dec_list[i].values(), label="Density "+format(P[i])+"%", linewidth=2)
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Number of decompositions', fontsize=15)
plt.title('SCC decompositions', fontsize=20)
plt.legend()
plt.savefig("SCC_decompositions.eps", format='eps', bbox_inches='tight')
plt.savefig("SCC_decompositions.png", format='png', bbox_inches='tight')
plt.savefig("SCC_decompositions.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Total variables reduction plot
for i in range(len(P)):
    plt.plot(var_reduction_list[i].keys(), var_reduction_list[i].values(), label="Density "+format(P[i])+"%", linewidth=2)
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Variable reduction [%]', fontsize=15)
plt.title('Total variable reduction', fontsize=20)
plt.legend()
plt.savefig("Var_reduction.eps", format='eps', bbox_inches='tight')
plt.savefig("Var_reduction.png", format='png', bbox_inches='tight')
plt.savefig("Var_reduction.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Total variables compared with dwave
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(1, 0, 4))
patches = []
for i in range(len(P)):
    plt.plot(var_reduction_list[i].keys(), var_reduction_list[i].values(), color=colors[i], linewidth=2)
    plt.plot(dwave_list[i].keys(), dwave_list[i].values(), color=colors[i], linewidth=2, linestyle='--')
    patches.append(mpatches.Patch(color=colors[i], label=format(P[i]) + "%"))
patches.append(mlines.Line2D([], [], color='0', linestyle='-', label="Qoolchain"))
patches.append(mlines.Line2D([], [], color='0', linestyle='--', label="Dwave"))
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Variable reduction [%]', fontsize=15)
plt.legend(title="Density", handles=patches)
plt.savefig("dwave_compare.eps", format='eps', bbox_inches='tight')
plt.savefig("dwave_compare.png", format='png', bbox_inches='tight')
plt.savefig("dwave_compare.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Shannon reduction compared with dwave
cmap = plt.get_cmap('magma')
colors = cmap(np.linspace(1, 0, 4))
patches = []
for i in range(len(P)):
    plt.plot(shannon_sizes[i].keys(), shannon_sizes[i].values(), color=colors[i], linewidth=2)
    plt.plot(dwave_list[i].keys(), dwave_list[i].values(), color=colors[i], linewidth=2, linestyle='--')
    patches.append(mpatches.Patch(color=colors[i], label=format(P[i])+"%"))
patches.append(mlines.Line2D([], [], color='0', linestyle='-', label="Qoolchain"))
patches.append(mlines.Line2D([], [], color='0', linestyle='--', label="Dwave"))
plt.xlabel('Number of nodes', fontsize=15)
plt.ylabel('Variable reduction [%]', fontsize=15)
plt.legend(title="Density", handles=patches)
plt.savefig("dwave-shannon.eps", format='eps', bbox_inches='tight')
plt.savefig("dwave-shannon.png", format='png', bbox_inches='tight')
plt.savefig("dwave-shannon.pdf", format='pdf', bbox_inches='tight')
plt.close()

# Plot of all techniques together
cmap = plt.get_cmap('plasma')
colors = cmap(np.linspace(1, 0, 4))
for i in range(len(P)):
    plt.plot(shannon_sizes[i].keys(), shannon_sizes[i].values(), color=colors[0], alpha=0.8, label='Shannon', linewidth=2)
    plt.fill_between(shannon_sizes[i].keys(), shannon_sizes[i].values(), color=colors[0], alpha=0.8)
    plt.plot(tot_var[i].keys(), tot_var[i].values(), color=colors[1], alpha=0.8, label='Decompositions', linewidth=2)
    plt.fill_between(tot_var[i].keys(), tot_var[i].values(), color=colors[1], alpha=0.8)
    plt.plot(probing_list[i].keys(), [x+y for x, y in zip(probing_list[i].values(), all_pers_list[i].values())], color=colors[2], alpha=0.8, label='Probing', linewidth=2)
    plt.fill_between(probing_list[i].keys(), [x+y for x, y in zip(probing_list[i].values(), all_pers_list[i].values())], color=colors[2], alpha=0.8)
    plt.plot(all_pers_list[i].keys(), all_pers_list[i].values(), color=colors[3], alpha=0.8, label='Persistencies', linewidth=2)
    plt.fill_between(all_pers_list[i].keys(), all_pers_list[i].values(), color=colors[3], alpha=0.8)
    plt.plot(dwave_list[i].keys(), dwave_list[i].values(), color='0', label='Dwave', linewidth=2, linestyle='--')
    plt.xlabel('Number of nodes', fontsize=15)
    plt.ylabel('Variable reduction [%]', fontsize=15)
    plt.legend()
    plt.savefig("Reduction_comparison_P" + format(P[i]) + ".eps", format='eps', bbox_inches='tight')
    plt.savefig("Reduction_comparison_P" + format(P[i]) + ".png", format='png', bbox_inches='tight')
    plt.savefig("Reduction_comparison_P" + format(P[i]) + ".pdf", format='pdf', bbox_inches='tight')
    plt.close()
