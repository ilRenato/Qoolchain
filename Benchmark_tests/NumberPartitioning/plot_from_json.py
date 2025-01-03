import json
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import rc
import numpy as np

matplotlib.rcParams.update(matplotlib.rcParamsDefault)
rc('text', usetex=True)
cmap = plt.get_cmap('RdYlGn')
colors = cmap(np.linspace(1, 0, 4))

with open("var_results.json", "r") as var_file:
    data = json.load(var_file)

    densities = list(data['Persistencies'].keys())
    nodes_max = int(float(list(data['Persistencies'][densities[0]].keys())[-1]))
    x = [int(float(num)) for num in data['Persistencies'][densities[0]].keys()]

    # Persistency plot
    i = 0
    for d in data['Persistencies'].keys():
        plt.plot(x, data['Persistencies'][d].values(), label=format(d), color=colors[i], linewidth=2)
        i += 1
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Persistencies $[\%]$', fontsize=16)
    plt.legend(title="Largest Numbernumber")
    # plt.savefig("Plots\\Persistencies.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\Persistencies.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\Persistencies.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Trivial decomposition plot
    i = 0
    for d in data['Trivial_decomposition'].keys():
        plt.plot(x, data['Trivial_decomposition'][d].values(), label=format(d), color=colors[i], linewidth=2)
        i += 1
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Number of decompositions', fontsize=16)
    plt.legend(title="Largest Number")
    # plt.savefig("Plots\\Trivial_decomposition.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\Trivial_decomposition.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\Trivial_decomposition.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # SCC decomposition plot
    i = 0
    for d in data['CSCC_decomposition'].keys():
        plt.plot(x, data['CSCC_decomposition'][d].values(), label=format(d), color=colors[i], linewidth=2)
        i += 1
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Number of decompositions', fontsize=16)
    plt.legend(title="Largest Number", loc='upper left')
    # plt.savefig("Plots\\SCC_decompositions.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\SCC_decompositions.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\SCC_decompositions.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Total variables reduction plot
    i = 0
    for d in data['Var_reduction'].keys():
        plt.plot(x, data['Var_reduction'][d].values(), label=format(d), color=colors[i], linewidth=2)
        i += 1
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Variable reduction $[\%]$', fontsize=16)
    plt.legend(title="Largest Number")
    # plt.savefig("Plots\\Var_reduction.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\Var_reduction.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\Var_reduction.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Total variables compared with dwave
    patches = []
    i = 0
    for d in data['Var_reduction'].keys():
        plt.plot(x, data['Var_reduction'][d].values(), color=colors[i], linewidth=2)
        plt.plot(x, data['Dwave'][d].values(), color=colors[i], linewidth=2, linestyle='--')
        patches.append(mpatches.Patch(color=colors[i], label=format(d)))
        i += 1
    patches.append(mlines.Line2D([], [], color='0', linestyle='-', label="Qoolchain"))
    patches.append(mlines.Line2D([], [], color='0', linestyle='--', label="Dwave"))
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Variable reduction $[\%]$', fontsize=16)
    plt.legend(title="Largest Number", handles=patches, loc='lower left')
    # plt.savefig("Plots\\dwave_compare.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\dwave_compare.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\dwave_compare.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Shannon reduction compared with dwave
    patches = []
    i = 0
    for d in data['Shannon_sizes'].keys():
        plt.plot(x, data['Shannon_sizes'][d].values(), color=colors[i], linewidth=2)
        plt.plot(x, data['Dwave'][d].values(), color=colors[i], linewidth=2, linestyle='--')
        patches.append(mpatches.Patch(color=colors[i], label=format(d)))
        i += 1
    patches.append(mlines.Line2D([], [], color='0', linestyle='-', label="Qoolchain"))
    patches.append(mlines.Line2D([], [], color='0', linestyle='--', label="Dwave"))
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Variable reduction $[\%]$', fontsize=16)
    plt.legend(title="Largest Number", handles=patches, loc='lower left')
    # plt.savefig("Plots\\dwave-shannon.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\dwave-shannon.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\dwave-shannon.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Plot of all techniques together
    cmap = plt.get_cmap('plasma')
    colors = cmap(np.linspace(1, 0, 4))
    for d in data['Shannon_sizes'].keys():
        plt.plot(x, data['Shannon_sizes'][d].values(), color=colors[0], alpha=0.8, label='Shannon', linewidth=2)
        plt.fill_between(x, data['Shannon_sizes'][d].values(), color=colors[0], alpha=0.8)
        plt.plot(x, data['Total_vars'][d].values(), color=colors[1], alpha=0.8, label='Decompositions', linewidth=2)
        plt.fill_between(x, data['Total_vars'][d].values(), color=colors[1], alpha=0.8)
        plt.plot(x, [x + y for x, y in zip(data['Probing'][d].values(), data['All_persistencies'][d].values())], color=colors[2], alpha=0.8, label='Probing', linewidth=2)
        plt.fill_between(x, [x + y for x, y in zip(data['Probing'][d].values(), data['All_persistencies'][d].values())], color=colors[2], alpha=0.8)
        plt.plot(x, data['All_persistencies'][d].values(), color=colors[3], alpha=0.8, label='Persistencies', linewidth=2)
        plt.fill_between(x, data['All_persistencies'][d].values(), color=colors[3], alpha=0.8)
        plt.plot(x, data['Dwave'][d].values(), color='0', label='Dwave', linewidth=2, linestyle='--')
        plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
        plt.xlabel('Set size', fontsize=16)
        plt.ylabel('Variable reduction $[\%]$', fontsize=16)
        plt.legend()
        # plt.savefig("Plots\\Reduction_comparison_P" + format(d) + ".eps", format='eps', bbox_inches='tight')
        plt.savefig("Plots\\Reduction_comparison_P" + format(d) + ".png", format='png', bbox_inches='tight')
        plt.savefig("Plots\\Reduction_comparison_P" + format(d) + ".pdf", format='pdf', bbox_inches='tight')
        plt.close()

with open("bounds_results.json", "r") as bounds_file:
    data = json.load(bounds_file)

    densities = list(data['Error'].keys())
    nodes_max = int(float(list(data['Error'][densities[0]].keys())[-1]))
    x = [int(float(num)) for num in data['Error'][densities[0]].keys()]

    cmap = plt.get_cmap('RdYlGn')
    colors = cmap(np.linspace(1, 0, 4))
    i = 0
    for d in data['Error'].keys():
        plt.plot(x, data['Error'][d].values(), label=format(d), linewidth=2, color=colors[i])
        i += 1
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Error over estimated range $[\%]$', fontsize=16)
    plt.legend(title="Largest Number")
    # plt.savefig("Plots\\Lower_bound.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\Lower_bound.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\Lower_bound.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Other methods plot
    patches = []
    i = 0
    for d in data['Error'].keys():
        plt.plot(x, data['Error'][d].values(), linewidth=2, linestyle='solid', color=colors[i])
        plt.plot(x, data['Naive'][d].values(), linewidth=2, linestyle='dashed', color=colors[i])
        plt.plot(x, data['Posneg'][d].values(), linewidth=2, linestyle='dotted', color=colors[i])
        patches.append(mpatches.Patch(color=colors[i], label=format(d)))
        i += 1
    patches.append(mlines.Line2D([], [], color='0', linestyle='solid', label="Qoolchain"))
    patches.append(mlines.Line2D([], [], color='0', linestyle='dashed', label="naive"))
    patches.append(mlines.Line2D([], [], color='0', linestyle='dotted', label="posneg"))
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Error over estimated range $[\%]$', fontsize=16)
    plt.legend(title="Largest Number", handles=patches)
    # plt.savefig("Plots\\Lower_bound_compare.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\Lower_bound_compare.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\Lower_bound_compare.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Dwave comparison
    patches = []
    i = 0
    for d in data['Error'].keys():
        plt.plot(x, data['Error'][d].values(), linewidth=2, color=colors[i])
        plt.plot(x, data['Dwave'][d].values(), linewidth=2, linestyle='dashed', color=colors[i])
        patches.append(mpatches.Patch(color=colors[i], label=format(d)))
        i += 1
    patches.append(mlines.Line2D([], [], color='0', linestyle='solid', label="Qoolchain"))
    patches.append(mlines.Line2D([], [], color='0', linestyle='dashed', label="Dwave"))
    plt.xticks(ticks=np.arange(10, nodes_max, 20), labels=np.arange(10, nodes_max, 20))
    plt.xlabel('Set size', fontsize=16)
    plt.ylabel('Error over estimated range $[\%]$', fontsize=16)
    plt.legend(title="Largest Number", handles=patches)
    # plt.savefig("Plots\\Range_dwave.eps", format='eps', bbox_inches='tight')
    plt.savefig("Plots\\Range_dwave.png", format='png', bbox_inches='tight')
    plt.savefig("Plots\\Range_dwave.pdf", format='pdf', bbox_inches='tight')
    plt.close()
