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

with open("bounds_results.json", "r") as bounds_file:
    data = json.load(bounds_file)

    x = np.arange(len(data['Error']))
    plt.bar(x, data['Error'].values())
    plt.xlabel('Benchmark', fontsize=16)
    plt.xticks(x, data['Error'].keys(), rotation=33)
    plt.ylabel('Error over estimated range $[\%]$', fontsize=16)
    plt.savefig("Lower_bound.png", format='png', bbox_inches='tight')
    plt.savefig("Lower_bound.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Other methods plot
    colors = cmap(np.linspace(1, 0, 4))
    width = 0.25
    plt.bar(x, data['Error'].values(), label='Qoolchain', color=colors[0], width=width)
    plt.bar(x + width, data['Posneg'].values(), label='posneg', color=colors[2], width=width)
    plt.bar(x + 2 * width, data['Naive'].values(), label='naive', color=colors[3], width=width)
    plt.xlabel('Benchmark', fontsize=16)
    plt.xticks(x + width, data['Error'].keys(), rotation=33)
    plt.ylabel('Error over estimated range $[\%]$', fontsize=16)
    plt.legend(loc='upper left')
    plt.savefig("Lower_bound_compare.png", format='png', bbox_inches='tight')
    plt.savefig("Lower_bound_compare.pdf", format='pdf', bbox_inches='tight')
    plt.close()

    # Dwave comparison
    colors = cmap(np.linspace(1, 0, 2))
    plt.bar(x, data['Error'].values(), label='Qoolchain', color=colors[0], width=width)
    plt.bar(x + width, data['Dwave'].values(), label='dwave', color=colors[1], width=width)
    plt.xlabel('Benchmark', fontsize=16)
    plt.xticks(x + width / 2, data['Error'].keys(), rotation=33)
    plt.ylabel('Error over estimated range $[\%]$', fontsize=16)
    plt.legend(loc='upper left')
    plt.savefig("Range_dwave.png", format='png', bbox_inches='tight')
    plt.savefig("Range_dwave.pdf", format='pdf', bbox_inches='tight')
    plt.close()
