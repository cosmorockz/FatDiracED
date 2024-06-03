import re
import glob
import numpy as np
import matplotlib.pyplot as plt

cutoff = 1
Q = 2
kinConst = 1e+2
directory = f'./Simulation_Data/cutoff_{cutoff}/Q_{Q}'

files = glob.glob(f'{directory}/*')

def lists_to_dict(LL, EE):
    spec_dict = {}
    for i, LzTot in enumerate(LL):
        if LzTot not in spec_dict.keys():
            spec_dict[LzTot] = []
        spec_dict[LzTot].append(EE[i])

    return spec_dict



spec_dict = {}
for filename in files:
    match = re.search(r'spec_(-?\d+\.?\d*)', filename)
    if match:
        intStr = float(match.group(1))
        LL, spec = np.loadtxt(filename, usecols=(0,1), unpack=True)
        dict2 = lists_to_dict(LL, spec)
        for key, vals in dict2.items():
            spec_dict.setdefault(key, {}).setdefault(intStr, [])
            spec_dict[key][intStr] = vals

# Plotting
markers = ['1', '2', '+', '_', 'x']  # Different markers for the lowest five values
colors = plt.cm.rainbow(np.linspace(0, 1, len(spec_dict)))  # Different colors for different keys

plt.figure(figsize=(10, 6))

for idx, key in enumerate(spec_dict.keys()):
    if key > -1e-8:
        color = colors[idx]
        markers_copy = markers[:]  # Make a copy of markers for each key
        for intStr, vals in spec_dict[key].items():
            vals_sorted = np.array(sorted(vals)[:5])/kinConst  # Sort and get the lowest five values
            marker = markers_copy.pop(0) if markers_copy else 'o'  # Use 'o' if markers_copy is empty
            for i in range(len(vals_sorted)):
                plt.scatter(intStr, vals_sorted[i], color=color, marker=markers[i], label=intStr, s=40, linewidths=0.5)
            # plt.scatter([intStr] * len(vals_sorted), vals_sorted, marker=marker, color=color, label=f'Key: {key}')


# Place annotations outside the plot area
for idx, (key, color) in enumerate(zip(spec_dict.keys(), colors)):
    if key > -1e-8:
        plt.text(1.05, 0.95 - idx*0.05, f'L=: {key}', color=color, transform=plt.gca().transAxes,
                 fontsize=10, verticalalignment='top', horizontalalignment='left')
plt.xlabel('$U_{int}$/$U_{kin}$')
plt.ylabel('$E$')
plt.title('Lowest five levels in each $L$ sector')
# plt.legend()
plt.grid(True)
plt.savefig(f'./Plots/cutoff_{cutoff}_Q_{Q}.pdf', bbox_inches='tight', dpi=300)
plt.grid(True)
plt.show()






















