import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

                                                                                
plt.rcParams.update({
    "font.family": "serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
})

                                                                                
data = {}
current = None
disp, ratio = [], []

with open("energy_all.txt") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        m = re.match(r"#(\S+)", line)
        if m:
            if current is not None and disp:
                data[current] = (np.array(disp), np.array(ratio))
            current = m.group(1)
            disp, ratio = [], []
        else:
            parts = line.split()
            if len(parts) == 2:
                try:
                    disp.append(float(parts[0]))
                    ratio.append(float(parts[1]))
                except ValueError:
                    pass
if current is not None and disp:
    data[current] = (np.array(disp), np.array(ratio))

                                                         
for key in data:
    t, r = data[key]
    data[key] = (t / t.max(), r)

designs = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']

PAIR_COLORS = {
    '0.1': '#4878CF', '0.9': '#4878CF',
    '0.2': '#D65F5F', '0.8': '#D65F5F',
    '0.3': '#3DA35D', '0.7': '#3DA35D',
    '0.4': '#E8882A', '0.6': '#E8882A',
    '0.5': '#8E6BB0',
}

LW = 1.6
QS = 0.05

                                                                                
fig, ax = plt.subplots(figsize=(7, 4.5))

ax.axhspan(0, QS, color="#c8e6c9", alpha=0.45, zorder=0)
ax.axhline(QS, color="#2e7d32", linewidth=1.0, linestyle='--', zorder=1)

for des in designs:
    if des not in data:
        continue
    t, r = data[des]
    ax.plot(t, r,
            color=PAIR_COLORS[des],
            linestyle='-',
            linewidth=LW,
            zorder=2)

ax.set_ylim(0, 1.05)
ax.set_xlim(0, 1.0)
ax.set_xlabel(r"Normalised time  $t\,/\,t_\mathrm{total}$", fontsize=10)
ax.set_ylabel("ALLKE / ALLIE", fontsize=10)
ax.grid(True, which='both', color='#e0e0e0', linewidth=0.5, zorder=0)
ax.tick_params(labelsize=9)

                                                                                
pair_handles = [
    Line2D([0], [0], color='#4878CF', lw=LW, ls='-', label='p = 0.1,  0.9'),
    Line2D([0], [0], color='#D65F5F', lw=LW, ls='-', label='p = 0.2,  0.8'),
    Line2D([0], [0], color='#3DA35D', lw=LW, ls='-', label='p = 0.3,  0.7'),
    Line2D([0], [0], color='#E8882A', lw=LW, ls='-', label='p = 0.4,  0.6'),
    Line2D([0], [0], color='#8E6BB0', lw=LW, ls='-', label='p = 0.5'),
]
qs_handle = mpatches.Patch(facecolor='#c8e6c9', edgecolor='#2e7d32',
                            label='QS zone (< 5%)')

ax.legend(
    handles=pair_handles + [qs_handle],
    fontsize=8.5, frameon=True, framealpha=0.9,
    edgecolor='#cccccc', loc='upper right',
    handlelength=2.2, labelspacing=0.4,
)

plt.savefig("energy_ratio_quasi_static.png", dpi=200, bbox_inches='tight')
print("Saved: energy_ratio_quasi_static.png")
plt.show()
