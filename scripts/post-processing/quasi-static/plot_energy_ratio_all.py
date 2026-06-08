import re
from pathlib import Path
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

HERE = Path(__file__).parent
QS   = 0.05
LW   = 1.6

                                                                               
def parse_t1(path: Path):
    sections = {}
    current_label = None
    t_buf, r_buf = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = re.match(r"#\s*(\S+)", line)
            if m:
                if current_label is not None and t_buf:
                    sections[current_label] = (np.array(t_buf), np.array(r_buf))
                current_label = m.group(1)
                t_buf, r_buf = [], []
            else:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        t_buf.append(float(parts[0]))
                        r_buf.append(float(parts[1]))
                    except ValueError:
                        pass
    if current_label is not None and t_buf:
        sections[current_label] = (np.array(t_buf), np.array(r_buf))

                               
    for key in ("1", "1.0"):
        if key in sections:
            t, r = sections[key]
            return t / t.max(), r
    raise KeyError(f"No t=1 section found in {path}. Keys: {list(sections)}")


def make_figure(series: list, outpath: Path, figsize=(7, 4.5)):
    all_r = np.concatenate([r for _, _, _, r in series])
    ymax = 1.05 if all_r.max() > 0.1 else 0.06
    fig, ax = plt.subplots(figsize=figsize)

    ax.axhspan(0, QS, color="#c8e6c9", alpha=0.45, zorder=0)
    ax.axhline(QS, color="#2e7d32", linewidth=1.0, linestyle='--', zorder=1)

    for label, color, t, r in series:
        ax.plot(t, r, color=color, linestyle='-', linewidth=LW,
                label=label, zorder=2)

    ax.set_ylim(0, ymax)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel(r"Normalised time  $t\,/\,t_\mathrm{total}$", fontsize=10)
    ax.set_ylabel("ALLKE / ALLIE", fontsize=10)
    ax.grid(True, which='both', color='#e0e0e0', linewidth=0.5, zorder=0)
    ax.tick_params(labelsize=9)

    legend_handles = [
        Line2D([0], [0], color=c, lw=LW, ls='-', label=lbl)
        for lbl, c, *_ in series
    ]
    qs_patch = mpatches.Patch(facecolor='#c8e6c9', edgecolor='#2e7d32',
                              label='QS zone (< 5%)')
    ax.legend(handles=legend_handles + [qs_patch],
              fontsize=8.5, frameon=True, framealpha=0.9,
              edgecolor='#cccccc', loc='upper right',
              handlelength=2.2, labelspacing=0.4)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath.name}")


baseline_series = [
    ("Baseline A", "#4878CF",
     *parse_t1(HERE / "Baseline/A-energy.txt")),
    ("Baseline B", "#D65F5F",
     *parse_t1(HERE / "Baseline/b-energy.txt")),
]
make_figure(baseline_series, HERE / "allke_allie_baseline.png")

                                                                                
g1_series = [
    ("G1-A", "#4878CF", *parse_t1(HERE / "G1/a.txt")),
    ("G1-B", "#D65F5F", *parse_t1(HERE / "G1/b.txt")),
    ("G1-C", "#3DA35D", *parse_t1(HERE / "G1/c.txt")),
]
make_figure(g1_series, HERE / "allke_allie_g1.png")

                                                                                
g2_series = [
    ("G2-A", "#4878CF", *parse_t1(HERE / "G2/a-energy.txt")),
    ("G2-B", "#D65F5F", *parse_t1(HERE / "G2/b-energy.txt")),
]
make_figure(g2_series, HERE / "allke_allie_g2.png", figsize=(6, 4.5))

                                                                                
g3_series = [
    ("G3", "#4878CF", *parse_t1(HERE / "G3/g3_energy.txt")),
]
make_figure(g3_series, HERE / "allke_allie_g3.png", figsize=(6, 4.5))

                                                                               
sym_data = {}
sym_path = HERE.parent / "Symmetric/energy_all.txt"
current = None
t_buf, r_buf = [], []
with open(sym_path) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        m = re.match(r"#(\S+)", line)
        if m:
            if current is not None and t_buf:
                sym_data[current] = (np.array(t_buf), np.array(r_buf))
            current = m.group(1)
            t_buf, r_buf = [], []
        else:
            parts = line.split()
            if len(parts) == 2:
                try:
                    t_buf.append(float(parts[0]))
                    r_buf.append(float(parts[1]))
                except ValueError:
                    pass
if current is not None and t_buf:
    sym_data[current] = (np.array(t_buf), np.array(r_buf))

                          
for k in sym_data:
    t, r = sym_data[k]
    sym_data[k] = (t / t.max(), r)

                                                 
ALL_COLORS = {
    "Baseline A": "#4878CF", "Baseline B": "#1a3f6f",
    "G1-A": "#D65F5F",       "G1-B": "#a02020",    "G1-C": "#f0a0a0",
    "G2-A": "#3DA35D",       "G2-B": "#1a5c30",
    "G3":   "#E8882A",
}
SYM_COLOR = "#8E6BB0"

all_series = []
for lbl, col in ALL_COLORS.items():
    family, var = lbl.split("-")[0], lbl.split("-")[-1] if "-" in lbl else ""
    if family == "Baseline":
        key = "A" if var == "A" else "b"
        path = HERE / f"Baseline/{key}-energy.txt"
        t, r = parse_t1(path)
    elif family == "G1":
        t, r = parse_t1(HERE / f"G1/{var.lower()}.txt")
    elif family == "G2":
        t, r = parse_t1(HERE / f"G2/{var.lower()}-energy.txt")
    elif family == "G3":
        t, r = parse_t1(HERE / "G3/g3_energy.txt")
    all_series.append((lbl, col, t, r))

                       
for des in ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9']:
    if des in sym_data:
        t, r = sym_data[des]
        all_series.append((f"Sym {des}", SYM_COLOR, t, r))

                                                                 
fig, ax = plt.subplots(figsize=(9, 5))

ax.axhspan(0, QS, color="#c8e6c9", alpha=0.45, zorder=0)
ax.axhline(QS, color="#2e7d32", linewidth=1.0, linestyle='--', zorder=1)

for label, color, t, r in all_series:
    alpha = 0.55 if "Sym" in label else 1.0
    ax.plot(t, r, color=color, linestyle='-', linewidth=LW if "Sym" not in label else 1.0,
            alpha=alpha, zorder=2)

ax.set_ylim(0, 1.05)
ax.set_xlim(0, 1.0)
ax.set_xlabel(r"Normalised time  $t\,/\,t_\mathrm{total}$", fontsize=10)
ax.set_ylabel("ALLKE / ALLIE", fontsize=10)
ax.grid(True, which='both', color='#e0e0e0', linewidth=0.5, zorder=0)
ax.tick_params(labelsize=9)

                        
legend_handles = [
    Line2D([0],[0], color="#4878CF", lw=LW, label="Baseline A"),
    Line2D([0],[0], color="#1a3f6f", lw=LW, label="Baseline B"),
    Line2D([0],[0], color="#D65F5F", lw=LW, label="G1-A"),
    Line2D([0],[0], color="#a02020", lw=LW, label="G1-B"),
    Line2D([0],[0], color="#f0a0a0", lw=LW, label="G1-C"),
    Line2D([0],[0], color="#3DA35D", lw=LW, label="G2-A"),
    Line2D([0],[0], color="#1a5c30", lw=LW, label="G2-B"),
    Line2D([0],[0], color="#E8882A", lw=LW, label="G3"),
    Line2D([0],[0], color=SYM_COLOR, lw=1.0, alpha=0.55, label="Symmetric (0.1–0.9)"),
    mpatches.Patch(facecolor='#c8e6c9', edgecolor='#2e7d32', label='QS zone (< 5%)'),
]
ax.legend(handles=legend_handles, fontsize=8.5, frameon=True, framealpha=0.9,
          edgecolor='#cccccc', loc='upper right', ncol=2,
          handlelength=2.2, labelspacing=0.4)

plt.tight_layout()
plt.savefig(HERE / "allke_allie_all.png", dpi=200, bbox_inches='tight')
plt.close()
print("Saved: allke_allie_all.png")
