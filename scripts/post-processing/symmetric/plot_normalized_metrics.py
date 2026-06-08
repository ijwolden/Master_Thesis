import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

                                                                             
E_s     = 200.0                                 
sigma_y = 5.0                                  
W       = 9.0                                                   
L0      = 9.0                                                     
d_ref   = 1.0                                                              

                                                                             
NUMBER_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")
file_path = Path(__file__).parent / "f-d-all.txt"
lines = [l.rstrip("\n") for l in file_path.read_text(encoding="utf-8").splitlines()
         if l.strip()]

header_cells = lines[0].split("\t")
job_columns = []
for idx, cell in enumerate(header_cells):
    m = NUMBER_PATTERN.search(cell.strip())
    if m:
        job_columns.append((float(m.group(0)), idx))

series = {a: {"d": [], "f": []} for a, _ in job_columns}
for line in lines[1:]:
    cells = line.split("\t")
    for a, col_idx in job_columns:
        if col_idx + 1 >= len(cells):
            continue
        dr, fr = cells[col_idx].strip(), cells[col_idx + 1].strip()
        if not dr or not fr:
            continue
        try:
            series[a]["d"].append(float(dr))
            series[a]["f"].append(float(fr))
        except ValueError:
            pass

                                                                             
def compute_stiffness(d_arr, f_arr, fraction=0.10):
    idx_max = np.argmax(f_arr)
    mask = d_arr <= fraction * d_arr[idx_max]
    if mask.sum() < 3:
        return np.nan
    return np.polyfit(d_arr[mask], f_arr[mask], 1)[0]

designs = sorted(series.keys())
a_vals, E_norm, s_norm, w_norm = [], [], [], []

for a in designs:
    d_arr = np.array(series[a]["d"])
    f_arr = np.array(series[a]["f"])
    if len(f_arr) < 5:
        continue

    idx_max = int(np.argmax(f_arr))
    F_max   = f_arr[idx_max]
    k       = compute_stiffness(d_arr, f_arr)
    W_peak  = float(np.trapezoid(f_arr[: idx_max + 1], d_arr[: idx_max + 1]))

                                                       
    E_star_norm   = k / E_s                                                         
    sigma_star_norm = F_max / (sigma_y * W * d_ref)                                
    w_star_norm   = W_peak / (sigma_y * L0**2 * d_ref)                                  

    a_vals.append(a)
    E_norm.append(E_star_norm)
    s_norm.append(sigma_star_norm)
    w_norm.append(w_star_norm)

a_vals = np.array(a_vals)
E_norm = np.array(E_norm)
s_norm = np.array(s_norm)
w_norm = np.array(w_norm)

                                                                             
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

plot_data = [
    (E_norm, r"$E^*/E_s$",                    "Normalised stiffness"),
    (s_norm, r"$\sigma^*/\sigma_y$",           "Normalised peak stress"),
    (w_norm, r"$w/(\sigma_y L_0^2 \, d)$",    "Normalised energy density"),
]
markers = ["o", "s", "^"]

for ax, (y, ylabel, title), mk in zip(axes, plot_data, markers):
    ax.plot(a_vals, y, marker=mk, linewidth=1.6, color="steelblue",
            markersize=6, markerfacecolor="white", markeredgewidth=1.5)
    ax.set_xlabel("$a$", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_title(title, fontsize=15)
    ax.set_xticks(a_vals)
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    ax.tick_params(labelsize=13)

plt.tight_layout(pad=1.5)

                                                                             
out_local = Path(__file__).parent / "normalized_metrics_Es_sigma.png"
out_thesis = Path(__file__).parent.parent.parent / "Images" / "results" / "normalized_metrics_Es_sigma.png"

for out in (out_local, out_thesis):
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")

plt.show()
