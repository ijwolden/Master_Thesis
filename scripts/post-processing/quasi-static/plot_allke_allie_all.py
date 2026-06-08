from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

ENERGY_DIR = Path(__file__).parent
METHODS_OUT = Path(__file__).parent
THRESHOLD = 0.05                           

FAMILIES = ["Baseline", "Symmetric", "G1", "G2", "G3"]


def parse_energy_file(path: Path) -> dict:
    sections = {}
    current_label = None
    rows = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if current_label is not None and rows:
                    sections[current_label] = np.array(rows, dtype=float)
                current_label = line.lstrip("#").strip()
                rows = []
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        rows.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        pass

    if current_label is not None and rows:
        sections[current_label] = np.array(rows, dtype=float)

    return sections


def _fd_to_zero(data: np.ndarray) -> np.ndarray:
    data = data[data[:, 0] >= 0]
    if data.size == 0:
        return data
    data = data[np.argsort(data[:, 0])]
    for i in range(1, len(data)):
        x0, y0 = data[i - 1]
        x1, y1 = data[i]
        if y0 > 0 and y1 <= 0 and x1 != x0:
            x_zero = x0 + (-y0) * (x1 - x0) / (y1 - y0)
            return np.vstack([data[:i], [x_zero, 0.0]])
    nonneg = data[data[:, 1] >= 0]
    if len(nonneg) >= 2:
        x0, y0 = nonneg[-2]
        x1, y1 = nonneg[-1]
        if x1 != x0:
            slope = (y1 - y0) / (x1 - x0)
            if slope < 0:
                x_zero = x1 - y1 / slope
                if np.isfinite(x_zero) and x_zero > x1:
                    return np.vstack([nonneg, [x_zero, 0.0]])
    return nonneg


def plot_fd_design(ax, sections: dict, design_name: str, add_legend: bool = True):
    try:
        keys_sorted = sorted(sections.keys(), key=float)
    except ValueError:
        keys_sorted = list(sections.keys())
    n = len(keys_sorted)
    colors = cm.tab10.colors[:n] if n <= 10 else [cm.tab10(i / n) for i in range(n)]
    for key, color in zip(keys_sorted, colors):
        data = _fd_to_zero(sections[key])
        if data.size == 0:
            continue
        ax.plot(data[:, 0], data[:, 1], color=color, linewidth=1.4,
                label=f"$t = {key}$ s")
    ax.set_title(design_name, fontsize=9)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Displacement (mm)", fontsize=8)
    ax.set_ylabel("Force (kN)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
    if add_legend:
        ax.legend(fontsize=7, loc="upper right")


def plot_fd_family(family: str):
    family_dir = ENERGY_DIR / family
    files = sorted(family_dir.glob("*-fd.txt"))
    if not files:
        return                                   
    n = len(files)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.5 * nrows),
                             squeeze=False)
    for idx, filepath in enumerate(files):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sections = parse_energy_file(filepath)               
        design_name = filepath.stem.replace("-fd", "")
        plot_fd_design(ax, sections, design_name)
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)
    fig.suptitle(f"F-D — {family}", fontsize=11, y=1.01)
    fig.tight_layout()
    out = ENERGY_DIR / f"fd_convergence_{family.lower()}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
                                                                                

def plot_design(ax, sections: dict, design_name: str, add_legend: bool = True):
    try:
        keys_sorted = sorted(sections.keys(), key=float)
    except ValueError:
        keys_sorted = list(sections.keys())

    n = len(keys_sorted)
    colors = cm.tab10.colors[:n] if n <= 10 else [cm.tab10(i / n) for i in range(n)]

    for (key, color) in zip(keys_sorted, colors):
        data = sections[key]
        try:
            t_total = float(key)
        except ValueError:
            t_total = 1.0
        t_norm = data[:, 0] / t_total
        ratio_pct = data[:, 1] * 100
                                                                            
        mask = ratio_pct > 0
        ax.semilogy(t_norm[mask], ratio_pct[mask], color=color, linewidth=1.4,
                    label=f"$t = {key}$ s")

    ax.axhline(THRESHOLD * 100, color="black", linestyle="--",
               linewidth=1.0, label="5% threshold")

    ax.set_title(design_name, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel("$t / t_{\\mathrm{total}}$", fontsize=8)
    ax.set_ylabel("ALLKE/ALLIE (%)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.6)

    if add_legend:
        ax.legend(fontsize=7, loc="upper right")


def plot_family(family: str):
    family_dir = ENERGY_DIR / family
    files = sorted(f for f in family_dir.glob("*.txt") if not f.name.endswith("-fd.txt"))

    if not files:
        print(f"  No energy files found in {family_dir} — skipping.")
        return

    n = len(files)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.5 * nrows),
                             squeeze=False)

    for idx, filepath in enumerate(files):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sections = parse_energy_file(filepath)
        design_name = filepath.stem.replace("-energy", "")
        plot_design(ax, sections, design_name)

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(f"Quasi-static verification — {family}", fontsize=11, y=1.01)
    fig.tight_layout()

    out = METHODS_OUT / f"allke_allie_{family.lower()}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_combined():
    all_files = []
    for family in FAMILIES:
        family_dir = ENERGY_DIR / family
        for f in sorted(family_dir.glob("*.txt")):
            if not f.name.endswith("-fd.txt"):
                all_files.append((family, f))

    if not all_files:
        print("No energy files found anywhere — skipping combined figure.")
        return

    n = len(all_files)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.0 * ncols, 3.2 * nrows),
                             squeeze=False)

    for idx, (family, filepath) in enumerate(all_files):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sections = parse_energy_file(filepath)
        design_name = f"{family} {filepath.stem.replace('-energy', '')}"
        plot_design(ax, sections, design_name, add_legend=(idx == 0))

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle("ALLKE/ALLIE vs. time",
                 fontsize=12, y=1.01)
    fig.tight_layout()

    out = METHODS_OUT / "allke_allie_all.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def main():
    ENERGY_DIR.mkdir(parents=True, exist_ok=True)

    for family in FAMILIES:
        print(f"Processing {family}...")
        plot_family(family)
        plot_fd_family(family)

    print("Generating combined figure...")
    plot_combined()
    print("Done.")


if __name__ == "__main__":
    main()
