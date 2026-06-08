import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

NUMBER_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_LOCAL  = SCRIPT_DIR / "intact_vs_precut_overlay.png"
OUT_THESIS = SCRIPT_DIR.parent.parent / "Images" / "results" / "intact_vs_precut_overlay.png"


def parse_file(file_path: Path, label_suffix: str = "") -> dict:
    lines = [l.rstrip("\n") for l in file_path.read_text(encoding="utf-8").splitlines() if l.strip()]
    if not lines:
        return {}

    header_cells = lines[0].split("\t")
    job_columns = []
    for idx, cell in enumerate(header_cells):
        m = NUMBER_PATTERN.search(cell.strip())
        if m:
            key = str(float(m.group(0))) + label_suffix
            job_columns.append((key, idx))

    series = {key: {"d": [], "f": []} for key, _ in job_columns}
    for line in lines[1:]:
        cells = line.split("\t")
        for key, col_idx in job_columns:
            if col_idx + 1 >= len(cells):
                continue
            dr, fr = cells[col_idx].strip(), cells[col_idx + 1].strip()
            if not dr or not fr:
                continue
            try:
                series[key]["d"].append(float(dr))
                series[key]["f"].append(float(fr))
            except ValueError:
                pass
    return series


def clip_to_zero(d, f):
    d, f = np.array(d), np.array(f)
    mask = f >= 0
    out_d, out_f = [], []
    for i, (di, fi) in enumerate(zip(d, f)):
        if fi >= 0:
            out_d.append(di)
            out_f.append(fi)
        elif i > 0 and f[i - 1] > 0:
            x0, y0 = d[i - 1], f[i - 1]
            x_zero = x0 + (-y0) * (di - x0) / (fi - y0)
            out_d.append(x_zero)
            out_f.append(0.0)
            break
    return np.array(out_d), np.array(out_f)


def main():
    intact  = parse_file(SCRIPT_DIR / "f-d-all.txt")
    precut  = parse_file(SCRIPT_DIR / "f-d-precut.txt", label_suffix=" precut")

    designs = sorted({k for k in intact}, key=float)

    fig, axes = plt.subplots(3, 3, figsize=(14, 11))
    fig.subplots_adjust(hspace=0.42, wspace=0.35)

    for ax, design in zip(axes.flat, designs):
        key_i = design
        key_p = design + " precut"

        if key_i in intact:
            d, f = clip_to_zero(intact[key_i]["d"], intact[key_i]["f"])
            ax.plot(d, f, color="steelblue", linewidth=2.0, label="Intact")

        if key_p in precut:
            d, f = clip_to_zero(precut[key_p]["d"], precut[key_p]["f"])
            ax.plot(d, f, color="firebrick", linewidth=2.0,
                    linestyle="--", label="Precut")

        ax.set_title(f"$a = {design}$", fontsize=14)
        ax.set_xlabel("Displacement", fontsize=13)
        ax.set_ylabel("Force", fontsize=13)
        ax.tick_params(labelsize=11)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
        ax.legend(fontsize=11, frameon=False)

    for out in (OUT_LOCAL, OUT_THESIS):
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
