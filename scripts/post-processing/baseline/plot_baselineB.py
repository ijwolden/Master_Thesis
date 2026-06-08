import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _to_float(value: str):
    value = value.strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_baseline_b(file_path: Path):
    with file_path.open(newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))

    if not rows:
        raise ValueError(f"No data found in {file_path}")

    header = rows[0]
    block_starts = [0, 6, 12, 18, 24]

    datasets = {}
    for start in block_starts:
        job_label = header[start + 2].strip()
        fd = []
        energy = []

        for row in rows[1:]:
            if len(row) <= start + 4:
                continue

            disp = _to_float(row[start])
            force = _to_float(row[start + 1])
            time = _to_float(row[start + 3])
            ratio = _to_float(row[start + 4])

            if disp is not None and force is not None:
                fd.append((disp, force))
            if time is not None and ratio is not None:
                energy.append((time, ratio))

        datasets[job_label] = {
            "fd": np.array(fd, dtype=float),
            "energy": np.array(energy, dtype=float),
        }

    return datasets


def force_curve_to_zero(data: np.ndarray) -> np.ndarray:
    if data.size == 0:
        return data

    data = data[data[:, 0] >= 0]
    if data.size == 0:
        return data

    data = data[np.argsort(data[:, 0])]

                                                                        
    for i in range(1, len(data)):
        x0, y0 = data[i - 1]
        x1, y1 = data[i]
        if y0 > 0 and y1 <= 0 and x1 != x0 and y1 != y0:
            x_zero = x0 + (-y0) * (x1 - x0) / (y1 - y0)
            return np.vstack([data[:i], [x_zero, 0.0]])

    nonnegative = data[data[:, 1] >= 0]
    if len(nonnegative) < 2:
        return nonnegative

    x0, y0 = nonnegative[-2]
    x1, y1 = nonnegative[-1]
    if x1 == x0:
        return nonnegative

    slope = (y1 - y0) / (x1 - x0)
    if slope < 0:
        x_zero = x1 - y1 / slope
        if np.isfinite(x_zero) and x_zero > x1:
            return np.vstack([nonnegative, [x_zero, 0.0]])

    return nonnegative


def main():
    base_dir = Path(__file__).parent
    data_file = base_dir / "baselineB.txt"

    datasets = parse_baseline_b(data_file)
    job_order = ["0.01", "0.1", "1", "1.5", "2"]

                                         
    plt.figure(figsize=(10, 6))
    for label in job_order:
        if label not in datasets:
            continue
        data = datasets[label]["fd"]
        if data.size == 0:
            continue
        data = force_curve_to_zero(data)
        if data.size == 0:
            continue
        plt.plot(data[:, 0], data[:, 1], linewidth=2, label=f"Job {label}")

    plt.xlabel("Displacement")
    plt.ylabel("Force")
    plt.title("Baseline B: Force-Displacement")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    fd_out = base_dir / "baselineB_force_displacement.png"
    plt.savefig(fd_out, dpi=300, bbox_inches="tight")

                                   
    plt.figure(figsize=(10, 6))
    for label in job_order:
        if label not in datasets:
            continue
        data = datasets[label]["energy"]
        if data.size == 0:
            continue
        plt.plot(data[:, 0], data[:, 1], linewidth=2, label=f"Job {label}")

    plt.xlabel("Time")
    plt.ylabel("ALLKE/ALLIE")
    plt.title("Baseline B: Energy Ratio vs Time")
    plt.yscale("log")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    e_out = base_dir / "baselineB_energy_ratio.png"
    plt.savefig(e_out, dpi=300, bbox_inches="tight")

    print(f"Saved: {fd_out.name}")
    print(f"Saved: {e_out.name}")


if __name__ == "__main__":
    main()
