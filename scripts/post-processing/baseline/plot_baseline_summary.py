from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_baseline_sections(file_path: Path):
    datasets = {}
    current_label = None
    current_data = []

    for raw_line in file_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            if current_label is not None and current_data:
                datasets[current_label] = np.array(current_data, dtype=float)
            current_label = line.lower().replace(" ", "")
            current_data = []
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        try:
            displacement = float(parts[0])
            force = float(parts[1])
        except ValueError:
            continue

        current_data.append((displacement, force))

    if current_label is not None and current_data:
        datasets[current_label] = np.array(current_data, dtype=float)

    return datasets


def force_curve_to_zero(data: np.ndarray) -> np.ndarray:
    if data.size == 0:
        return data

    data = data[data[:, 0] >= 0]
    if data.size == 0:
        return data

    data = data[np.argsort(data[:, 0])]

    for index in range(1, len(data)):
        x0, y0 = data[index - 1]
        x1, y1 = data[index]
        if y0 > 0 and y1 <= 0 and x1 != x0 and y1 != y0:
            x_zero = x0 + (-y0) * (x1 - x0) / (y1 - y0)
            return np.vstack([data[:index], [x_zero, 0.0]])

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


def extract_metrics(data: np.ndarray):
    if data.size == 0:
        return None

    data = data[data[:, 0] >= 0]
    data = data[np.argsort(data[:, 0])]
    data = data[data[:, 1] >= 0]
    if len(data) == 0:
        return None

    peak_index = int(np.argmax(data[:, 1]))
    fmax = float(data[peak_index, 1])
    work_peak = float(np.trapezoid(data[: peak_index + 1, 1], data[: peak_index + 1, 0]))
    return fmax, work_peak


def main():
    base_dir = Path(__file__).parent
    data_file = base_dir / "Baseline.txt"
    datasets = parse_baseline_sections(data_file)

    labels = [
        ("#a", "A"),
        ("#a-precut", "A-precut"),
        ("#a-precut-2", "A-precut-2"),
        ("#b", "B"),
        ("#bprecut", "B-precut"),
        ("#b-precut-2", "B-precut-2"),
    ]

    summary = []
    for key, display in labels:
        if key not in datasets:
            continue
        metrics = extract_metrics(datasets[key])
        if metrics is None:
            continue
        summary.append((display, datasets[key], metrics[0], metrics[1]))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

                                                         
    ax = axes[0, 0]
    for key, display in labels[:3]:
        if key in datasets:
            data = force_curve_to_zero(datasets[key])
            if data.size:
                ax.plot(data[:, 0], data[:, 1], linewidth=2, label=display)
    ax.set_xlabel("Displacement")
    ax.set_ylabel("Force")
    ax.set_title("Baseline A")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for key, display in labels[3:]:
        if key in datasets:
            data = force_curve_to_zero(datasets[key])
            if data.size:
                ax.plot(data[:, 0], data[:, 1], linewidth=2, label=display)
    ax.set_xlabel("Displacement")
    ax.set_ylabel("Force")
    ax.set_title("Baseline B")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

                                                   
    names = [item[0] for item in summary]
    fmax_values = [item[2] for item in summary]
    work_values = [item[3] for item in summary]
    x = np.arange(len(names))
    width = 0.35

    ax = axes[1, 0]
    ax.bar(x, fmax_values, color="#4c78a8")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_xlabel("Case")
    ax.set_ylabel(r"$F_{\max}$")
    ax.set_title("Peak force")
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(bottom=0)
    for xpos, value in zip(x, fmax_values):
        ax.text(xpos, value, f"{value:.3f}", ha="center", va="bottom", fontsize=7)

    ax = axes[1, 1]
    ax.bar(x, work_values, color="#f58518")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_xlabel("Case")
    ax.set_ylabel(r"$W_{\mathrm{peak}}$")
    ax.set_title("Work to peak")
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(bottom=0)
    for xpos, value in zip(x, work_values):
        ax.text(xpos, value, f"{value:.3f}", ha="center", va="bottom", fontsize=7)

    fig.suptitle("Baseline Results Summary", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_file = base_dir / "baseline_summary_metrics.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_file.name}")


if __name__ == "__main__":
    main()
