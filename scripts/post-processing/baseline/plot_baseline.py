import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict, Iterable


def parse_datasets(file_path: str) -> Dict[str, np.ndarray]:
    with open(file_path, "r") as f:
        lines = f.readlines()

    datasets = {}
    current_label = None
    current_data = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith("#"):
            if current_label and current_data:
                datasets[current_label] = np.array(current_data)
            current_label = line
            current_data = []
            continue

        try:
            displacement, force = map(float, line.split()[:2])
            current_data.append([displacement, force])
        except ValueError:
            continue

    if current_label and current_data:
        datasets[current_label] = np.array(current_data)

    return datasets


def get_dataset(normalized: Dict[str, np.ndarray], candidates: Iterable[str]) -> np.ndarray | None:
    for key in candidates:
        if key in normalized:
            return normalized[key]
    return None


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


data_file = Path(__file__).with_name("Baseline.txt")
datasets = parse_datasets(str(data_file))

                                                                            
normalized = {k.lower().replace(" ", ""): v for k, v in datasets.items()}
a_data = get_dataset(normalized, ["#a"])
a_precut_data = get_dataset(normalized, ["#aprecut", "#a-precut"])
a_precut2_data = get_dataset(normalized, ["#aprecut2", "#a-precut-2", "#a_precut_2"])
b_data = get_dataset(normalized, ["#b"])
b_precut_data = get_dataset(normalized, ["#bprecut", "#b-precut"])
b_precut2_data = get_dataset(normalized, ["#bprecut2", "#b-precut-2", "#b_precut_2"])

if any(d is None for d in (a_data, a_precut_data, a_precut2_data, b_data, b_precut_data, b_precut2_data)):
    missing = []
    if a_data is None:
        missing.append("A")
    if a_precut_data is None:
        missing.append("A-precut")
    if a_precut2_data is None:
        missing.append("A-precut-2")
    if b_data is None:
        missing.append("B")
    if b_precut_data is None:
        missing.append("B-precut")
    if b_precut2_data is None:
        missing.append("B-precut-2")
    available = ", ".join(sorted(datasets.keys()))
    raise ValueError(
        f"Missing baseline dataset(s): {', '.join(missing)}. Available sections: {available}"
    )

plt.figure(figsize=(10, 6))
plt.plot(*force_curve_to_zero(a_data).T, linewidth=2, color="#1f77b4", label="A")
plt.plot(
    *force_curve_to_zero(a_precut_data).T,
    linewidth=2,
    linestyle="--",
    color="#1f77b4",
    label="A-precut",
)
plt.plot(
    *force_curve_to_zero(a_precut2_data).T,
    linewidth=2,
    linestyle=":",
    color="#1f77b4",
    label="A-precut-2",
)
plt.plot(*force_curve_to_zero(b_data).T, linewidth=2, color="#ff7f0e", label="B")
plt.plot(
    *force_curve_to_zero(b_precut_data).T,
    linewidth=2,
    linestyle="--",
    color="#ff7f0e",
    label="B-precut",
)
plt.plot(
    *force_curve_to_zero(b_precut2_data).T,
    linewidth=2,
    linestyle=":",
    color="#ff7f0e",
    label="B-precut-2",
)

plt.xlabel("Displacement", fontsize=12)
plt.ylabel("Force", fontsize=12)
plt.title("Force Displacement of Baseline Design", fontsize=14, fontweight="bold")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig("baseline_ab_plot.png", dpi=300, bbox_inches="tight")
print("Plot saved as 'baseline_ab_plot.png'")
plt.show()
