from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "cut.txt"
OUT_FILE = BASE_DIR / "cut_force_displacement.png"


def parse_sections(file_path: Path) -> dict[str, np.ndarray]:
    datasets: dict[str, np.ndarray] = {}
    current_label = None
    current_data: list[tuple[float, float]] = []

    for raw_line in file_path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            if current_label and current_data:
                datasets[current_label] = np.array(current_data, dtype=float)
            current_label = line[1:].strip()
            current_data = []
            continue

        parts = line.split()
        if len(parts) >= 2:
            current_data.append((float(parts[0]), float(parts[1])))

    if current_label and current_data:
        datasets[current_label] = np.array(current_data, dtype=float)

    return datasets


def main() -> None:
    datasets = parse_sections(DATA_FILE)

    label_map = {
        "loc3": "Right (a=5)",
        "middle": "Middle (a=4.5)",
        "loc2": "Left (a=4)",
    }
    color_map = {
        "loc3": "#1f3a5f",
        "middle": "#b64926",
        "loc2": "#3d6b35",
    }

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharey=True)

    for axis, strut_count in zip(axes, ["1", "2", "3"]):
        for prefix in ["loc3", "middle", "loc2"]:
            key = f"{prefix}-{strut_count}"
            data = datasets[key]
            axis.plot(
                data[:, 0],
                data[:, 1],
                linewidth=2.2,
                label=label_map[prefix],
                color=color_map[prefix],
            )

        axis.set_title(f"Cut Through {strut_count} Strut{'s' if strut_count != '1' else ''}")
        axis.set_xlabel("Displacement [mm]")
        axis.grid(True, alpha=0.25)
        axis.ticklabel_format(axis="both", style="sci", scilimits=(-3, 3))

    axes[0].set_ylabel("Force [kN]")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.05),
    )
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(OUT_FILE, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
