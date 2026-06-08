import re
from pathlib import Path

import matplotlib.pyplot as plt


NUMBER_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")


def read_combined_force_displacement(file_path: Path, label_suffix: str = ""):
    lines = [line.rstrip("\n") for line in file_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError("Input file is empty.")

    header_cells = lines[0].split("\t")
    job_columns = []
    for idx, cell in enumerate(header_cells):
        label = cell.strip()
        match = NUMBER_PATTERN.search(label)
        if match:
            base_label = str(float(match.group(0)))
            job_columns.append((base_label + label_suffix, idx, base_label))

    if not job_columns:
        raise ValueError("Could not find job labels in first line.")

    series = {
        label: {"displacement": [], "force": [], "base_label": base_label}
        for label, _, base_label in job_columns
    }
    expected_cols = len(header_cells)

    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) < expected_cols:
            cells.extend([""] * (expected_cols - len(cells)))

        for label, col_idx, _ in job_columns:
            if col_idx + 1 >= len(cells):
                continue

            disp_raw = cells[col_idx].strip()
            force_raw = cells[col_idx + 1].strip()
            if not disp_raw or not force_raw:
                continue

            try:
                displacement = float(disp_raw)
                force = float(force_raw)
            except ValueError:
                continue

            series[label]["displacement"].append(displacement)
            series[label]["force"].append(force)

    return series


def merge_series(*series_dicts):
    merged = {}
    for series in series_dicts:
        merged.update(series)
    return merged


def _filter_non_negative_points(displacement, force):
    points = []
    paired = list(zip(displacement, force))

    for idx, (x_curr, y_curr) in enumerate(paired):
        if y_curr >= 0.0:
            points.append((x_curr, y_curr))
            continue

        if idx == 0:
            continue

        x_prev, y_prev = paired[idx - 1]
        if y_prev < 0.0:
            continue

                                                                  
        if y_prev == y_curr:
            x_zero = x_curr
        else:
            x_zero = x_prev + (0.0 - y_prev) * (x_curr - x_prev) / (y_curr - y_prev)

        points.append((x_zero, 0.0))
        break

    return points


def _truncate_primary_loading(points, drop_fraction: float = 0.2):
    if len(points) < 3:
        return points

    y_values = [y for _, y in points]
    peak_idx = max(range(len(y_values)), key=y_values.__getitem__)
    if peak_idx >= len(points) - 1:
        return points

    cutoff_idx = len(points)
    for idx in range(peak_idx + 1, len(points)):
        previous_force = points[idx - 1][1]
        current_force = points[idx][1]
        if previous_force <= 0:
            continue

        relative_drop = (previous_force - current_force) / previous_force
        if relative_drop >= drop_fraction:
            cutoff_idx = idx
            break

    return points[:cutoff_idx]


def _build_color_map(series):
    base_labels = sorted({data["base_label"] for data in series.values()}, key=float)
    cmap = plt.get_cmap("tab10")
    return {base_label: cmap(idx % 10) for idx, base_label in enumerate(base_labels)}


def plot_force_displacement(series, output_file: Path, title: str, primary_loading_only: bool = False):
    fig, ax = plt.subplots(figsize=(12, 7))
    color_map = _build_color_map(series)

    labels_sorted = sorted(
        series.keys(),
        key=lambda lbl: (float(series[lbl]["base_label"]), 1 if lbl.endswith(" precut") else 0),
    )

    for label in labels_sorted:
        x = series[label]["displacement"]
        y = series[label]["force"]
        filtered_points = _filter_non_negative_points(x, y)
        if primary_loading_only:
            filtered_points = _truncate_primary_loading(filtered_points)
        if not filtered_points:
            continue
        x_filtered, y_filtered = zip(*filtered_points)
        base_label = series[label]["base_label"]
        is_precut = label.endswith(" precut")
        ax.plot(
            x_filtered,
            y_filtered,
            linewidth=2,
            linestyle=":" if is_precut else "-",
            color=color_map[base_label],
            label=label,
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Displacement")
    ax.set_ylabel("Force")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _build_shared_color_map(series_intact: dict) -> dict:
    base_labels = sorted(series_intact.keys(), key=float)
    cmap = plt.get_cmap("tab10")
    return {lbl: cmap(i % 10) for i, lbl in enumerate(base_labels)}


def plot_single_panel(
    source: dict,
    label_suffix: str,
    panel_title: str,
    output_file: Path,
    color_map: dict,
) -> None:
    base_labels = sorted(
        {k.removesuffix(label_suffix) for k in source}, key=float
    )

    fig, ax = plt.subplots(figsize=(7, 5.5), constrained_layout=True)

    for lbl in base_labels:
        key = lbl + label_suffix
        if key not in source:
            continue
        x = source[key]["displacement"]
        y = source[key]["force"]
        pts = _filter_non_negative_points(x, y)
        if not pts:
            continue
        xf, yf = zip(*pts)
        ax.plot(
            xf, yf,
            linewidth=1.8,
            color=color_map[lbl],
            label=f"$a = {lbl}$",
        )

    ax.set_title(panel_title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Displacement", fontsize=11)
    ax.set_ylabel("Force", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=9, frameon=False)

    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    script_dir = Path(__file__).resolve().parent
    input_file = script_dir / "f-d-all.txt"
    input_file_precut = script_dir / "f-d-precut.txt"
    output_file_full = script_dir / "f-d-all.png"
    output_file_primary = script_dir / "f-d-all-primary-loading.png"
    output_file_intact = script_dir / "f-d-intact.png"
    output_file_precut_out = script_dir / "f-d-precut-plot.png"

    series_regular = read_combined_force_displacement(input_file)
    series_precut = read_combined_force_displacement(input_file_precut, label_suffix=" precut")
    series = merge_series(series_regular, series_precut)
    plot_force_displacement(
        series,
        output_file_full,
        title="Force-Displacement Curves for symmetric designs",
        primary_loading_only=False,
    )
    plot_force_displacement(
        series,
        output_file_primary,
        title="Primary Loading Branch",
        primary_loading_only=True,
    )
    color_map = _build_shared_color_map(series_regular)
    plot_single_panel(series_regular, "", "Intact", output_file_intact, color_map)
    plot_single_panel(series_precut, " precut", "Precut", output_file_precut_out, color_map)
    print(f"Saved plot: {output_file_full}")
    print(f"Saved plot: {output_file_primary}")
    print(f"Saved plot: {output_file_intact}")
    print(f"Saved plot: {output_file_precut_out}")


if __name__ == "__main__":
    main()
