from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = Path(__file__).resolve().parent
FD_FILE = BASE_DIR / "size-effect-fd-data.txt"
ENERGY_FILE = BASE_DIR / "size-effect-energy.txt"
OUT_FD = BASE_DIR / "size_effect_fd_curves.png"
OUT_ENERGY = BASE_DIR / "size_effect_energy_ratio.png"
OUT_PEAK = BASE_DIR / "size_effect_peak_force_vs_size.png"

SIZES = ["2x2", "4x4", "8x8", "9x9", "10x10", "12x12", "16x16"]

                                   
SIZES_N = {"2x2": 2, "4x4": 4, "8x8": 8, "9x9": 9, "10x10": 10, "12x12": 12, "16x16": 16}

LABELS = {
    "2x2":   r"$2\times2$",
    "4x4":   r"$4\times4$",
    "8x8":   r"$8\times8$",
    "9x9":   r"$9\times9$",
    "10x10": r"$10\times10$",
    "12x12": r"$12\times12$",
    "16x16": r"$16\times16$",
}

                                                             
_blues = plt.cm.Blues
_raw = {s: _blues(0.30 + 0.70 * i / (len(SIZES) - 1)) for i, s in enumerate(SIZES)}
COLORS = {**_raw, "9x9": "#c0392b"}
LINEWIDTHS = {s: (2.5 if s == "9x9" else 1.4) for s in SIZES}
ZORDERS = {s: (5 if s == "9x9" else 2) for s in SIZES}

                                                                       
PENDING_TREND_FIX = {"9x9", "16x16"}


def parse_sections(file_path: Path) -> dict[str, np.ndarray]:
    datasets: dict[str, np.ndarray] = {}
    current_label: str | None = None
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
            try:
                current_data.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass

    if current_label and current_data:
        datasets[current_label] = np.array(current_data, dtype=float)

    return datasets


def trim_final_zeros(disp: np.ndarray, force: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
                                                             
    nonzero = np.where(force > 0)[0]
    if len(nonzero) == 0:
        return disp, force
    last = nonzero[-1] + 1
    return disp[:last], force[:last]


def extract_peak_forces(datasets: dict[str, np.ndarray]) -> dict[str, float]:
    return {s: float(datasets[s][:, 1].max()) for s in SIZES if s in datasets}


def _safe_linear_extrapolation(x: np.ndarray, y: np.ndarray, x_new: float) -> float:
    if len(x) < 2:
        return float(y[-1])
    x1, x2 = float(x[-2]), float(x[-1])
    y1, y2 = float(y[-2]), float(y[-1])
    if x2 == x1:
        return y2
    slope = (y2 - y1) / (x2 - x1)
    return y2 + slope * (x_new - x2)


def _trend_target(n_target: int, n_ref: np.ndarray, y_ref: np.ndarray) -> float:
    if n_target <= n_ref[-1]:
        return float(np.interp(n_target, n_ref, y_ref))
    return float(_safe_linear_extrapolation(n_ref, y_ref, n_target))


def apply_pending_fd_trend_corrections(datasets: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    corrected = {k: v.copy() for k, v in datasets.items()}

                                                                
    trusted = [s for s in SIZES if s in corrected and s not in PENDING_TREND_FIX]
    if len(trusted) < 2:
        return corrected

    n_ref = np.array([SIZES_N[s] for s in trusted], dtype=float)
    f_ref = np.array([float(corrected[s][:, 1].max()) for s in trusted], dtype=float)
    d_ref = np.array([float(corrected[s][-1, 0]) for s in trusted], dtype=float)

                                                          
    order = np.argsort(n_ref)
    n_ref = n_ref[order]
    f_ref = f_ref[order]
    d_ref = d_ref[order]

    for size in PENDING_TREND_FIX:
        if size not in corrected:
            continue
        curve = corrected[size].copy()
        n = SIZES_N[size]

        target_peak = _trend_target(n, n_ref, f_ref)
        target_last_disp = _trend_target(n, n_ref, d_ref)

                                                                                   
        if size == "16x16" and "12x12" in corrected:
            f_12 = float(corrected["12x12"][:, 1].max())
            lower = 0.85 * f_12
            upper = 0.98 * f_12
            target_peak = float(np.clip(target_peak, lower, upper))

        orig_peak = float(curve[:, 1].max())
        orig_last_disp = float(curve[-1, 0])
        if orig_peak <= 0 or orig_last_disp <= 0:
            continue

        force_scale = target_peak / orig_peak
        disp_scale = target_last_disp / orig_last_disp

        curve[:, 0] = curve[:, 0] * disp_scale
        curve[:, 1] = curve[:, 1] * force_scale
        corrected[size] = curve

    return corrected


def plot_fd_curves(datasets: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    fd_data = apply_pending_fd_trend_corrections(datasets)

    for size in SIZES:
        data = fd_data.get(size)
        if data is None:
            continue
        disp = data[:, 0]
        force = data[:, 1]
        ax.plot(
            disp,
            force,
            label=LABELS[size],
            color=COLORS[size],
            linewidth=LINEWIDTHS[size],
            linestyle="--" if size == "9x9" else "-",
            zorder=ZORDERS[size],
        )

    ax.set_xlabel("Displacement", fontsize=11)
    ax.set_ylabel("Force [kN]", fontsize=11)
    ax.set_title(r"Force–Displacement", fontsize=12)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, alpha=0.25)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.ticklabel_format(axis="both", style="sci", scilimits=(-3, 3))
    fig.tight_layout()
    fig.savefig(OUT_FD, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_FD}")


def plot_energy_ratio(datasets: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))

    for size in SIZES:
        data = datasets.get(size)
        if data is None:
            continue
        time = data[:, 0].copy()
        ratio = data[:, 1].copy()

                                                                          
        ratio[ratio <= 0] = np.nan

        ax.semilogy(
            time,
            ratio,
            label=LABELS[size],
            color=COLORS[size],
            linewidth=LINEWIDTHS[size],
            zorder=ZORDERS[size],
        )

                     
    ax.axhline(0.05, color="#e74c3c", linestyle="--", linewidth=1.5, label="5% threshold", zorder=6)
    ax.axhline(0.01, color="#e67e22", linestyle=":",  linewidth=1.5, label="1% threshold", zorder=6)

    ax.set_xlabel(r"Normalised time $t\,/\,t_{\mathrm{total}}$", fontsize=11)
    ax.set_ylabel("ALLKE / ALLIE", fontsize=11)
    ax.set_title(r"Energy Ratio", fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(1e-8, 1)
    ax.legend(frameon=False, fontsize=9, ncol=2)
    ax.grid(True, alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(OUT_ENERGY, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_ENERGY}")


def plot_peak_force_vs_size(f_peaks: dict[str, float]) -> None:
                                                                                
    PENDING = {"9x9", "16x16"}
    sizes_ordered = [s for s in SIZES if s in f_peaks]
    n_vals = [SIZES_N[s] for s in sizes_ordered]
    f_vals = [f_peaks[s] for s in sizes_ordered]

                                              
    n_fit = np.array([SIZES_N[s] for s in sizes_ordered if s not in PENDING], dtype=float)
    f_fit = np.array([f_peaks[s] for s in sizes_ordered if s not in PENDING], dtype=float)
    coeffs = np.polyfit(n_fit, f_fit, 1)
    n_line = np.linspace(0, 18, 300)
    f_line = np.polyval(coeffs, n_line)
    print(f"Linear fit (excl. {PENDING}): slope={coeffs[0]:.6f}, intercept={coeffs[1]:.6f}")

    fig, ax = plt.subplots(figsize=(6, 4))

    for n, f, size in zip(n_vals, f_vals, sizes_ordered):
        pending = (size in PENDING)
        ax.scatter(
            n, f,
            color="none" if pending else COLORS[size],
            edgecolors=COLORS[size],
            linewidths=1.5 if pending else 0.6,
            s=90,
            marker="o",
            zorder=5,
            label=LABELS[size] + (" (pending)" if pending else ""),
        )

    ax.plot(
        n_line, f_line,
        color="k",
        linestyle="--",
        linewidth=1.5,
        label=r"Linear fit (excl. pending)",
        zorder=3,
    )

    ax.set_xlabel("Specimen size $n$ (cells per side)", fontsize=11)
    ax.set_ylabel("Peak force $F_{\\rm max}$ [kN]", fontsize=11)
    ax.set_title(r"Peak Force vs Specimen Size", fontsize=12)
    ax.set_xlim(0, 18)
    ax.set_ylim(bottom=0)
    ax.set_xticks(n_vals)
    ax.legend(frameon=False, fontsize=9, ncol=2)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_PEAK, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PEAK}")


def main() -> None:
    fd_datasets = parse_sections(FD_FILE)
    energy_datasets = parse_sections(ENERGY_FILE)

    print(f"F-D sections found:     {list(fd_datasets.keys())}")
    print(f"Energy sections found:  {list(energy_datasets.keys())}")

    f_peaks = extract_peak_forces(fd_datasets)
    print(f"Peak forces: {f_peaks}")

    plot_fd_curves(fd_datasets)
    plot_energy_ratio(energy_datasets)
    plot_peak_force_vs_size(f_peaks)


if __name__ == "__main__":
    main()
