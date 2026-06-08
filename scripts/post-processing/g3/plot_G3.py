import matplotlib.pyplot as plt
import numpy as np


def read_g3_data(filename):
    data = {"G3": ([], []), "G3-Precut": ([], [])}
    current_section = None

    with open(filename, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith("#"):
                section = line[1:]
                if section in data:
                    current_section = section
                else:
                    current_section = None
                continue

            if current_section is None:
                continue

            parts = line.split()
            if len(parts) != 2:
                continue

            try:
                displacement = float(parts[0])
                force = float(parts[1])
            except ValueError:
                continue

            data[current_section][0].append(displacement)
            data[current_section][1].append(force)

    g3 = (np.array(data["G3"][0]), np.array(data["G3"][1]))
    g3_precut = (np.array(data["G3-Precut"][0]), np.array(data["G3-Precut"][1]))
    return g3, g3_precut


def main():
    g3, g3_precut = read_g3_data("g3-f-d.txt")

    plt.figure(figsize=(10, 6))
    plt.plot(g3[0], g3[1], label="G3", linewidth=2, color="C0")
    plt.plot(g3_precut[0], g3_precut[1], label="G3-precut", linewidth=2, linestyle="--", color="C1")

    plt.xlabel("Displacement", fontsize=15)
    plt.ylabel("Force", fontsize=15)
    plt.title("G3", fontsize=15)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=15)
    plt.tick_params(labelsize=13)
    plt.tight_layout()

    output_name = "G3_force_displacement.png"
    plt.savefig(output_name, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Plot saved as '{output_name}'")


if __name__ == "__main__":
    main()
