import sys, os
import numpy as np
import matplotlib.pyplot as plt

                                                                             
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT,
    'Analytical-framework', 'single-cell-analysis', 'core-properties'))
from compliance_matrix_FINAL import compliance_matrix_final

DATA_FILE = os.path.join(
    REPO_ROOT, 'master-thesis', 'Parametric-Abaqus',
    'Symmetric', 'f-d-all.txt')

                                                                             
N_CELLS   = 9                           
b         = 1.0                         
h         = 1.0                          
d_strut   = 0.1                                       
t_incl    = 0.1                                                
W         = N_CELLS * b                                           
H         = N_CELLS * h                        

                                                                             
E_s_ref = 200.0       


def read_symmetric_data(filepath):
    with open(filepath) as f:
        header = f.readline().strip().split()                            
        n = len(header)
        rows = []
        for line in f:
            parts = line.split()
            if len(parts) == 2 * n:
                try:
                    rows.append([float(p) for p in parts])
                except ValueError:
                    pass

    data = np.array(rows)

    results = {}
    for i, name in enumerate(header):
        x_apex = float(name)
        disp  = data[:, i * 2]
        force = data[:, i * 2 + 1]

        idx_max = int(np.argmax(force))
        F_max   = float(force[idx_max])
        d_max   = float(disp[idx_max])

                                                        
        mask = disp <= 0.10 * d_max
        if mask.sum() >= 3:
            k = float(np.polyfit(disp[mask], force[mask], 1)[0])
        else:
            k = np.nan

        results[x_apex] = {'k': k, 'F_max': F_max, 'd_max': d_max}

    return results


def k_to_Estar(k):
    return k * W / (H * d_strut)


def analytical_E11(x_apex, E_s):
    theta_deg = np.degrees(np.arctan(h / (x_apex * b)))
    E11, *_ = compliance_matrix_final(E_s, t_incl, d_strut, b, h, theta_deg)
    return E11


def calibrate_E_s(fe_data, x_ref=0.5):
    k_ref   = fe_data[x_ref]['k']
    E_star_FE = k_to_Estar(k_ref)
                                                          
    E11_at_Es1 = analytical_E11(x_ref, 1.0)                         
    E_s_cal = E_star_FE / E11_at_Es1
    return E_s_cal


def main():
    print("=" * 70)
    print("ANALYTICAL–FE STIFFNESS COMPARISON  (Symmetric Design Group)")
    print("=" * 70)

                                                                             
    fe = read_symmetric_data(DATA_FILE)
    x_values = sorted(fe.keys())

                                                                             
    E_s_cal = calibrate_E_s(fe, x_ref=0.5)
    print(f"\nBack-calculated E_s from FE baseline (x=0.5): {E_s_cal:.1f} GPa")
    print(f"  Expected ≈ 200 GPa (FE solid modulus).")
    print(f"  Relative deviation from 200 GPa: {100*(E_s_cal-200)/200:+.1f}%")
    print(f"  Note: the comparison is interpreted in the same kN-mm-GPa")
    print(f"  unit system used for the thesis post-processing.\n")

                                                                             
    header = (f"{'x':>5}  {'θ (°)':>7}  {'k_FE':>8}  {'E*_FE':>8}  "
              f"{'E11_cal':>9}  {'k_pred':>8}  {'err%':>7}  "
              f"{'E*/Es_FE':>9}  {'E*/Es_an':>9}")
    print(header)
    print("-" * len(header))

    rows = []
    for x in x_values:
        theta_deg = np.degrees(np.arctan(h / (x * b)))
        k_FE      = fe[x]['k']
        E_star_FE = k_to_Estar(k_FE)                    
        E11_cal   = analytical_E11(x, E_s_cal)          
        k_pred    = E11_cal * H * d_strut / W             
        err_pct   = 100 * (k_pred - k_FE) / k_FE

                                                                           
        norm_FE  = E_star_FE / E_s_ref
        norm_an  = analytical_E11(x, E_s_ref) / E_s_ref

        rows.append((x, theta_deg, k_FE, E_star_FE, E11_cal, k_pred,
                     err_pct, norm_FE, norm_an))
        print(f"{x:5.1f}  {theta_deg:7.2f}  {k_FE:8.5f}  "
              f"{E_star_FE:8.3f}  {E11_cal:9.3f}  {k_pred:8.5f}  "
              f"{err_pct:+7.1f}  {norm_FE:9.5f}  {norm_an:9.5f}")

    rows = np.array(rows)

                                                                             
    print("\n--- Directional Trend Comparison (normalised to x = 0.5) ---")
    idx05 = list(x_values).index(0.5)
    ref_k = rows[idx05, 2]
    ref_an = rows[idx05, 5]
    print(f"{'x':>5}  {'k_FE / k(0.5)':>14}  {'k_pred / k(0.5)':>16}")
    print("-" * 40)
    for i, x in enumerate(x_values):
        k_ratio_FE = rows[i, 2] / ref_k
        k_ratio_an = rows[i, 5] / ref_an
        print(f"{x:5.1f}  {k_ratio_FE:14.4f}  {k_ratio_an:16.4f}")

    print(f"\nMean absolute prediction error: "
          f"{np.mean(np.abs(rows[:, 6])):.1f}%")
    print(f"Max absolute prediction error:  "
          f"{np.max(np.abs(rows[:, 6])):.1f}%")

                                                                            
    x_arr         = rows[:, 0]
    E_star_FE_arr = rows[:, 3]                          
    E_star_an_arr = rows[:, 4]                                              

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x_arr, E_star_FE_arr, 'o-', color='steelblue',
            label=r'FE  ($E^* = k_\mathrm{FE}\,W/Hd$)',
            linewidth=1.8, markersize=7)
    ax.plot(x_arr, E_star_an_arr, 's--', color='firebrick',
            label=r'Analytical (homogenisation)',
            linewidth=1.8, markersize=7)

    ax.set_xlabel('$a$', fontsize=12)
    ax.set_ylabel(r'Effective modulus $E^*$ [MPa]', fontsize=12)
    ax.set_title('Effect of $a$ on Effective Modulus\n'
                 'Symmetric design', fontsize=11)
    ax.set_xticks(x_arr)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out_path = os.path.join(SCRIPT_DIR, 'stiffness_comparison.png')
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"\nSaved figure → {out_path}")
    plt.show()

    return rows, E_s_cal


if __name__ == '__main__':
    results, E_s_calibrated = main()
