import numpy as np
import matplotlib.pyplot as plt

def compute_analytical_moduli(a_vals, E_s=1.0, t=0.1, x=1.0, z=1.0, d=1.0):
    E11_list = []
    E22_list = []
    G12_list = []
    nu12_list = []
    nu21_list = []
    
    coeff = (E_s * t * d) / x

    for a in a_vals:
                                                    
        if a <= 0.001 or a >= 0.999:                              
            E11_list.append(np.nan)
            E22_list.append(np.nan)
            G12_list.append(np.nan)
            nu12_list.append(np.nan)
            nu21_list.append(np.nan)
            continue
            
        theta = np.arctan2(z, a * x)
        alpha = np.arctan2(z, (1.0 - a) * x)

        c_t = np.cos(theta); s_t = np.sin(theta)
        c_a = np.cos(alpha); s_a = np.sin(alpha)

                                                                
        C11_star = 2.0 + (c_t**4 / s_t) + (c_a**4 / s_a)
        C22_star = s_t**3 + s_a**3
        C12_star = (c_t**2 * s_t) + (c_a**2 * s_a)
        C66_star = C12_star                              
        C16_star = c_t**3 - c_a**3
        C26_star = (c_t * s_t**2) - (c_a * s_a**2)

                                     
        C = coeff * np.array([
            [C11_star, C12_star, C16_star],
            [C12_star, C22_star, C26_star],
            [C16_star, C26_star, C66_star]
        ])

                                           
        try:
            S = np.linalg.inv(C)
                                
            E11 = 1.0 / S[0, 0]
            E22 = 1.0 / S[1, 1]
            G12 = 1.0 / S[2, 2]
            nu12 = -S[0, 1] / S[0, 0]
            nu21 = -S[0, 1] / S[1, 1]
        except np.linalg.LinAlgError:
            E11, E22, G12, nu12, nu21 = np.nan, np.nan, np.nan, np.nan, np.nan

        E11_list.append(E11)
        E22_list.append(E22)
        G12_list.append(G12)
        nu12_list.append(nu12)
        nu21_list.append(nu21)

    return np.array(E11_list), np.array(E22_list), np.array(G12_list), np.array(nu12_list), np.array(nu21_list)

if __name__ == "__main__":
    a_vals = np.linspace(0.1, 0.9, 100)
    E11, E22, G12, nu12, nu21 = compute_analytical_moduli(a_vals)

                            
    print(f"{'a':<6} | {'E11':<7} | {'E22':<7} | {'G12':<7} | {'nu12':<7}")
    print("-" * 45)
    for a in [0.1, 0.3, 0.5, 0.7, 0.9]:
        idx = np.argmin(np.abs(a_vals - a))
        print(f"{a_vals[idx]:.4f} | {E11[idx]:.4f}  | {E22[idx]:.4f}  | {G12[idx]:.4f}  | {nu12[idx]:.4f}")

                 
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].plot(a_vals, E11, label='$E_{11}$ (Horizontal)', lw=2)
    axes[0].plot(a_vals, E22, label='$E_{22}$ (Vertical)', lw=2)
    axes[0].plot(a_vals, G12, label='$G_{12}$ (Shear)', lw=2)
    axes[0].axvline(0.5, color='gray', linestyle='--', alpha=0.5, label='Symmetric Baseline')
    axes[0].set_xlabel('$a$')
    axes[0].set_ylabel('Effective Modulus [normalized]')
    axes[0].set_title('Analytical Stiffness vs. Apex Position', pad=15)
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(a_vals, nu12, label='$\\nu_{12}$', color='purple', lw=2)
    axes[1].plot(a_vals, nu21, label='$\\nu_{21}$', color='magenta', lw=2, linestyle='--')
    axes[1].axvline(0.5, color='gray', linestyle='--', alpha=0.5)
    axes[1].set_xlabel('$a$')
    axes[1].set_ylabel("Poisson's Ratio")
    axes[1].set_title("Poisson Effect vs. Apex Position", pad=15)
    axes[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('analytical_moduli_corrected.png', dpi=300)
    print("\nSaved plot to analytical_moduli_corrected.png")
