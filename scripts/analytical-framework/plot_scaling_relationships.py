import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import sys

                                        
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core-properties'))
from compliance_matrix_FINAL import compliance_matrix_final

                                                                             
E_s      = 200_000                                
sigma_ys = 250                                  
sigma_f  = 250                                       
d        = 0.1                                              
x        = 1.0                          
z        = 1.0                           

                                                                       
theta_values = [63.43, 70.0, 75.0, 80.0, 85.0]

                                                                     
t_values = np.linspace(0.01, 0.22, 80)


def compute_K_IC(t, d, x, z, theta_deg, E_s, sigma_f):
    theta = np.radians(theta_deg)
    cot_theta = 1 / np.tan(theta)
    cot_alpha = (x / z) - cot_theta
    alpha = np.arctan(1 / cot_alpha)

    A_h = 2 * t * d
    A_i = t * d
    L_h = x
    L_L = z / np.sin(theta)
    L_R = z / np.sin(alpha)

    E_11, E_22, nu_12, nu_21, G_12, rho_raw, alpha_deg =\
        compliance_matrix_final(E_s, t, d, x, z, theta_deg)

    rho = rho_raw / d                                                     

                                                         
    eps_11 = 1.0 / E_11
    eps_22 = -nu_12 * eps_11
    sigma_h = abs(E_s * A_h / L_h * L_h * eps_11) / A_h
    sigma_L = abs(E_s * A_i / L_L * L_L * (
        eps_11 * np.cos(theta)**2 + eps_22 * np.sin(theta)**2)) / A_i
    sigma_R = abs(E_s * A_i / L_R * L_R * (
        eps_11 * np.cos(np.pi / 2 + alpha)**2 +
        eps_22 * np.sin(np.pi / 2 + alpha)**2)) / A_i

    stresses = [sigma_h, sigma_L, sigma_R]
    sigma_mean = np.mean(stresses)
    uniformity = max(0.0, 1 - np.std(stresses) / sigma_mean)
    conc = max(stresses) / sigma_mean
    angle_factor = 1.0 - abs(theta_deg - 60) / 90
    D = 0.4 * 1.9 * max(uniformity, 0.1)**0.5 +\
        0.4 * 1.9 / conc**0.3 +\
        0.2 * (1.0 + 0.9 * angle_factor)

    K_IC = D * rho * sigma_f * np.sqrt(x * 1e-3)                                  
    return K_IC, rho


all_rho = []
all_E   = []
all_sig = []
all_K   = []

for theta_deg in theta_values:
    for t in t_values:
        E_11, E_22, nu_12, nu_21, G_12, rho_raw, _ =\
            compliance_matrix_final(E_s, t, d, x, z, theta_deg)
        rho = rho_raw / d                        

                                                          
        theta = np.radians(theta_deg)
        cot_theta = 1 / np.tan(theta)
        cot_alpha = (x / z) - cot_theta
        alpha = np.arctan(1 / cot_alpha)
        A_i = t * d
        L_L = z / np.sin(theta)
        eps_11 = 1.0 / E_11
        eps_22 = -nu_12 * eps_11
        sigma_struts = [
            abs(E_s * (2 * t * d) / x * x * eps_11) / (2 * t * d),
            abs(E_s * A_i / L_L * L_L * (
                eps_11 * np.cos(theta)**2 + eps_22 * np.sin(theta)**2)) / A_i,
        ]
        sigma_eff = 1.0 * (sigma_ys / max(sigma_struts))

        K_IC, _ = compute_K_IC(t, d, x, z, theta_deg, E_s, sigma_f)

        all_rho.append(rho)
        all_E.append(E_11 / E_s)
        all_sig.append(sigma_eff / sigma_ys)
        all_K.append(K_IC)

rho    = np.array(all_rho)
E_norm = np.array(all_E)
s_norm = np.array(all_sig)
K_arr  = np.array(all_K)

                                                                             
def power_law_fit(rho_arr, y_arr):
    log_rho = np.log(rho_arr)
    log_y   = np.log(y_arr)
    n, log_C = np.polyfit(log_rho, log_y, 1)
    r2 = float(np.corrcoef(log_rho, log_y)[0, 1] ** 2)
    return n, np.exp(log_C), r2

n_E,   C_E,   r2_E   = power_law_fit(rho, E_norm)
n_sig, C_sig, r2_sig = power_law_fit(rho, s_norm)
n_K,   C_K,   r2_K   = power_law_fit(rho, K_arr)

print(f"\nPower-law fits (t-variation, multiple theta):")
print(f"  E*/E_s    = {C_E:.4g}  x rho^{n_E:.2f}   R^2 = {r2_E:.3f}")
print(f"  s*/sy     = {C_sig:.4g} x rho^{n_sig:.2f}   R^2 = {r2_sig:.3f}")
print(f"  K_Ic*     = {C_K:.4g} x rho^{n_K:.2f}   R^2 = {r2_K:.3f}")

rho_fit  = np.linspace(rho.min() * 0.9, rho.max() * 1.1, 300)
E_fit    = C_E   * rho_fit ** n_E
sig_fit  = C_sig * rho_fit ** n_sig
K_fit    = C_K   * rho_fit ** n_K

                                                                             
fig, axes = plt.subplots(1, 3, figsize=(12, 4.0))

scatter_kw = dict(s=12, alpha=0.55, color='steelblue', zorder=3, label='_nolegend_')
line_kw    = dict(color='crimson', linewidth=2.0, zorder=4)

panels = [
    (E_norm,  E_fit,  r2_E,  n_E,
     r'$E_{11}^* / E_s$'),
    (s_norm, sig_fit, r2_sig, n_sig,
     r'$\sigma^* / \sigma_y$'),
    (K_arr,   K_fit,  r2_K,  n_K,
     r'$K_{Ic}^*$ (MPa$\sqrt{\mathrm{m}}$)'),
]
panel_labels = ['(a)', '(b)', '(c)']

for ax, plabel, (y_data, y_fit, r2, n, ylabel) in zip(axes, panel_labels, panels):
    ax.scatter(rho, y_data, **scatter_kw)
    fit_label = (r'$\propto \bar{\rho}^{' + f'{n:.2f}' + r'}$'
                 + f'\n$R^2 = {r2:.2f}$')
    ax.plot(rho_fit, y_fit, **line_kw, label=fit_label)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\bar{\rho}$', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(plabel, fontsize=12, loc='left', fontweight='bold')
    ax.legend(fontsize=9.5, framealpha=0.9)
    ax.grid(True, which='both', linestyle='--', linewidth=0.4, alpha=0.5)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

plt.tight_layout(pad=1.5)

                                                                             
script_dir = os.path.dirname(os.path.abspath(__file__))
local_out  = os.path.join(script_dir, '..', 'results', 'plots')
thesis_out = os.path.join(script_dir,
                          '..', '..', '..', '..', 'master-thesis',
                          'Images', 'results')

for out_dir in (local_out, thesis_out):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, 'scaling_relationships.png')
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f"Saved: {os.path.normpath(path)}")

plt.show()
