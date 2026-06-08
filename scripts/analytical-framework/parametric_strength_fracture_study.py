import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core-properties'))
from compliance_matrix_FINAL import compliance_matrix_final

def calculate_effective_strength(E_s, t, d, x, z, theta_deg, sigma_ys):
    theta = np.radians(theta_deg)
    cot_theta = 1/np.tan(theta)
    cot_alpha = (x/z) - cot_theta
    alpha = np.arctan(1/cot_alpha)
    
    A_h = 2 * t * d                                      
    A_i = t * d                                         
    V = x * z
    
                   
    L_h = x
    L_L = z / np.sin(theta)
    L_R = z / np.sin(alpha)
    
                                           
    E_11, E_22, nu_12, nu_21, G_12, rho, alpha_deg = compliance_matrix_final(E_s, t, d, x, z, theta_deg)
    
                             
    sigma_11 = 1.0       
    eps_11 = sigma_11 / E_11
    eps_22 = -nu_12 * eps_11
    
                  
    phi_h = 0.0
    phi_L = theta
    phi_R = np.pi/2 + alpha                                     
    
                              
    struts = [
        {'name': 'horizontal', 'L': L_h, 'phi': phi_h, 'A': A_h},
        {'name': 'left inclined', 'L': L_L, 'phi': phi_L, 'A': A_i},
        {'name': 'right inclined', 'L': L_R, 'phi': phi_R, 'A': A_i}
    ]
    
    strut_stresses = []
    for strut in struts:
        L = strut['L']
        phi = strut['phi']
        A = strut['A']
        c = np.cos(phi)
        s = np.sin(phi)
        
        delta_L = L * (eps_11 * c**2 + eps_22 * s**2)
        F = (E_s * A / L) * delta_L
        sigma_strut = abs(F / A)
        strut_stresses.append(sigma_strut)
    
                         
    max_stress = max(strut_stresses)
    critical_idx = np.argmax(strut_stresses)
    critical_strut = struts[critical_idx]['name']
    
                        
    sigma_eff = sigma_11 * (sigma_ys / max_stress)
    sigma_normalized = sigma_eff / sigma_ys
    
    return {
        'sigma_eff': sigma_eff,
        'sigma_norm': sigma_normalized,
        'rho': rho,
        'alpha_deg': alpha_deg,
        'critical_strut': critical_strut,
        'strut_stresses': strut_stresses,
        'max_stress': max_stress,
        'E_11': E_11
    }


def calculate_fracture_toughness(t, d, x, z, theta_deg, sigma_f, E_s):
    theta = np.radians(theta_deg)
    cot_theta = 1/np.tan(theta)
    cot_alpha = (x/z) - cot_theta
    alpha = np.arctan(1/cot_alpha)
    
    A_h = 2 * t * d                                      
    A_i = t * d                                         
    V = x * z
    
                                                 
    E_11, E_22, nu_12, nu_21, G_12, rho, alpha_deg = compliance_matrix_final(E_s, t, d, x, z, theta_deg)
    
                                              
    l_char = x
    
                                                             
    L_h = x
    L_L = z / np.sin(theta)
    L_R = z / np.sin(alpha)
    
    sigma_11 = 1.0
    eps_11 = sigma_11 / E_11
    eps_22 = -nu_12 * eps_11
    
                    
    sigma_h = abs((E_s * A_h / L_h) * L_h * eps_11) / A_h
    sigma_L = abs((E_s * A_i / L_L) * L_L * (eps_11 * np.cos(theta)**2 + eps_22 * np.sin(theta)**2)) / A_i
    sigma_R = abs((E_s * A_i / L_R) * L_R * (eps_11 * np.cos(np.pi/2 + alpha)**2 + eps_22 * np.sin(np.pi/2 + alpha)**2)) / A_i
    
                                       
    strut_stresses = [sigma_h, sigma_L, sigma_R]
    sigma_mean = np.mean(strut_stresses)
    sigma_std = np.std(strut_stresses)
    uniformity = max(0, 1 - (sigma_std / sigma_mean))
    
    max_stress = max(strut_stresses)
    concentration_factor = max_stress / sigma_mean
    
                                                
    D_equilateral = 1.9                                             
    D_uniformity = D_equilateral * max(uniformity, 0.1)**0.5
    D_concentration = D_equilateral / concentration_factor**0.3
    
    angle_factor = 1.0 - abs(theta_deg - 60) / 90
    D_geometric = 1.0 + 0.9 * angle_factor
    
    D_estimate = 0.4 * D_uniformity + 0.4 * D_concentration + 0.2 * D_geometric
    
                     
    K_IC = D_estimate * rho * sigma_f * np.sqrt(l_char)
    
                                                       
    K_IC_conservative = 0.5 * rho * sigma_f * np.sqrt(l_char)
    K_IC_literature = 1.9 * rho * sigma_f * np.sqrt(l_char)
    
    return {
        'K_IC': K_IC,
        'K_IC_conservative': K_IC_conservative,
        'K_IC_literature': K_IC_literature,
        'D_estimate': D_estimate,
        'rho': rho,
        'alpha_deg': alpha_deg,
        'uniformity': uniformity,
        'concentration_factor': concentration_factor
    }


def parametric_study(t=0.1, d=0.1, x=1.0, z=1.0, E_s=200000, sigma_ys=250, sigma_f=250, n_points=100):
    
                       
    theta_min = np.degrees(np.arctan(2))                        
    theta_max = 90.0                   
    
    theta_range = np.linspace(theta_min, theta_max - 0.1, n_points)                     
    
             
    results = {
        'theta': [],
        'alpha': [],
        'rho': [],
        'E_11': [],
        'sigma_eff': [],
        'sigma_norm': [],
        'K_IC': [],
        'K_IC_conservative': [],
        'K_IC_literature': [],
        'D_estimate': [],
        'critical_strut': [],
        'max_strut_stress': [],
        'uniformity': [],
        'concentration': []
    }
    
    print("=" * 80)
    print("PARAMETRIC STUDY: STRENGTH & FRACTURE TOUGHNESS")
    print("=" * 80)
    print(f"\nGeometry: b={x} mm, h={z} mm, t×d={t}×{d} mm²")
    print(f"Material: E_s={E_s} MPa, σ_ys={sigma_ys} MPa, σ_f={sigma_f} MPa")
    print(f"Angle range: θ = [{theta_min:.2f}°, {theta_max:.2f}°]")
    print(f"Number of points: {n_points}\n")
    
    print("Computing properties for all angles...")
    
    for theta_deg in theta_range:
                  
        strength_data = calculate_effective_strength(E_s, t, d, x, z, theta_deg, sigma_ys)
        
                            
        fracture_data = calculate_fracture_toughness(t, d, x, z, theta_deg, sigma_f, E_s)
        
                       
        results['theta'].append(theta_deg)
        results['alpha'].append(strength_data['alpha_deg'])
        results['rho'].append(strength_data['rho'])
        results['E_11'].append(strength_data['E_11'])
        results['sigma_eff'].append(strength_data['sigma_eff'])
        results['sigma_norm'].append(strength_data['sigma_norm'])
        results['K_IC'].append(fracture_data['K_IC'])
        results['K_IC_conservative'].append(fracture_data['K_IC_conservative'])
        results['K_IC_literature'].append(fracture_data['K_IC_literature'])
        results['D_estimate'].append(fracture_data['D_estimate'])
        results['critical_strut'].append(strength_data['critical_strut'])
        results['max_strut_stress'].append(strength_data['max_stress'])
        results['uniformity'].append(fracture_data['uniformity'])
        results['concentration'].append(fracture_data['concentration_factor'])
    
                             
    for key in results:
        if key != 'critical_strut':
            results[key] = np.array(results[key])
    
    print("✓ Complete!\n")
    
                              
    print("=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print(f"Relative density: ρ̄ ∈ [{results['rho'].min():.4f}, {results['rho'].max():.4f}]")
    print(f"Effective stiffness: E* ∈ [{results['E_11'].min():.0f}, {results['E_11'].max():.0f}] MPa")
    print(f"Effective strength: σ* ∈ [{results['sigma_eff'].min():.2f}, {results['sigma_eff'].max():.2f}] MPa")
    print(f"Normalized strength: σ*/σ_ys ∈ [{results['sigma_norm'].min():.4f}, {results['sigma_norm'].max():.4f}]")
    print(f"Fracture toughness: K_IC* ∈ [{results['K_IC'].min():.3f}, {results['K_IC'].max():.3f}] MPa√m")
    print(f"Topology factor: D ∈ [{results['D_estimate'].min():.3f}, {results['D_estimate'].max():.3f}]")
    
                                   
    critical_counts = {}
    for strut in results['critical_strut']:
        critical_counts[strut] = critical_counts.get(strut, 0) + 1
    
    print(f"\nCritical strut statistics:")
    for strut, count in critical_counts.items():
        print(f"  {strut}: {count}/{n_points} ({100*count/n_points:.1f}%)")
    
    return results


def create_plots(results):
    
    fig = plt.figure(figsize=(18, 12))
    
                                        
    rho_percent = results['rho'] * 100
    
                                    
    ax1 = plt.subplot(3, 3, 1)
    ax1.plot(results['theta'], results['sigma_eff'], 'b-', linewidth=2)
    ax1.set_xlabel('θ', fontsize=11)
    ax1.set_ylabel('σ* (MPa)', fontsize=11)
    ax1.set_title('Effective Strength vs Angle', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
                                     
    ax2 = plt.subplot(3, 3, 2)
    ax2.plot(results['theta'], results['sigma_norm'], 'g-', linewidth=2)
    ax2.set_xlabel('θ', fontsize=11)
    ax2.set_ylabel('σ*/σ_ys', fontsize=11)
    ax2.set_title('Normalized Strength vs Angle', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
                            
    ax3 = plt.subplot(3, 3, 3)
    ax3.plot(rho_percent, results['sigma_eff'], 'r-', linewidth=2)
    ax3.set_xlabel('ρ̄ (%)', fontsize=11)
    ax3.set_ylabel('σ* (MPa)', fontsize=11)
    ax3.set_title('Strength Scaling with Density', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
                                    
    ax4 = plt.subplot(3, 3, 4)
    ax4.plot(results['theta'], results['K_IC'], 'purple', linewidth=2.5, label='Estimated D')
    ax4.plot(results['theta'], results['K_IC_conservative'], 'orange', linewidth=1.5, 
             linestyle='--', alpha=0.7, label='Conservative (D=0.5)')
    ax4.plot(results['theta'], results['K_IC_literature'], 'cyan', linewidth=1.5, 
             linestyle='--', alpha=0.7, label='Literature (D=1.9)')
    ax4.set_xlabel('θ', fontsize=11)
    ax4.set_ylabel('K_IC* (MPa√m)', fontsize=11)
    ax4.set_title('Fracture Toughness vs Angle', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)
    
                        
    ax5 = plt.subplot(3, 3, 5)
    ax5.plot(rho_percent, results['K_IC'], 'purple', linewidth=2)
    ax5.set_xlabel('ρ̄ (%)', fontsize=11)
    ax5.set_ylabel('K_IC* (MPa√m)', fontsize=11)
    ax5.set_title('K_IC* Scaling with Density', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    
                                   
    ax6 = plt.subplot(3, 3, 6)
    ax6.plot(results['theta'], results['D_estimate'], 'm-', linewidth=2)
    ax6.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='Conservative')
    ax6.axhline(y=1.9, color='cyan', linestyle='--', alpha=0.5, label='Equilateral')
    ax6.set_xlabel('θ', fontsize=11)
    ax6.set_ylabel('Topology Factor D', fontsize=11)
    ax6.set_title('D Factor vs Angle', fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    
                           
    ax7 = plt.subplot(3, 3, 7)
    ax7.plot(results['theta'], results['E_11'], 'darkblue', linewidth=2)
    ax7.set_xlabel('θ', fontsize=11)
    ax7.set_ylabel('E* (MPa)', fontsize=11)
    ax7.set_title('(g) Effective Stiffness vs Angle', fontsize=12, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    
                                  
    ax8 = plt.subplot(3, 3, 8)
    ax8.plot(results['theta'], results['max_strut_stress'], 'darkred', linewidth=2)
    ax8.set_xlabel('θ', fontsize=11)
    ax8.set_ylabel('Max Strut Stress (MPa)', fontsize=11)
    ax8.set_title('(h) Strut Stress for σ_11 = 1 MPa', fontsize=12, fontweight='bold')
    ax8.grid(True, alpha=0.3)
    
                                    
    ax9 = plt.subplot(3, 3, 9)
    ax9.plot(results['theta'], results['concentration'], 'darkorange', linewidth=2)
    ax9.set_xlabel('θ', fontsize=11)
    ax9.set_ylabel('Stress Concentration Factor', fontsize=11)
    ax9.set_title('(i) σ_max / σ_mean vs Angle', fontsize=12, fontweight='bold')
    ax9.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'plots')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'parametric_strength_fracture_study.png'), dpi=300, bbox_inches='tight')
    print("\n✓ Saved: results/plots/parametric_strength_fracture_study.png")
    
                                                 
    fig2, axes = plt.subplots(2, 2, figsize=(12, 10))
    
                              
    ax = axes[0, 0]
    ax.scatter(rho_percent, results['sigma_eff'], c=results['theta'], cmap='viridis', s=50, alpha=0.7)
    ax.set_xlabel('ρ̄ (%)', fontsize=11)
    ax.set_ylabel('σ* (MPa)', fontsize=11)
    ax.set_title('Strength-Density Scaling (color = θ)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('θ', fontsize=10)
    
                          
    ax = axes[0, 1]
    ax.scatter(rho_percent, results['K_IC'], c=results['theta'], cmap='plasma', s=50, alpha=0.7)
    ax.set_xlabel('ρ̄ (%)', fontsize=11)
    ax.set_ylabel('K_IC* (MPa√m)', fontsize=11)
    ax.set_title('K_IC*-Density Scaling (color = θ)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('θ', fontsize=10)
    
                           
    ax = axes[1, 0]
    ax.scatter(results['E_11'], results['sigma_eff'], c=results['theta'], cmap='coolwarm', s=50, alpha=0.7)
    ax.set_xlabel('E* (MPa)', fontsize=11)
    ax.set_ylabel('σ* (MPa)', fontsize=11)
    ax.set_title('Strength vs Stiffness Trade-off', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('θ', fontsize=10)
    
                      
    ax = axes[1, 1]
    ax.scatter(results['sigma_eff'], results['K_IC'], c=results['theta'], cmap='spring', s=50, alpha=0.7)
    ax.set_xlabel('σ* (MPa)', fontsize=11)
    ax.set_ylabel('K_IC* (MPa√m)', fontsize=11)
    ax.set_title('Toughness vs Strength Trade-off', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('θ', fontsize=10)
    
    plt.tight_layout()
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'plots')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'scaling_relationships.png'), dpi=300, bbox_inches='tight')
    print("✓ Saved: results/plots/scaling_relationships.png")
    
    plt.show()


if __name__ == "__main__":
                          
    results = parametric_study(
        t=0.1,
        d=0.1,
        x=1.0,
        z=1.0,
        E_s=200000,
        sigma_ys=250,
        sigma_f=250,
        n_points=100
    )
    
                  
    create_plots(results)
    
    print("\n" + "=" * 80)
    print("PARAMETRIC STUDY COMPLETE")
    print("=" * 80)
