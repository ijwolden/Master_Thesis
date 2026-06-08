import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core-properties'))

def relative_density(t, d, x, z, theta_deg):
    theta = np.radians(theta_deg)
    cot_theta = 1/np.tan(theta)
    cot_alpha = (x/z) - cot_theta
    alpha = np.arctan(1/cot_alpha)
    
    A_h = 2 * t * d                                      
    A_i = t * d                                         
    
                        
    total_volume = A_h * x + A_i * (z / np.sin(theta)) + A_i * (z / np.sin(alpha))
    rho = total_volume / (x * z)
    
    return rho, np.degrees(alpha)

def calculate_strut_stresses(E_s, t, d, x, z, theta_deg):
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
    
                                                               
    E_11_eff = (1/V) * E_s * (A_h * L_h * 1.0 + 
                              A_i * L_L * np.cos(theta)**4 + 
                              A_i * L_R * np.cos(np.pi - alpha)**4)
    
                             
    sigma_11 = 1.0       
    eps_11 = sigma_11 / E_11_eff
    
                                
    nu_12 = 0.22
    eps_22 = -nu_12 * eps_11
    
                                             
    struts = [
        {'name': 'horizontal', 'L': L_h, 'phi': 0.0, 'A': A_h},
        {'name': 'left inclined', 'L': L_L, 'phi': theta, 'A': A_i},
        {'name': 'right inclined', 'L': L_R, 'phi': np.pi - alpha, 'A': A_i}
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
    
    return strut_stresses, L_h, L_L, L_R

def estimate_D_from_stress(strut_stresses, theta_deg):
    sigma_h, sigma_L, sigma_R = strut_stresses
    
                                                                             
    D_equilateral = 1.9
    
                                             
    sigma_mean = np.mean([sigma_h, sigma_L, sigma_R])
    sigma_std = np.std([sigma_h, sigma_L, sigma_R])
    uniformity = max(0, 1 - (sigma_std / sigma_mean))                  
    
                                           
    max_stress = max(strut_stresses)
    concentration_factor = max_stress / sigma_mean
    
                                                     
    D_uniformity = D_equilateral * max(uniformity, 0.1)**0.5                                   
    D_concentration = D_equilateral / concentration_factor**0.3                
    
                                       
    angle_factor = 1.0 - abs(theta_deg - 60) / 90                
    D_geometric = 1.0 + 0.9 * angle_factor
    
                                 
    D_estimate = 0.4 * D_uniformity + 0.4 * D_concentration + 0.2 * D_geometric
    
    return {
        'D_estimate': D_estimate,
        'D_uniformity': D_uniformity,
        'D_concentration': D_concentration,
        'D_geometric': D_geometric,
        'uniformity': uniformity,
        'concentration_factor': concentration_factor,
        'max_stress': max_stress,
        'mean_stress': sigma_mean
    }

def characteristic_length(x, z, theta_deg):
    theta = np.radians(theta_deg)
    cot_theta = 1/np.tan(theta)
    cot_alpha = (x/z) - cot_theta
    alpha = np.arctan(1/cot_alpha)
    
                   
    L_h = x
    L_L = z / np.sin(theta)
    L_R = z / np.sin(alpha)
    
                                                  
    l_char = L_h                       
    
    return l_char, L_h, L_L, L_R

def effective_fracture_toughness(D, rho_bar, sigma_f, l):
    K_IC = D * rho_bar * sigma_f * np.sqrt(l)
    return K_IC

def calculate_fracture_properties(t, d, x, z, theta_deg, sigma_f, D, E_s=200000):
                          
    rho, alpha_deg = relative_density(t, d, x, z, theta_deg)
    
                        
    strut_stresses, L_h, L_L, L_R = calculate_strut_stresses(E_s, t, d, x, z, theta_deg)
    
                                         
    D_analysis = estimate_D_from_stress(strut_stresses, theta_deg)
    
                                                                        
    l_char, _, _, _ = characteristic_length(x, z, theta_deg)
    
                                  
    K_IC = effective_fracture_toughness(D, rho, sigma_f, l_char)
    
                                     
    K_IC_estimated = effective_fracture_toughness(D_analysis['D_estimate'], rho, sigma_f, l_char)
    
                                         
    K_IC_MPa_sqrtm = K_IC / np.sqrt(1000)
    K_IC_estimated_MPa_sqrtm = K_IC_estimated / np.sqrt(1000)
    
    return {
        'rho': rho,
        'alpha_deg': alpha_deg,
        'l_char': l_char,
        'L_h': L_h,
        'L_L': L_L,
        'L_R': L_R,
        'K_IC': K_IC,
        'K_IC_MPa_sqrtm': K_IC_MPa_sqrtm,
        'K_IC_estimated': K_IC_estimated,
        'K_IC_estimated_MPa_sqrtm': K_IC_estimated_MPa_sqrtm,
        'D_estimated': D_analysis['D_estimate'],
        'D_uniformity': D_analysis['D_uniformity'],
        'D_concentration': D_analysis['D_concentration'],
        'uniformity': D_analysis['uniformity'],
        'concentration_factor': D_analysis['concentration_factor'],
        'strut_stresses': strut_stresses
    }

if __name__ == "__main__":
                
    t = 0.1       
    d = 0.1       
    x = 1.0       
    z = 1.0       
    sigma_f = 250                                                      
    
                                 
    D_conservative = 0.5                          
    D_moderate = 1.0                          
    D_literature = 1.9                                                      
    
    print("=" * 80)
    print("EFFECTIVE FRACTURE TOUGHNESS ANALYSIS")
    print("=" * 80)
    print(f"\nGeometry: b={x} mm, h={z} mm, t×d={t}×{d} mm²")
    print(f"Material: σ_f={sigma_f} MPa (fracture strength)\n")
    
    print("Topology factors:")
    print(f"  D = {D_conservative:.1f} (conservative estimate for triangular lattices)")
    print(f"  D = {D_moderate:.1f} (moderate estimate)")
    print(f"  D = {D_literature:.1f} (equilateral triangle, Romijn & Fleck 2007)")
    print("\n" + "=" * 80)
    
                 
    test_angles = [70, 75, 80, 85]
    
    print("\n" + "=" * 80)
    print("STRESS-BASED D ESTIMATION")
    print("=" * 80)
    print(f"\n{'θ (°)':>6} {'α (°)':>7} {'D est.':>8} {'Uniform':>8} {'Concen.':>8} {'Stress Ratio':>14}")
    print(f"{'':>6} {'':>7} {'':>8} {'method':>8} {'method':>8} {'max/mean':>14}")
    print("-" * 80)
    
    for theta_deg in test_angles:
        props = calculate_fracture_properties(t, d, x, z, theta_deg, sigma_f, D_moderate)
        print(f"{theta_deg:6.0f} {props['alpha_deg']:7.2f} {props['D_estimated']:8.3f} {props['D_uniformity']:8.3f} "
              f"{props['D_concentration']:8.3f} {props['concentration_factor']:14.2f}")
    
    print("\n" + "=" * 80)
    print("FRACTURE TOUGHNESS RESULTS")
    print("=" * 80)
    print(f"\n{'θ (°)':>6} {'ρ̄ (%)':>8} {'ℓ (mm)':>9} {'D (est)':>10} {'K_IC* (est)':>14} {'K_IC* (D=1.9)':>15}")
    print(f"{'':>6} {'':>8} {'':>9} {'':>10} {'MPa√m':>14} {'MPa√m':>15}")
    print("-" * 80)
    
    results = {'theta': [], 'rho': [], 'l_char': [], 
               'K_IC_cons': [], 'K_IC_mod': [], 'K_IC_lit': [],
               'K_IC_est': [], 'D_est': [], 'concentration': []}
    
    for theta_deg in test_angles:
        props_cons = calculate_fracture_properties(t, d, x, z, theta_deg, sigma_f, D_conservative)
        props_lit = calculate_fracture_properties(t, d, x, z, theta_deg, sigma_f, D_literature)
        
        results['theta'].append(theta_deg)
        results['rho'].append(props_cons['rho'] * 100)
        results['l_char'].append(props_cons['l_char'])
        results['K_IC_cons'].append(props_cons['K_IC_MPa_sqrtm'])
        results['K_IC_mod'].append(props_cons['K_IC_estimated_MPa_sqrtm'])                 
        results['K_IC_lit'].append(props_lit['K_IC_MPa_sqrtm'])
        results['K_IC_est'].append(props_cons['K_IC_estimated_MPa_sqrtm'])
        results['D_est'].append(props_cons['D_estimated'])
        results['concentration'].append(props_cons['concentration_factor'])
        
        print(f"{theta_deg:6.0f} {props_cons['rho']*100:8.3f} {props_cons['l_char']:9.4f} {props_cons['D_estimated']:10.3f} "
              f"{props_cons['K_IC_estimated_MPa_sqrtm']:14.4f} {props_lit['K_IC_MPa_sqrtm']:15.4f}")
    
    print("\n" + "=" * 80)
    print("COMPARISON OF D VALUES:")
    print(f"  Stress-based estimate: D = {np.mean(results['D_est']):.3f} ± {np.std(results['D_est']):.3f}")
    print(f"  Conservative bound:    D = {D_conservative:.1f}")
    print(f"  Literature (θ=60°):    D = {D_literature:.1f}")
    print("\nRecommendation: Use D ≈ {:.2f} for your geometry (based on stress analysis)".format(np.mean(results['D_est'])))
    print("=" * 80)
    
                  
    theta_range = np.linspace(63.43, 90, 50)
    K_IC_cons_range = []
    K_IC_est_range = []
    K_IC_lit_range = []
    D_est_range = []
    rho_range = []
    concentration_range = []
    
    for theta_deg in theta_range:
        props_cons = calculate_fracture_properties(t, d, x, z, theta_deg, sigma_f, D_conservative)
        props_lit = calculate_fracture_properties(t, d, x, z, theta_deg, sigma_f, D_literature)
        
        K_IC_cons_range.append(props_cons['K_IC_MPa_sqrtm'])
        K_IC_est_range.append(props_cons['K_IC_estimated_MPa_sqrtm'])
        K_IC_lit_range.append(props_lit['K_IC_MPa_sqrtm'])
        D_est_range.append(props_cons['D_estimated'])
        rho_range.append(props_cons['rho'] * 100)
        concentration_range.append(props_cons['concentration_factor'])
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
                                                   
    ax = axes[0, 0]
    ax.plot(theta_range, K_IC_cons_range, 'b-', linewidth=2, label=f'D = {D_conservative}')
    ax.plot(theta_range, K_IC_est_range, 'g-', linewidth=2, label=f'D = estimated')
    ax.plot(theta_range, K_IC_lit_range, 'r-', linewidth=2, label=f'D = {D_literature}')
    ax.set_xlabel('θ')
    ax.set_ylabel('K$_{IC}^*$ (MPa√m)')
    ax.set_title('Fracture Toughness vs Angle')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
                                    
    ax = axes[0, 1]
    ax.plot(theta_range, D_est_range, 'g-', linewidth=2, label='Estimated D')
    ax.axhline(y=D_conservative, color='b', linestyle='--', linewidth=1.5, label=f'D = {D_conservative}')
    ax.axhline(y=D_literature, color='r', linestyle='--', linewidth=1.5, label=f'D = {D_literature}')
    ax.set_xlabel('θ')
    ax.set_ylabel('Topology Factor D')
    ax.set_title('Estimated D vs Angle')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
                                                  
    ax = axes[1, 0]
    ax.plot(theta_range, concentration_range, 'purple', linewidth=2)
    ax.set_xlabel('θ')
    ax.set_ylabel('Stress Concentration Factor')
    ax.set_title('Stress Concentration vs Angle')
    ax.grid(True, alpha=0.3)
    
                                               
    ax = axes[1, 1]
    ax.plot(rho_range, K_IC_est_range, 'g-', linewidth=2, label='Estimated D')
    ax.plot(rho_range, K_IC_lit_range, 'r--', linewidth=1.5, label=f'D = {D_literature}')
    ax.set_xlabel('ρ̄ (%)')
    ax.set_ylabel('K$_{IC}^*$ (MPa√m)')
    ax.set_title('Toughness Scaling with Density')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.3f}'))
    
    plt.tight_layout()
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'plots')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, 'fracture_toughness_analysis.png'), dpi=300, bbox_inches='tight')
    print("\nPlot saved: results/plots/fracture_toughness_analysis.png")
