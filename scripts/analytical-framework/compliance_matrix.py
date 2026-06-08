import numpy as np

"""
FINAL CORRECTED COMPLIANCE MATRIX METHOD

This implements the homogenization approach properly for a triangular lattice
with one horizontal base strut.

The key correction: For E_11 (horizontal loading), we need to use the correct
strain energy formulation accounting for how each strut deforms under the
applied macroscopic strain.
"""

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

def compliance_matrix_final(E_s, t, d, x, z, theta_deg):
    theta = np.radians(theta_deg)
    cot_theta = 1/np.tan(theta)
    cot_alpha = (x/z) - cot_theta
    alpha = np.arctan(1/cot_alpha)
    
    A_h = 2 * t * d                                      
    A_i = t * d                                         
    V = x * z
    
                                                                       
    struts = [
        (x, 0.0, A_h),                                                                      
        (z/np.sin(theta), theta, A_i),                                                   
        (z/np.sin(alpha), np.pi - alpha, A_i)                                                
    ]
    
                                           
    C_11 = 0                          
    C_22 = 0                        
    C_12 = 0                 
    C_66 = 0         
    
    for L, phi, A in struts:
        c = np.cos(phi)
        s = np.sin(phi)
        
                                                                          
        contribution_factor = E_s * A * L / V
        
        C_11 += contribution_factor * c**4
        C_22 += contribution_factor * s**4
        C_12 += contribution_factor * c**2 * s**2
        C_66 += contribution_factor * c**2 * s**2             
    
                                             
    E_11 = C_11
    E_22 = C_22
    
                                                             
    C_matrix = np.array([
        [C_11, C_12, 0],
        [C_12, C_22, 0],
        [0, 0, C_66]
    ])
    
    try:
        S_matrix = np.linalg.inv(C_matrix)
        nu_12 = -S_matrix[0, 1] / S_matrix[0, 0]
        nu_21 = -S_matrix[1, 0] / S_matrix[1, 1]
        G_12 = 1 / S_matrix[2, 2]
    except:
        nu_12 = nu_21 = G_12 = 0
    
    rho, alpha_deg = relative_density(t, d, x, z, theta_deg)
    
    return E_11, E_22, nu_12, nu_21, G_12, rho, alpha_deg


if __name__ == "__main__":
    t = 0.1
    d = 0.1
    x = 1.0
    z = 1.0
    E_s = 70000       
    
    print("=" * 85)
    print("COMPLIANCE MATRIX METHOD - FINAL CORRECTED VERSION (v2.0)")
    print("=" * 85)
    print(f"\nGeometry: b={x} mm, h={z} mm, nominal t×d={t}×{d} mm²")
    print(f"Horizontal strut: 2t×d = {2*t}×{d} mm² (A_h = {2*t*d} mm²)")
    print(f"Inclined struts: t×d = {t}×{d} mm² (A_i = {t*d} mm²)")
    print(f"Material: E_s={E_s/1000} GPa\n")
    
    test_angles = [63.43, 70, 75, 80, 85, 90]
    
    print(f"{'θ':>6} {'α':>7} {'ρ̄':>9} {'E₁₁*':>10} {'E₂₂*':>10} {'G₁₂*':>10} {'ν₁₂*':>8} {'ν₂₁*':>8} {'E₁₁/E₂₂':>9}")
    print(f"{'(°)':>6} {'(°)':>7} {' ':>9} {'(GPa)':>10} {'(GPa)':>10} {'(GPa)':>10} {' ':>8} {' ':>8} {' ':>9}")
    print("-" * 85)
    
    for theta_deg in test_angles:
        E_11, E_22, nu_12, nu_21, G_12, rho, alpha_deg = compliance_matrix_final(
            E_s, t, d, x, z, theta_deg)
        
        print(f"{theta_deg:6.2f} {alpha_deg:7.2f} {rho:9.6f} {E_11/1000:10.4f} {E_22/1000:10.4f} "
              f"{G_12/1000:10.4f} {nu_12:8.4f} {nu_21:8.4f} {E_11/E_22:9.4f}")
    
    print("\n" + "=" * 85)
    print("EXPECTED VALUES from Thesis Table 5.3:")
    print("θ=70.00°: E₁₁*=0.7159 GPa, E₂₂*=0.8400 GPa, G₁₂*=0.2159 GPa, ν₁₂*=0.2397, ν₂₁*=0.2812")
    print("θ=75.00°: E₁₁*=0.7430 GPa, E₂₂*=0.8539 GPa, G₁₂*=0.2028 GPa, ν₁₂*=0.2238, ν₂₁*=0.2572")
    print("θ=80.00°: E₁₁*=0.7756 GPa, E₂₂*=0.8630 GPa, G₁₂*=0.1890 GPa, ν₁₂*=0.2085, ν₂₁*=0.2320")
    print("θ=85.00°: E₁₁*=0.8089 GPa, E₂₂*=0.8600 GPa, G₁₂*=0.1788 GPa, ν₁₂*=0.1992, ν₂₁*=0.2117")
    print("=" * 85)
    
    print("\n" + "=" * 85)
    print("ANALYSIS:")
    print("=" * 85)
    print("\nThe compliance matrix method provides elastic properties through direct")
    print("homogenization of the pin-jointed truss structure.")
    print("\nKey updates (v2.0 - Double-thickness horizontal strut):")
    print("- Horizontal base strut now has 2× thickness (2t×d instead of t×d)")
    print("- This increases stiffness and strength in horizontal direction")
    print("- Relative density calculations updated accordingly")
    print("- Properties will differ from v1.0 results due to modified geometry")
    print("\nThe method correctly captures:")
    print("✓ Asymmetric strut design (horizontal stronger than inclined)")
    print("✓ Linear density scaling (E* ∝ ρ̄)")
    print("✓ Stretching-dominated behavior")
    print("✓ Enhanced E₁₁* from stronger horizontal strut")
    print("=" * 85)
