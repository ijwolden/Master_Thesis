import numpy as np
import pandas as pd
from pathlib import Path
import json

def load_force_displacement(filepath):
    data = pd.read_csv(filepath, sep='\t', header=0)
                                                 
    disp = data.iloc[:, 0].values
    force = data.iloc[:, 1].values
                                      
    valid = ~(np.isnan(disp) | np.isnan(force))
    return disp[valid], force[valid]

def compute_initial_stiffness(disp, force, fraction=0.1):
               
    idx_max = np.argmax(force)
    disp_max = disp[idx_max]
    
                                           
    linear_region = disp <= fraction * disp_max
    
    if np.sum(linear_region) < 3:
        return np.nan
    
                       
    coeffs = np.polyfit(disp[linear_region], force[linear_region], 1)
    stiffness = coeffs[0]         
    
    return stiffness

def compute_peak_metrics(disp, force):
    idx_max = np.argmax(force)
    F_max = force[idx_max]
    delta_max = disp[idx_max]
    
                                                
    W_peak = np.trapezoid(force[:idx_max+1], disp[:idx_max+1])
    
    return F_max, delta_max, W_peak

def compute_ductility_index(disp, force, threshold=0.5):
    idx_max = np.argmax(force)
    F_max = force[idx_max]
    delta_max = disp[idx_max]
    
                      
    post_peak = disp > delta_max
    if not np.any(post_peak):
        return np.nan
    
    disp_post = disp[post_peak]
    force_post = force[post_peak]
    
                                            
    below_threshold = force_post < threshold * F_max
    if not np.any(below_threshold):
                                                            
        return np.nan
    
    idx_failure = np.where(below_threshold)[0][0]
    delta_failure = disp_post[idx_failure]
    
    psi = delta_failure / delta_max
    return psi

def compute_softening_slope(disp, force):
    idx_max = np.argmax(force)
    F_max = force[idx_max]
    
                      
    post_peak_mask = (disp > disp[idx_max]) & (force > 0.5 * F_max) & (force < F_max)
    
    if np.sum(post_peak_mask) < 3:
        return np.nan
    
    disp_post = disp[post_peak_mask]
    force_post = force[post_peak_mask]
    
                       
    coeffs = np.polyfit(disp_post, force_post, 1)
    softening_slope = coeffs[0]                                      
    
    return softening_slope

def compute_residual_strength(disp, force):
    idx_max = np.argmax(force)
    delta_max = disp[idx_max]
    target_disp = 2.0 * delta_max
    
                                     
    if disp[-1] < target_disp:
        return np.nan
    
    idx_residual = np.argmin(np.abs(disp - target_disp))
    residual_force = force[idx_residual]
    
                        
    return residual_force / force[idx_max]


def process_baseline():
    base_dir = Path(__file__).parent / "Baseline"
    
    results = {}
    
                  
    files = {
        'A': 'BaselineA.txt',
        'B': 'BaselineB.txt'
    }
    
    for case, filename in files.items():
        filepath = base_dir / filename
        if not filepath.exists():
            print(f"Warning: {filepath} not found")
            continue
        
        disp, force = load_force_displacement(filepath)
        
                                                                                       
        data = pd.read_csv(filepath, sep='\t', header=0)
        n_cols = len(data.columns)
        
                                              
        case_results = []
        for i in range(0, n_cols, 2):
            if i+1 >= n_cols:
                break
            
            disp = data.iloc[:, i].values
            force = data.iloc[:, i+1].values
            valid = ~(np.isnan(disp) | np.isnan(force))
            disp = disp[valid]
            force = force[valid]
            
            if len(disp) < 10:
                continue
            
            k = compute_initial_stiffness(disp, force)
            F_max, delta_max, W_peak = compute_peak_metrics(disp, force)
            psi = compute_ductility_index(disp, force)
            soft_slope = compute_softening_slope(disp, force)
            residual = compute_residual_strength(disp, force)
            
            case_results.append({
                'k': k,
                'F_max': F_max,
                'delta_max': delta_max,
                'W_peak': W_peak,
                'psi': psi,
                'softening_slope': soft_slope,
                'residual_strength': residual
            })
        
        results[case] = case_results
    
    return results


def process_symmetric():
    base_dir = Path(__file__).parent / "Symmetric"
    
    results = {}
    
                  
    intact_file = base_dir / "f-d-all.txt"
    if intact_file.exists():
        data = pd.read_csv(intact_file, sep='\t', header=0)
        
                                                     
        designs = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']
        
        for i, design in enumerate(designs):
            col_idx = i * 2                                           
            if col_idx + 1 >= len(data.columns):
                break
            
            disp = data.iloc[:, col_idx].values
            force = data.iloc[:, col_idx + 1].values
            valid = ~(np.isnan(disp) | np.isnan(force))
            disp = disp[valid]
            force = force[valid]
            
            if len(disp) < 10:
                continue
            
            k = compute_initial_stiffness(disp, force)
            F_max, delta_max, W_peak = compute_peak_metrics(disp, force)
            psi = compute_ductility_index(disp, force)
            soft_slope = compute_softening_slope(disp, force)
            
            results[f'{design}_intact'] = {
                'k': k,
                'F_max': F_max,
                'delta_max': delta_max,
                'W_peak': W_peak,
                'psi': psi,
                'softening_slope': soft_slope
            }
    
                  
    precut_file = base_dir / "f-d-precut.txt"
    if precut_file.exists():
        data = pd.read_csv(precut_file, sep='\t', header=0)
        
        designs = ['0.1', '0.4', '0.5', '0.9']                                     
        
        for i, design in enumerate(designs):
            col_idx = i * 2
            if col_idx + 1 >= len(data.columns):
                break
            
            disp = data.iloc[:, col_idx].values
            force = data.iloc[:, col_idx + 1].values
            valid = ~(np.isnan(disp) | np.isnan(force))
            disp = disp[valid]
            force = force[valid]
            
            if len(disp) < 10:
                continue
            
            k = compute_initial_stiffness(disp, force)
            F_max, delta_max, W_peak = compute_peak_metrics(disp, force)
            psi = compute_ductility_index(disp, force)
            
            results[f'{design}_precut'] = {
                'k': k,
                'F_max': F_max,
                'delta_max': delta_max,
                'W_peak': W_peak,
                'psi': psi
            }
    
    return results


def compute_notch_sensitivity(results):
    eta_values = {}
    
              
    for case in ['A', 'B']:
        if case in results['baseline']:
            case_data = results['baseline'][case]
            if len(case_data) >= 2:                                       
                F_intact = case_data[0]['F_max']
                F_precut = case_data[1]['F_max']
                eta_values[f'Baseline_{case}'] = F_intact / F_precut
    
               
    sym_results = results.get('symmetric', {})
    for design in ['0.1', '0.4', '0.5', '0.9']:
        intact_key = f'{design}_intact'
        precut_key = f'{design}_precut'
        if intact_key in sym_results and precut_key in sym_results:
            F_intact = sym_results[intact_key]['F_max']
            F_precut = sym_results[precut_key]['F_max']
            eta_values[f'Symmetric_{design}'] = F_intact / F_precut
    
    return eta_values


def main():
    print("Computing all metrics...")
    
    all_results = {}
    
                         
    print("\n1. Processing Baseline...")
    all_results['baseline'] = process_baseline()
    
    print("\n2. Processing Symmetric...")
    all_results['symmetric'] = process_symmetric()
    
    print("\n3. Computing notch sensitivity ratios...")
    all_results['notch_sensitivity'] = compute_notch_sensitivity(all_results)
    
                  
    output_file = Path(__file__).parent / "computed_metrics.json"
    
                                                                       
    def convert_numpy(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(item) for item in obj]
        return obj
    
    all_results_clean = convert_numpy(all_results)
    
    with open(output_file, 'w') as f:
        json.dump(all_results_clean, f, indent=2)
    
    print(f"\n✓ Results saved to: {output_file}")
    
                   
    print("\n" + "="*60)
    print("SUMMARY OF COMPUTED METRICS")
    print("="*60)
    
    print("\nBASELINE:")
    for case, data in all_results['baseline'].items():
        print(f"\n  {case}:")
        for i, variant in enumerate(data):
            label = ['intact', 'precut', 'precut-2'][i] if i < 3 else f'variant_{i}'
            print(f"    {label}: k={variant['k']:.5f}, F_max={variant['F_max']:.5f}, psi={variant['psi']}")
    
    print("\nSYMMETRIC:")
    for design, data in all_results['symmetric'].items():
        print(f"  {design}: k={data['k']:.5f}, F_max={data['F_max']:.5f}, psi={data.get('psi', 'N/A')}")
    
    print("\nNOTCH SENSITIVITY:")
    for design, eta in all_results['notch_sensitivity'].items():
        print(f"  {design}: η = {eta:.3f}")
    
    print("\n" + "="*60)
    print("✓ All metrics computed successfully!")
    print("="*60)
    
    return all_results


if __name__ == "__main__":
    results = main()
