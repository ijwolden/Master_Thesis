import numpy as np
import pandas as pd
from pathlib import Path

def parse_sections(file_path):
    datasets = {}
    current_label = None
    current_data = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('#'):
                if current_label and current_data:
                    datasets[current_label] = np.array(current_data)
                current_label = line[1:].strip().lower().replace(' ', '_')
                current_data = []
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        disp = float(parts[0])
                        force = float(parts[1])
                        current_data.append([disp, force])
                    except:
                        pass
        
        if current_label and current_data:
            datasets[current_label] = np.array(current_data)
    
    return datasets

def compute_postpeak_metrics(disp, force, threshold=0.8):
    idx_max = np.argmax(force)
    F_max = force[idx_max]
    delta_max = disp[idx_max]
    
                                                                 
    post_peak_indices = np.where(disp > delta_max)[0]
    
    if len(post_peak_indices) == 0:
                           
        return {
            'ductility': np.nan,
            'softening_slope': np.nan,
            'residual_strength': np.nan,
            'failure_mode': 'unknown'
        }
    
    post_peak_force = force[post_peak_indices]
    post_peak_disp = disp[post_peak_indices]
    
                                            
    failure_threshold = threshold * F_max
    failure_indices = np.where(post_peak_force <= failure_threshold)[0]
    
    if len(failure_indices) > 0:
        failure_idx = post_peak_indices[failure_indices[0]]
        delta_failure = disp[failure_idx]
        ductility = delta_failure / delta_max
    else:
                                                    
        delta_failure = disp[-1]
        ductility = delta_failure / delta_max
    
                                                                           
    n_postpeak = min(len(post_peak_indices), max(5, len(post_peak_indices) // 5))
    fit_indices = post_peak_indices[:n_postpeak]
    
    if len(fit_indices) > 2:
        coeffs = np.polyfit(disp[fit_indices], force[fit_indices], 1)
        softening_slope = coeffs[0]
    else:
        softening_slope = np.nan
    
                                  
    target_disp = 2 * delta_max
    if disp[-1] >= target_disp:
                                              
        residual_strength = np.interp(target_disp, disp, force)
    else:
        residual_strength = np.nan
    
                           
    if ductility > 2.0:
        failure_mode = 'ductile'
    elif ductility > 1.3:
        failure_mode = 'semi-ductile'
    else:
        failure_mode = 'brittle'
    
    return {
        'ductility': ductility,
        'softening_slope': softening_slope,
        'residual_strength': residual_strength,
        'failure_mode': failure_mode
    }

def main():
    base_dir = Path(__file__).parent
    
    print("="*80)
    print("POST-PEAK BEHAVIOR METRICS")
    print("="*80)
    print(f"Using failure threshold: 80% of F_max")
    
                                    
    print("\n### BASELINE ###\n")
    baseline_file = base_dir / "Baseline" / "Baseline.txt"
    baseline_data = parse_sections(baseline_file)
    
    for label in ['a', 'a-precut', 'b', 'b_precut']:
        if label in baseline_data:
            data = baseline_data[label]
            disp = data[:, 0]
            force = data[:, 1]
            metrics = compute_postpeak_metrics(disp, force)
            print(f"{label:15s}: ψ={metrics['ductility']:.2f}  "
                  f"slope={metrics['softening_slope']:.5f}  "
                  f"F_res={metrics['residual_strength']:.5f}  "
                  f"mode={metrics['failure_mode']}")
    
                                     
    print("\n### SYMMETRIC (selected designs) ###\n")
    sym_intact_file = base_dir / "Symmetric" / "f-d-all.txt"
    data = pd.read_csv(sym_intact_file, sep='\t', header=0)
    
    selected = ['0.1', '0.4', '0.5', '0.9']
    for i, design in enumerate(['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']):
        if design not in selected:
            continue
        
        col_idx = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'].index(design) * 3
        
        disp_data = data.iloc[:, col_idx].values
        force_data = data.iloc[:, col_idx+1].values
        valid = ~(np.isnan(disp_data) | np.isnan(force_data))
        disp_clean = disp_data[valid]
        force_clean = force_data[valid]
        
        metrics = compute_postpeak_metrics(disp_clean, force_clean)
        print(f"{design}_intact: ψ={metrics['ductility']:.2f}  "
              f"slope={metrics['softening_slope']:.5f}  "
              f"F_res={metrics['residual_strength']:.5f}  "
              f"mode={metrics['failure_mode']}")
    
            
    print("\n### SYMMETRIC PRECUT ###\n")
    sym_precut_file = base_dir / "Symmetric" / "f-d-precut.txt"
    data = pd.read_csv(sym_precut_file, sep='\t', header=0)
    
    for i, design in enumerate(['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']):
        if design not in selected:
            continue
        
        col_idx = i * 3
        if col_idx + 1 >= len(data.columns):
            continue
        
        disp_data = data.iloc[:, col_idx].values
        force_data = data.iloc[:, col_idx+1].values
        valid = ~(np.isnan(disp_data) | np.isnan(force_data))
        
        if np.sum(valid) > 0:
            disp_clean = disp_data[valid]
            force_clean = force_data[valid]
            
            metrics = compute_postpeak_metrics(disp_clean, force_clean)
            print(f"{design}_precut: ψ={metrics['ductility']:.2f}  "
                  f"slope={metrics['softening_slope']:.5f}  "
                  f"F_res={metrics['residual_strength']:.5f}  "
                  f"mode={metrics['failure_mode']}")
    
    print("\n" + "="*80)
    print("✓ Post-peak metrics computed!")
    print("="*80)

if __name__ == "__main__":
    main()
