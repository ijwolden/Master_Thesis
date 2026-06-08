import numpy as np
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

def compute_stiffness(disp, force, fraction=0.1):
    idx_max = np.argmax(force)
    disp_max = disp[idx_max]
    
                   
    mask = disp <= fraction * disp_max
    if np.sum(mask) < 3:
        return np.nan
    
                
    coeffs = np.polyfit(disp[mask], force[mask], 1)
    return coeffs[0]

def main():
    base_dir = Path(__file__).parent
    
    print("="*70)
    print("COMPUTING INITIAL STIFFNESS FOR ALL DESIGNS")
    print("="*70)
    
                                    
    print("\n### BASELINE ###\n")
    baseline_file = base_dir / "Baseline" / "Baseline.txt"
    baseline_data = parse_sections(baseline_file)
    
    baseline_results = {}
    for label, data in baseline_data.items():
        disp = data[:, 0]
        force = data[:, 1]
        k = compute_stiffness(disp, force)
        F_max = np.max(force)
        baseline_results[label] = {'k': k, 'F_max': F_max}
        print(f"{label:15s}:  k = {k:.5f},  F_max = {F_max:.5f}")
    
                                     
    print("\n### SYMMETRIC ###\n")
    sym_intact_file = base_dir / "Symmetric" / "f-d-all.txt"
    sym_precut_file = base_dir / "Symmetric" / "f-d-precut.txt"
    
    sym_results = {}
    
                    
    with open(sym_intact_file, 'r') as f:
        lines = f.readlines()
    
                                            
    header = lines[0].strip().split('\t')
    designs_intact = [h.strip() for h in header if h.strip() and not h.strip().startswith('Unnamed')]
    
               
    data = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        row = []
        for p in parts:
            try:
                row.append(float(p))
            except:
                row.append(np.nan)
        data.append(row)
    data = np.array(data)
    
                                                        
    for i, design in enumerate(designs_intact):
        col_idx = i * 2
        if col_idx + 1 >= data.shape[1]:
            break
        
        disp = data[:, col_idx]
        force = data[:, col_idx + 1]
        
                    
        valid = ~(np.isnan(disp) | np.isnan(force))
        disp = disp[valid]
        force = force[valid]
        
        if len(disp) < 10:
            continue
        
        k = compute_stiffness(disp, force)
        F_max = np.max(force)
        sym_results[f'{design}_intact'] = {'k': k, 'F_max': F_max}
        print(f"{design}_intact:  k = {k:.5f},  F_max = {F_max:.5f}")
    
                    
    if sym_precut_file.exists():
        with open(sym_precut_file, 'r') as f:
            lines = f.readlines()
        
        header = lines[0].strip().split('\t')
        designs_precut = [h.strip() for h in header if h.strip() and not h.strip().startswith('Unnamed')]
        
        data = []
        for line in lines[1:]:
            parts = line.strip().split('\t')
            row = []
            for p in parts:
                try:
                    row.append(float(p))
                except:
                    row.append(np.nan)
            data.append(row)
        data = np.array(data)
        
        for i, design in enumerate(designs_precut):
            col_idx = i * 2
            if col_idx + 1 >= data.shape[1]:
                break
            
            disp = data[:, col_idx]
            force = data[:, col_idx + 1]
            
            valid = ~(np.isnan(disp) | np.isnan(force))
            disp = disp[valid]
            force = force[valid]
            
            if len(disp) < 10:
                continue
            
            k = compute_stiffness(disp, force)
            F_max = np.max(force)
            sym_results[f'{design}_precut'] = {'k': k, 'F_max': F_max}
            print(f"{design}_precut:  k = {k:.5f},  F_max = {F_max:.5f}")
    
                                             
    print("\n### NOTCH SENSITIVITY (η) ###\n")
    
                
    eta_a = baseline_results['a']['F_max'] / baseline_results['a-precut']['F_max']
    print(f"Baseline A:  η = {eta_a:.3f}")
    
                
    eta_b = baseline_results['b']['F_max'] / baseline_results['b_precut']['F_max']
    print(f"Baseline B:  η = {eta_b:.3f}")
    
               
    for design in ['0.1', '0.4', '0.5', '0.9']:
        intact_key = f'{design}_intact'
        precut_key = f'{design}_precut'
        if intact_key in sym_results and precut_key in sym_results:
            eta = sym_results[intact_key]['F_max'] / sym_results[precut_key]['F_max']
            print(f"Symmetric {design}:  η = {eta:.3f}")
    
    print("\n" + "="*70)
    print("✓ All stiffness values computed!")
    print("="*70)
    
                            
    return baseline_results, sym_results

if __name__ == "__main__":
    baseline, symmetric = main()
