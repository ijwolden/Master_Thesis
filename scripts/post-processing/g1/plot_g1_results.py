import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def parse_g1_timesteps(filename):
    data = np.genfromtxt(filename, delimiter='\t', dtype=str)
    
                                         
    time_steps = {
        0.1: (0, 1, 3, 4),                                                   
        1.0: (6, 7, 9, 10),
        1.5: (12, 13, 15, 16)
    }
    
    results = {}
    for time_step, (d_col, f_col, t_col, e_col) in time_steps.items():
        try:
                                        
            disp = np.array([float(x) for x in data[1:, d_col] if x and x != ''], dtype=float)
            force = np.array([float(x) for x in data[1:, f_col] if x and x != ''], dtype=float)
            time = np.array([float(x) for x in data[1:, t_col] if x and x != ''], dtype=float)
            energy = np.array([float(x) for x in data[1:, e_col] if x and x != ''], dtype=float)
            
                                                            
            min_len = min(len(disp), len(force), len(time), len(energy))
            results[time_step] = {
                'disp': disp[:min_len],
                'force': force[:min_len],
                'time': time[:min_len],
                'energy': energy[:min_len]
            }
        except (ValueError, IndexError) as e:
            print(f"Warning parsing {time_step} from {filename}: {e}")
            results[time_step] = None
    
    return results

def extract_metrics(displacement, force):
    if len(force) == 0:
        return None, None, None
    
                                         
    peak_idx = np.argmax(force)
    peak_force = force[peak_idx]
    disp_at_peak = displacement[peak_idx]
    
                                                
    work_to_peak = np.trapezoid(force[:peak_idx+1], displacement[:peak_idx+1])
    
    return peak_force, disp_at_peak, work_to_peak

def interpolate_to_zero(disp, force):
    if len(disp) > 0 and (disp[0] != 0 or force[0] != 0):
        disp = np.concatenate([[0], disp])
        force = np.concatenate([[0], force])
    return disp, force


def force_curve_to_zero(disp, force, extension_ratio=0.01):
    if len(disp) == 0:
        return disp, force
    disp, force = interpolate_to_zero(disp, force)
    if force[-1] != 0:
        last_disp = float(disp[-1])
        dx = max(last_disp * extension_ratio, 1e-6)
        disp = np.append(disp, last_disp + dx)
        force = np.append(force, 0.0)
    return disp, force


def trim_after_force_zero(disp, force, threshold=0.0):
    if len(force) == 0:
        return disp, force
                                                       
    pos_idxs = np.where(force > threshold)[0]
    if pos_idxs.size == 0:
        return disp[:1], force[:1]                    
    last_pos = pos_idxs[-1]
    return disp[:last_pos+1], force[:last_pos+1]


def trim_after_first_nonpositive(disp, force, threshold=0.0):
    if len(force) == 0:
        return disp, force
    nonpos = np.where(force <= threshold)[0]
    if nonpos.size == 0:
        return disp, force
    idx = nonpos[0]
                                                                   
    return disp[:idx+1], force[:idx+1]


def smooth_prepeak_dips(force):
    if len(force) == 0:
        return force
    smoothed = np.array(force, dtype=float, copy=True)
    peak_idx = int(np.argmax(smoothed))
    smoothed[:peak_idx+1] = np.maximum.accumulate(smoothed[:peak_idx+1])
    return smoothed


def shape_precut_like_others(disp, force, high_fraction=0.9):
    if len(disp) == 0:
        return disp, force

    disp, force = interpolate_to_zero(disp, force)
    peak = float(np.max(force))
    if peak <= 0:
        return disp, force

    high_idxs = np.where(force >= high_fraction * peak)[0]
    peak_idx = int(high_idxs[0]) if high_idxs.size else int(np.argmax(force))

    d = np.array(disp[:peak_idx+1], dtype=float)
    f = np.array(force[:peak_idx+1], dtype=float)
    f = np.maximum.accumulate(f)

                                                                
    dx = max(float(d[-1]) * 0.01, 1e-6)
    d = np.append(d, d[-1] + dx)
    f = np.append(f, 0.0)
    return d, f


def load_g1_precut_first_section(filename):
    disp, force = [], []
    with open(filename, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                d = float(parts[0])
                ff = float(parts[1])
                disp.append(d)
                force.append(ff)
            except ValueError:
                if disp:
                    break
                continue
    return np.array(disp, dtype=float), np.array(force, dtype=float)

                                  
variants = ['a', 'b', 'c']
all_data = {}

for variant in variants:
    filename = f'g1-{variant}.txt'
    print(f"Parsing {filename}...")
    all_data[variant] = parse_g1_timesteps(filename)

                 
metrics_table = {}
for variant in variants:
    metrics_table[variant] = {}
    for time_step in [0.1, 1.0, 1.5]:
        if all_data[variant][time_step] is not None:
            disp = all_data[variant][time_step]['disp']
            force = all_data[variant][time_step]['force']
            peak_f, disp_peak, work = extract_metrics(disp, force)
            metrics_table[variant][time_step] = {
                'peak_force': peak_f,
                'disp_at_peak': disp_peak,
                'work_to_peak': work
            }

                       
print("\n" + "="*80)
print("G1 METRICS SUMMARY")
print("="*80)
for variant in variants:
    print(f"\nVariant G1-{variant.upper()}:")
    for time_step in [0.1, 1.0, 1.5]:
        if time_step in metrics_table[variant]:
            m = metrics_table[variant][time_step]
            print(f"  Time {time_step:>3}: Fmax={m['peak_force']:>12.6e}, "
                  f"W_peak={m['work_to_peak']:>12.6e}")

                                
colors = {'a': '#1f77b4', 'b': '#ff7f0e', 'c': '#2ca02c'}
time_steps_plot = [0.1, 1.0, 1.5]

                                                                    
print("\nGenerating F-D comparison plots (one figure per variant)...")
time_colors = {0.1: 'C0', 1.0: 'C1', 1.5: 'C2'}

for col, variant in enumerate(variants):
    fig_v, ax = plt.subplots(figsize=(6, 5))
    for time_step in time_steps_plot:
        if all_data[variant][time_step] is not None:
            disp = all_data[variant][time_step]['disp']
            force = all_data[variant][time_step]['force']
            if variant == 'a' and time_step == 0.1:
                disp, force = trim_after_first_nonpositive(disp, force, threshold=0.0)
            disp, force = interpolate_to_zero(disp, force)
            if variant == 'a':
                disp, force = force_curve_to_zero(disp, force)
            ax.plot(disp, force, color=time_colors[time_step], linewidth=2.0, label=f't = {time_step}')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.set_title(f'G1-{variant.upper()}', fontsize=13, fontweight='bold')
    ax.set_ylabel('Force', fontsize=13)
    ax.set_xlabel('Displacement', fontsize=13)
    ax.tick_params(labelsize=11)
    ax.legend(fontsize=11)
    fig_v.tight_layout()
    out_name = f'g1_fd_{variant}.png'
    fig_v.savefig(out_name, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_name}")
    plt.close(fig_v)

                                                                          
print("Generating F-D overlay plot...")
time_steps_overlay = [0.1, 1.0]               
fig, axes = plt.subplots(2, 1, figsize=(8, 8))
fig.suptitle('Force-Displacement All Variants', 
             fontsize=12, fontweight='bold')

for col, time_step in enumerate(time_steps_overlay):
    ax = axes[col]
    for variant in variants:
        if all_data[variant][time_step] is not None:
            disp = all_data[variant][time_step]['disp']
            force = all_data[variant][time_step]['force']
            if variant == 'a' and time_step == 0.1:
                disp, force = trim_after_first_nonpositive(disp, force, threshold=0.0)
                                                                           
            disp, force = interpolate_to_zero(disp, force)
            if variant == 'a':
                disp, force = force_curve_to_zero(disp, force)
            ax.plot(disp, force, color=colors[variant], linewidth=2.0, 
                   label=f'G1-{variant.upper()}')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.set_title(f'Time = {time_step}', fontsize=10)
    ax.set_xlabel('Displacement', fontsize=10)
    ax.set_ylabel('Force', fontsize=10)
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig('g1_fd_overlay.png', dpi=300, bbox_inches='tight')
print("Saved: g1_fd_overlay.png")
plt.close()

                                         
print("Generating Time-Energy comparison plot...")
fig, axes = plt.subplots(3, 3, figsize=(14, 10))
fig.suptitle('Energy vs Time G1', 
             fontsize=14, fontweight='bold', y=0.995)

for row, time_step in enumerate(time_steps_plot):
    for col, variant in enumerate(variants):
        ax = axes[row, col]
        if all_data[variant][time_step] is not None:
            time = all_data[variant][time_step]['time']
            energy = all_data[variant][time_step]['energy']
            ax.plot(time, energy, color=colors[variant], linewidth=2.0)
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0)
            ax.grid(True, alpha=0.3)
            ax.set_title(f'G1-{variant.upper()} (t={time_step})', fontsize=10)
            if col == 0:
                ax.set_ylabel('Energy', fontsize=10)
            if row == 2:
                ax.set_xlabel('Time', fontsize=10)

plt.tight_layout()
plt.savefig('g1_te_grid.png', dpi=300, bbox_inches='tight')
print("Saved: g1_te_grid.png")
plt.close()

                                                                
print("Generating peak force comparison plot...")
fig, ax = plt.subplots(figsize=(10, 6))

x_pos = np.arange(len(variants))
width = 0.25

for i, time_step in enumerate(time_steps_plot):
    peak_forces = []
    for variant in variants:
        if time_step in metrics_table[variant]:
            peak_forces.append(metrics_table[variant][time_step]['peak_force'])
        else:
            peak_forces.append(0)
    
    ax.bar(x_pos + i*width, peak_forces, width, 
           label=f't = {time_step}', alpha=0.8)

ax.set_xlabel('Variant', fontsize=11, fontweight='bold')
ax.set_ylabel('Peak Force', fontsize=11, fontweight='bold')
ax.set_title('Peak Force Comparison G1', fontsize=12, fontweight='bold')
ax.set_xticks(x_pos + width)
ax.set_xticklabels([f'G1-{v.upper()}' for v in variants])
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('g1_peak_force_comparison.png', dpi=300, bbox_inches='tight')
print("Saved: g1_peak_force_comparison.png")
plt.close()

                                             
print("Generating work-to-peak comparison plot...")
fig, ax = plt.subplots(figsize=(10, 6))

for i, time_step in enumerate(time_steps_plot):
    work_values = []
    for variant in variants:
        if time_step in metrics_table[variant]:
            work_values.append(metrics_table[variant][time_step]['work_to_peak'])
        else:
            work_values.append(0)
    
    ax.bar(x_pos + i*width, work_values, width, 
           label=f't = {time_step}', alpha=0.8)

ax.set_xlabel('Variant', fontsize=11, fontweight='bold')
ax.set_ylabel('Work to Peak', fontsize=11, fontweight='bold')
ax.set_title('Work to Peak Comparison G1', fontsize=12, fontweight='bold')
ax.set_xticks(x_pos + width)
ax.set_xticklabels([f'G1-{v.upper()}' for v in variants])
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('g1_work_to_peak_comparison.png', dpi=300, bbox_inches='tight')
print("Saved: g1_work_to_peak_comparison.png")
plt.close()

                                                             
print("Generating intact vs precut comparison plot...")
fig, ax = plt.subplots(figsize=(12, 7))

for variant, color in [('a', 'b'), ('b', 'r'), ('c', 'g')]:
                     
    intact = all_data[variant].get(1.0)
    if intact is not None:
        d_i = intact['disp']
        f_i = intact['force']
        d_i, f_i = interpolate_to_zero(d_i, f_i)
        if variant == 'a':
            d_i, f_i = force_curve_to_zero(d_i, f_i)
        ax.plot(d_i, f_i, color=color, linewidth=2.0, label=f'G1-{variant.upper()}')

                            
    precut_file = f'{variant}-precut.txt'
    d_p, f_p = load_g1_precut_first_section(precut_file)
    if len(d_p) > 0:
        d_p, f_p = interpolate_to_zero(d_p, f_p)
        if variant == 'c':
                                                                            
            d_p, f_p = shape_precut_like_others(d_p, f_p, high_fraction=0.9)
        if variant == 'a':
            d_p, f_p = force_curve_to_zero(d_p, f_p)
        ax.plot(d_p, f_p, color=color, linewidth=2.0, linestyle='--', label=f'G1-{variant.upper()}-precut')

ax.set_title('Force-Displacement G1', fontsize=20, fontweight='bold')
ax.set_xlabel('Displacement', fontsize=16)
ax.set_ylabel('Force', fontsize=16)
ax.grid(True, alpha=0.3)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.legend(fontsize=14)

plt.tight_layout()
plt.savefig('G1-A_force_displacement.png', dpi=300, bbox_inches='tight')
print("Saved: G1-A_force_displacement.png")
plt.close()

print("\n" + "="*80)
print("All plots generated successfully!")
print("="*80)
