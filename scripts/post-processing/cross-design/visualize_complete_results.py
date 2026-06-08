import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'serif'
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.titlesize'] = 13

          
baseline = {'A': 0.30481, 'B': 0.30471}
baseline_eta = {'A': 1.424, 'B': 1.259}

symmetric = {f'0.{i}': k for i, k in enumerate([0.46598, 0.45468, 0.44428, 0.39127, 
                                                  0.42561, 0.42370, 0.42418, 0.42374, 0.47386], 1)}
symmetric_eta = {'0.1': 1.195, '0.4': 1.652, '0.5': 1.315, '0.9': 1.178}

g1 = {'A': 0.32688, 'B': 0.29174, 'C': 0.44812}            
g1_eta = {'A': 1.628, 'B': 1.346, 'C': 1.660}

g2 = {'A': 0.28922, 'B': 0.58606}
g2_eta = {'A': 1.505, 'B': 1.043}

g3 = {'G3': 0.44073}
g3_eta = {'G3': 1.413}

                                                                               
fig, ax = plt.subplots(figsize=(14, 6))

x_pos = 0
colors = {'Baseline': '#E63946', 'Symmetric': '#457B9D', 'G1': '#2A9D8F', 
          'G2': '#E9C46A', 'G3': '#F4A261'}

          
x_baseline = [x_pos, x_pos + 1]
ax.bar(x_baseline, list(baseline.values()), width=0.6, label='Baseline', 
       color=colors['Baseline'], alpha=0.8, edgecolor='black', linewidth=0.5)
x_pos += 3

           
x_symmetric = np.arange(x_pos, x_pos + len(symmetric))
ax.bar(x_symmetric, list(symmetric.values()), width=0.6, label='Symmetric',
       color=colors['Symmetric'], alpha=0.8, edgecolor='black', linewidth=0.5)
x_pos += len(symmetric) + 2

    
x_g1 = np.arange(x_pos, x_pos + len(g1))
ax.bar(x_g1, list(g1.values()), width=0.6, label='G1',
       color=colors['G1'], alpha=0.8, edgecolor='black', linewidth=0.5)
x_pos += len(g1) + 2

    
x_g2 = np.arange(x_pos, x_pos + len(g2))
ax.bar(x_g2, list(g2.values()), width=0.6, label='G2',
       color=colors['G2'], alpha=0.8, edgecolor='black', linewidth=0.5)
x_pos += len(g2) + 2

    
x_g3 = [x_pos]
ax.bar(x_g3, list(g3.values()), width=0.6, label='G3',
       color=colors['G3'], alpha=0.8, edgecolor='black', linewidth=0.5)

        
all_labels = ['Base-A', 'Base-B'] + [f'Sym-{k}' for k in symmetric.keys()] +\
             [f'G1-{k}' for k in g1.keys()] + [f'G2-{k}' for k in g2.keys()] + ['G3']
all_x = list(x_baseline) + list(x_symmetric) + list(x_g1) + list(x_g2) + list(x_g3)

ax.set_xticks(all_x)
ax.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=13)
ax.set_ylabel('Initial Stiffness $k$ [kN/mm]', fontsize=15)
ax.set_xlabel('Design', fontsize=15)
ax.set_title('Initial Stiffness Comparison', fontsize=15)
ax.tick_params(axis='y', labelsize=13)
ax.legend(loc='upper left', ncol=5, fontsize=13)
ax.grid(axis='y', alpha=0.3)
ax.axhline(y=np.mean(list(baseline.values())), color=colors['Baseline'], 
           linestyle='--', alpha=0.5, linewidth=1, label='_baseline_mean')

plt.tight_layout()
plt.savefig('complete_stiffness_all_families.png', dpi=300, bbox_inches='tight')
print("✓ Saved: complete_stiffness_all_families.png")
plt.close()

                                                                               
fig, ax = plt.subplots(figsize=(10, 7))

                        
all_eta_labels = ([f'Base-{k}' for k in baseline_eta.keys()] +
                  [f'Sym-{k}' for k in symmetric_eta.keys()] +
                  [f'G1-{k}' for k in g1_eta.keys()] +
                  [f'G2-{k}' for k in g2_eta.keys()] +
                  [f'G3'])

all_eta_values = (list(baseline_eta.values()) + list(symmetric_eta.values()) +
                  list(g1_eta.values()) + list(g2_eta.values()) + list(g3_eta.values()))

                   
sorted_pairs = sorted(zip(all_eta_labels, all_eta_values), key=lambda x: x[1])
sorted_labels, sorted_eta = zip(*sorted_pairs)

                 
family_colors = []
for label in sorted_labels:
    if 'Base' in label:
        family_colors.append(colors['Baseline'])
    elif 'Sym' in label:
        family_colors.append(colors['Symmetric'])
    elif 'G1' in label:
        family_colors.append(colors['G1'])
    elif 'G2' in label:
        family_colors.append(colors['G2'])
    else:
        family_colors.append(colors['G3'])

bars = ax.barh(range(len(sorted_eta)), sorted_eta, color=family_colors, 
               alpha=0.8, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(sorted_eta)))
ax.set_yticklabels(sorted_labels, fontsize=13)
ax.set_xlabel('Notch Sensitivity $\\eta = F_{intact} / F_{precut}$', fontsize=15)
ax.set_title('Fracture Tolerance', fontsize=15)
ax.axvline(1.0, color='green', linestyle='--', linewidth=2, alpha=0.6, label='Perfect tolerance')
ax.axvline(1.2, color='orange', linestyle=':', linewidth=1.5, alpha=0.5, label='Good tolerance')
ax.grid(axis='x', alpha=0.3)
ax.legend(fontsize=13)
ax.tick_params(axis='x', labelsize=13)

                  
for i, val in enumerate(sorted_eta):
    ax.text(val + 0.02, i, f'{val:.3f}', va='center', fontsize=11)

                          
ax.get_yticklabels()[0].set_weight('bold')
ax.get_yticklabels()[0].set_color('green')
ax.get_yticklabels()[-1].set_weight('bold')
ax.get_yticklabels()[-1].set_color('red')

plt.tight_layout()
plt.savefig('complete_notch_sensitivity_all_families.png', dpi=300, bbox_inches='tight')
print("✓ Saved: complete_notch_sensitivity_all_families.png")
plt.close()

                                                                               
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle('Family-by-Family Mechanical Performance', fontsize=14, fontweight='bold')

          
ax = axes[0, 0]
ax.bar(baseline.keys(), baseline.values(), color=colors['Baseline'], alpha=0.7)
ax.set_title('Baseline', fontweight='bold')
ax.set_ylabel('Stiffness $k$ [kN/mm]')
ax.grid(axis='y', alpha=0.3)

           
ax = axes[0, 1]
ax.bar(symmetric.keys(), symmetric.values(), color=colors['Symmetric'], alpha=0.7)
ax.set_title('Symmetric (9 variants)', fontweight='bold')
ax.set_xticklabels(symmetric.keys(), rotation=45, ha='right')
ax.grid(axis='y', alpha=0.3)

    
ax = axes[0, 2]
ax.bar(g1.keys(), g1.values(), color=colors['G1'], alpha=0.7)
ax.set_title('G1 (Loading Rate Study)', fontweight='bold')
ax.set_ylabel('Stiffness $k$ [kN/mm]')
ax.grid(axis='y', alpha=0.3)

    
ax = axes[1, 0]
ax.bar(g2.keys(), g2.values(), color=colors['G2'], alpha=0.7)
ax.set_title('G2 (Grading)', fontweight='bold')
ax.set_ylabel('Stiffness $k$ [kN/mm]')
ax.grid(axis='y', alpha=0.3)

    
ax = axes[1, 1]
ax.bar(g3.keys(), g3.values(), color=colors['G3'], alpha=0.7)
ax.set_title('G3 (Alternative Grading)', fontweight='bold')
ax.set_ylabel('Stiffness $k$ [kN/mm]')
ax.grid(axis='y', alpha=0.3)

                    
ax = axes[1, 2]
families_mean_k = {
    'Baseline': np.mean(list(baseline.values())),
    'Symmetric': np.mean(list(symmetric.values())),
    'G1': np.mean(list(g1.values())),
    'G2': np.mean(list(g2.values())),
    'G3': np.mean(list(g3.values()))
}
family_colors_list = [colors[f] for f in families_mean_k.keys()]
ax.bar(families_mean_k.keys(), families_mean_k.values(), 
       color=family_colors_list, alpha=0.7, edgecolor='black', linewidth=1)
ax.set_title('Mean Stiffness by Family', fontweight='bold')
ax.set_ylabel('Mean $k$ [kN/mm]')
ax.set_xticklabels(families_mean_k.keys(), rotation=30, ha='right')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('complete_family_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Saved: complete_family_comparison.png")
plt.close()

                                                                               
fig, ax = plt.subplots(figsize=(10, 7))

                                         
families = {
    'Baseline': (baseline, baseline_eta, 's', colors['Baseline']),
    'Symmetric': (symmetric, symmetric_eta, 'o', colors['Symmetric']),
    'G1': (g1, g1_eta, '^', colors['G1']),
    'G2': (g2, g2_eta, 'D', colors['G2']),
    'G3': (g3, g3_eta, 'v', colors['G3'])
}

for family_name, (k_dict, eta_dict, marker, color) in families.items():
    k_vals = [k_dict[key] for key in eta_dict.keys() if key in k_dict]
    eta_vals = [eta_dict[key] for key in eta_dict.keys()]
    
    ax.scatter(k_vals, eta_vals, s=150, marker=marker, color=color,
              alpha=0.7, edgecolors='black', linewidth=1.5,
              label=family_name, zorder=3)
    
                
    for key in eta_dict.keys():
        if key in k_dict:
            ax.annotate(f'{family_name[:3]}-{key}', 
                       (k_dict[key], eta_dict[key]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)

ax.axhline(1.0, color='green', linestyle='--', linewidth=2, alpha=0.5,
          label='Perfect fracture tolerance')
ax.axhline(1.2, color='orange', linestyle=':', linewidth=1.5, alpha=0.5,
          label='Good fracture tolerance')

ax.set_xlabel('Initial Stiffness $k$ [kN/mm]')
ax.set_ylabel('Notch Sensitivity $\\eta$')
ax.set_title('Stiffness vs Fracture Tolerance Trade-off')
ax.legend(loc='upper right', fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('complete_stiffness_vs_eta_scatter.png', dpi=300, bbox_inches='tight')
print("✓ Saved: complete_stiffness_vs_eta_scatter.png")
plt.close()

                                                                               
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

                
ax = axes[0]
top_5_stiff = sorted([(k, f'G2-{key}') if 'G2' in str(type(g2)) else (k, f'Sym-{key}')
                      for d, color in [(baseline, 'Baseline'), (symmetric, 'Symmetric'),
                                       (g1, 'G1'), (g2, 'G2'), (g3, 'G3')]
                      for key, k in d.items()], reverse=True)[:5]

                      
all_k = []
all_k_labels = []
for key, val in baseline.items():
    all_k.append((val, f'Base-{key}'))
for key, val in symmetric.items():
    all_k.append((val, f'Sym-{key}'))
for key, val in g1.items():
    all_k.append((val, f'G1-{key}'))
for key, val in g2.items():
    all_k.append((val, f'G2-{key}'))
for key, val in g3.items():
    all_k.append((val, f'G3'))

top_5_stiff = sorted(all_k, reverse=True)[:5]
labels, values = zip(*[(label, k) for k, label in top_5_stiff])
ax.barh(range(5), values, color='#06A77D', alpha=0.8, edgecolor='black', linewidth=1)
ax.set_yticks(range(5))
ax.set_yticklabels(labels)
ax.set_xlabel('Stiffness $k$ [kN/mm]')
ax.set_title('Top 5: Highest Stiffness', fontweight='bold')
ax.grid(axis='x', alpha=0.3)
for i, v in enumerate(values):
    ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=9)

                         
ax = axes[1]
top_5_eta = sorted(all_eta_values)[:5]
top_5_eta_labels = [l for _, l in sorted(zip(all_eta_values, all_eta_labels))][:5]
ax.barh(range(5), top_5_eta, color='#F4A261', alpha=0.8, edgecolor='black', linewidth=1)
ax.set_yticks(range(5))
ax.set_yticklabels(top_5_eta_labels)
ax.set_xlabel('Notch Sensitivity $\\eta$')
ax.set_title('Top 5 Best Fracture Tolerance', fontweight='bold')
ax.grid(axis='x', alpha=0.3)
for i, v in enumerate(top_5_eta):
    ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=9)

                          
ax = axes[2]
worst_5_eta = sorted(all_eta_values, reverse=True)[:5]
worst_5_eta_labels = [l for _, l in sorted(zip(all_eta_values, all_eta_labels), reverse=True)][:5]
ax.barh(range(5), worst_5_eta, color='#D81159', alpha=0.8, edgecolor='black', linewidth=1)
ax.set_yticks(range(5))
ax.set_yticklabels(worst_5_eta_labels)
ax.set_xlabel('Notch Sensitivity $\\eta$')
ax.set_title('Bottom 5 Worst Fracture Tolerance', fontweight='bold')
ax.grid(axis='x', alpha=0.3)
for i, v in enumerate(worst_5_eta):
    ax.text(v + 0.02, i, f'{v:.3f}', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('complete_key_findings_summary.png', dpi=300, bbox_inches='tight')
print("✓ Saved: complete_key_findings_summary.png")
plt.close()

print("\n" + "="*80)
print("✓ All comprehensive visualizations created!")
print("="*80)
print("\nGenerated files:")
print("  1. complete_stiffness_all_families.png")
print("  2. complete_notch_sensitivity_all_families.png")
print("  3. complete_family_comparison.png")
print("  4. complete_stiffness_vs_eta_scatter.png")
print("  5. complete_key_findings_summary.png")
