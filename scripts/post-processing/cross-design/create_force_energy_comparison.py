import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

           
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150

                          
families_data = {
    'Baseline': {
        'designs': ['A', 'B'],
        'F_intact': [0.06319, 0.06416],
        'F_precut': [0.04437, 0.05095],
        'W_intact': [0.003633, 0.003890],
        'W_precut': [0.001545, 0.002159],
        'color': '#E63946'
    },
    'Symmetric': {
        'designs': ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9'],
        'F_intact': [0.05156, 0.05824, 0.06441, 0.07218, 0.07244, 0.07105, 0.06922, 0.06723, 0.06505],
        'F_precut': [0.04317, 0.04725, 0.05124, 0.04370, 0.05507, 0.05685, 0.05823, 0.06012, 0.05520],
        'W_intact': [0.001473, 0.001754, 0.002120, 0.002591, 0.002718, 0.002627, 0.002469, 0.002282, 0.002002],
        'W_precut': [0.001060, 0.001234, 0.001478, 0.000953, 0.001461, 0.001597, 0.001689, 0.001812, 0.001441],
        'color': '#457B9D'
    },
    'G1': {
        'designs': ['A', 'B', 'C'],
        'F_intact': [0.06411, 0.05482, 0.07130],            
        'F_precut': [0.03937, 0.04073, 0.04296],
        'W_intact': [0.007214, 0.006144, 0.009249],
        'W_precut': [0.002021, 0.002227, 0.002386],
        'color': '#2A9D8F'
    },
    'G2': {
        'designs': ['A', 'B'],
        'F_intact': [0.03940, 0.06719],
        'F_precut': [0.02617, 0.06439],
        'W_intact': [0.001292, 0.004502],
        'W_precut': [0.000541, 0.003641],
        'color': '#E9C46A'
    },
    'G3': {
        'designs': ['G3'],
        'F_intact': [0.05950],
        'F_precut': [0.04211],
        'W_intact': [0.003108],
        'W_precut': [0.001629],
        'color': '#F4A261'
    }
}

                               
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

                           
x_positions = []
intact_forces = []
precut_forces = []
intact_energies = []
precut_energies = []
colors = []
labels = []

x_pos = 0
family_boundaries = []
family_positions = []
family_labels_list = []

for family_name, data in families_data.items():
    n_designs = len(data['designs'])
    
                                                
    family_start = x_pos
    
    for i, design in enumerate(data['designs']):
        x_positions.append(x_pos)
        intact_forces.append(data['F_intact'][i])
        precut_forces.append(data['F_precut'][i])
        intact_energies.append(data['W_intact'][i])
        precut_energies.append(data['W_precut'][i])
        colors.append(data['color'])
        
                      
        if family_name == 'Symmetric':
            labels.append(design)
        elif family_name in ['G1', 'G2']:
            labels.append(f"{family_name}-{design}")
        elif family_name == 'G3':
            labels.append('G3')
        else:
            labels.append(f"{family_name}\n{design}")
        
        x_pos += 1
    
                                            
    family_center = (family_start + x_pos - 1) / 2
    family_positions.append(family_center)
    family_labels_list.append(family_name)
    
    family_boundaries.append(x_pos - 0.5)
    x_pos += 0.5                        

x_positions = np.array(x_positions)

                    
width = 0.35
ax1.bar(x_positions - width/2, intact_forces, width, label='Intact', 
        color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
ax1.bar(x_positions + width/2, precut_forces, width, label='Precut',
        color=colors, alpha=0.5, edgecolor='black', linewidth=0.5, hatch='//')

ax1.set_ylabel('Peak Force $F_{\mathrm{max}}$ (kN)', fontsize=15, fontweight='bold')
ax1.set_title('Comparison of Peak Force', fontsize=15, fontweight='bold')
ax1.set_xticks(x_positions)
ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=13)
ax1.legend(fontsize=13, loc='upper left')
ax1.grid(axis='y', alpha=0.3)

                       
for boundary in family_boundaries[:-1]:
    ax1.axvline(boundary, color='gray', linestyle='--', linewidth=1, alpha=0.5)

                   
for pos, label in zip(family_positions, family_labels_list):
    ax1.text(pos, ax1.get_ylim()[1] * 0.95, label, 
             ha='center', fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

                      
ax2.bar(x_positions - width/2, intact_energies, width, label='Intact',
        color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
ax2.bar(x_positions + width/2, precut_energies, width, label='Precut',
        color=colors, alpha=0.5, edgecolor='black', linewidth=0.5, hatch='//')

ax2.set_ylabel('Work to Peak $W_{\mathrm{peak}}$ (J)', fontsize=15, fontweight='bold')
ax2.set_title('Comparison of Energy Absorption to Peak Load', fontsize=15, fontweight='bold')
ax2.set_xticks(x_positions)
ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=13)
ax2.legend(fontsize=13, loc='upper left')
ax2.grid(axis='y', alpha=0.3)

                       
for boundary in family_boundaries[:-1]:
    ax2.axvline(boundary, color='gray', linestyle='--', linewidth=1, alpha=0.5)

                   
for pos, label in zip(family_positions, family_labels_list):
    ax2.text(pos, ax2.get_ylim()[1] * 0.95, label,
             ha='center', fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('comparative_force_energy_all_families.png', dpi=150, bbox_inches='tight')
print("✓ Saved: comparative_force_energy_all_families.png")
plt.close()

print("\nKey observations:")
print(f"  Highest intact force: Symmetric 0.5 (F_max = 0.07244 kN)")
print(f"  Highest precut force: G2-B (F_max = 0.06439 kN)")
print(f"  Highest intact energy: G1-C (W_peak = 0.009249 J)")
print(f"  Highest precut energy: G2-B (W_peak = 0.003641 J)")
print(f"  Most damage-tolerant (force): G2-B (4.1% reduction)")
print(f"  Least damage-tolerant (force): G1-C (39.7% reduction)")
