import matplotlib.pyplot as plt
import numpy as np

def read_force_displacement(filename):
    regular_x, regular_y = [], []
    precut_x, precut_y = [], []
    reading_precut = False
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
                                                        
            if line.startswith('#Precut'):
                reading_precut = True
                continue
                              
            if not line:
                continue
                        
            try:
                values = line.split()
                if len(values) == 2:
                    x, y = float(values[0]), float(values[1])
                    if reading_precut:
                        precut_x.append(x)
                        precut_y.append(y)
                    else:
                        regular_x.append(x)
                        regular_y.append(y)
            except ValueError as e:
                print(f"Error parsing line '{line}': {e}")
                continue
    
    print(f"Read {filename}: regular={len(regular_x)}, precut={len(precut_x)}")
    return (np.array(regular_x), np.array(regular_y)), (np.array(precut_x), np.array(precut_y))

                           
g2a_regular, g2a_precut = read_force_displacement('G2-A.txt')
g2b_regular, g2b_precut = read_force_displacement('G2-B.txt')

                           
print(f"G2-A regular: {len(g2a_regular[0])} points")
print(f"G2-A precut: {len(g2a_precut[0])} points")
print(f"G2-B regular: {len(g2b_regular[0])} points")
print(f"G2-B precut: {len(g2b_precut[0])} points")

                 
plt.figure(figsize=(10, 6))

plt.plot(g2a_regular[0], g2a_regular[1], label='G2-A', linewidth=2, color='C0')
plt.plot(g2a_precut[0], g2a_precut[1], label='G2-A-precut', linewidth=2, linestyle='--', color='C0')
plt.plot(g2b_regular[0], g2b_regular[1], label='G2-B', linewidth=2, color='C1')
plt.plot(g2b_precut[0], g2b_precut[1], label='G2-B-precut', linewidth=2, linestyle='--', color='C1')

plt.xlabel('Displacement', fontsize=12)
plt.ylabel('Force', fontsize=12)
plt.title('G2', fontsize=14)
plt.legend(fontsize=10)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.grid(True, alpha=0.3)
plt.tight_layout()

               
plt.savefig('G2_force_displacement.png', dpi=300, bbox_inches='tight')
plt.show()

print("Plot saved as 'G2_force_displacement.png'")
