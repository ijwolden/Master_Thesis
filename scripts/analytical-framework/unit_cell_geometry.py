import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

                                                                                
plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'serif',
    'font.size': 10,
})

                
STRUT_COL = '#1c2b3a'                                   
DIM_COL   = '#606060'                                        
BLUE      = '#2563eb'            
RED       = '#b91c1c'            

                                                                     
LW_BASE  = 5.5                                         
LW_INC   = 2.5                                         


def draw_unit_cell(ax, a, x=1.0, z=1.0):

                      
    BL = np.array([0.0, 0.0])                
    BR = np.array([x,   0.0])                 
    AP = np.array([a*x, z  ])         

    theta_deg = np.degrees(np.arctan2(z, a * x))
    alpha_deg = np.degrees(np.arctan2(z, (1 - a) * x))

                                                                             
    ax.plot([BL[0], BR[0]], [BL[1], BR[1]],
            color=STRUT_COL, lw=LW_BASE, solid_capstyle='round', zorder=3)
    ax.plot([BL[0], AP[0]], [BL[1], AP[1]],
            color=STRUT_COL, lw=LW_INC,  solid_capstyle='round', zorder=3)
    ax.plot([BR[0], AP[0]], [BR[1], AP[1]],
            color=STRUT_COL, lw=LW_INC,  solid_capstyle='round', zorder=3)

                                                                             
    for p in (BL, BR, AP):
        ax.plot(*p, 'o', ms=6, color=STRUT_COL,
                markerfacecolor='white', markeredgewidth=1.5, zorder=5)

                                                                             
    r = 0.20               

    ax.add_patch(Arc(BL, 2*r, 2*r, angle=0,
                     theta1=0, theta2=theta_deg,
                     color=BLUE, lw=1.8, zorder=4))
    ax.add_patch(Arc(BR, 2*r, 2*r, angle=0,
                     theta1=180 - alpha_deg, theta2=180,
                     color=RED,  lw=1.8, zorder=4))

                                 
    tr = np.radians(theta_deg / 2)
    ar = np.radians(alpha_deg / 2)
    ax.text(BL[0] + 1.65*r*np.cos(tr),
            BL[1] + 1.65*r*np.sin(tr),
            r'$\theta$', ha='center', va='center', fontsize=13, color=BLUE)
    ax.text(BR[0] - 1.65*r*np.cos(ar),
            BR[1] + 1.65*r*np.sin(ar),
            r'$\alpha$', ha='center', va='center', fontsize=13, color=RED)

                                                                              
    arw = dict(arrowstyle='<->', color=DIM_COL, lw=1.0, mutation_scale=8)
    pad = 0.17

                                     
    yd = -pad
    for px in (BL[0], BR[0]):
        ax.plot([px, px], [0, yd], '--', color=DIM_COL, lw=0.7, alpha=0.45)
    ax.annotate('', xy=(BR[0], yd), xytext=(BL[0], yd), arrowprops=arw)
    ax.text(x / 2, yd - 0.07, r'$x$',
            ha='center', va='top', fontsize=12, color=DIM_COL)

                                  
    xd = -pad
    for py in (BL[1], AP[1]):
        ax.plot([0, xd], [py, py], '--', color=DIM_COL, lw=0.7, alpha=0.45)
    ax.annotate('', xy=(xd, AP[1]), xytext=(xd, BL[1]), arrowprops=arw)
    ax.text(xd - 0.07, z / 2, r'$z$',
            ha='right', va='center', fontsize=12, color=DIM_COL)

                                   
    yt = z + pad
    ax.plot([AP[0], AP[0]], [z, yt + 0.03], '--', color=DIM_COL, lw=0.7, alpha=0.45)
    for px in (BL[0], BR[0]):
        ax.plot([px, px], [0, yt], '--', color=DIM_COL, lw=0.7, alpha=0.45)
    ax.annotate('', xy=(AP[0], yt), xytext=(BL[0], yt), arrowprops=arw)
    ax.annotate('', xy=(BR[0], yt), xytext=(AP[0], yt), arrowprops=arw)
    ax.text(a * x / 2,       yt + 0.07, r'$ax$',
            ha='center', va='bottom', fontsize=11, color=DIM_COL)
    ax.text((a * x + x) / 2, yt + 0.07, r'$(1-a)x$',
            ha='center', va='bottom', fontsize=11, color=DIM_COL)

                                                                            
    ax.plot([AP[0], AP[0]], [0, z], ':', color='#aaaaaa', lw=0.9, zorder=1)

                                                                              
    def strut_label(p1, p2, text, side):
        mid = (p1 + p2) / 2.0
        d   = p2 - p1
        n   = np.array([-d[1], d[0]])
        n   = n / np.linalg.norm(n) * 0.13 * side
        ax.text(mid[0] + n[0], mid[1] + n[1], text,
                ha='center', va='center', fontsize=10, color=STRUT_COL,
                bbox=dict(boxstyle='round,pad=0.12',
                          fc='white', ec='none', alpha=0.9))

                                                                              
    strut_label(BL, BR, r'$2t$', side=+1)
                                                                                
    strut_label(BL, AP, r'$t$',  side=+1)
                                                                                  
    strut_label(BR, AP, r'$t$',  side=-1)

                                                                             
    ax.set_xlim(-0.50, 1.52)
    ax.set_ylim(-0.48, 1.58)
    ax.set_aspect('equal')
    ax.axis('off')


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.6))

draw_unit_cell(ax1, a=0.35)
ax1.set_title(r'Asymmetric:\ $a = 0.35$,\ $\theta \neq \alpha$',
              fontsize=11, pad=16)

draw_unit_cell(ax2, a=0.50)
ax2.set_title(r'Symmetric:\ $a = 0.50$,\ $\theta = \alpha \approx 63.4^\circ$',
              fontsize=11, pad=16)

plt.tight_layout(pad=1.5)

                                                                               
out_dir = os.path.dirname(os.path.abspath(__file__))
for ext in ('pdf', 'png'):
    path = os.path.join(out_dir, f'unit_cell_geometry.{ext}')
    plt.savefig(path, bbox_inches='tight', dpi=150)
    print(f'Saved  {path}')

plt.show()
