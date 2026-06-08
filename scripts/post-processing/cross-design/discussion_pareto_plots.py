import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
from pathlib import Path

BASE = Path(__file__).parent

                                                                               
def parse_sections(path):
    datasets, label, rows = {}, None, []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('#'):
                if label and rows:
                    datasets[label] = np.array(rows, float)
                label = line[1:].strip().lower().replace(' ', '_')
                rows = []
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        rows.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        pass
    if label and rows:
        datasets[label] = np.array(rows, float)
    return datasets


def stiffness(d, f, frac=0.1):
    i_max = np.argmax(f)
    mask = d <= frac * d[i_max]
    if mask.sum() < 3:
        return np.nan
    return np.polyfit(d[mask], f[mask], 1)[0]


def peak_metrics(d, f):
    k = stiffness(d, f)
    return k, float(np.max(f))


def eta_ratio(F_intact, F_precut):
    return F_intact / F_precut if F_precut > 0 else np.nan


records = []


bl = parse_sections(BASE / 'Baseline/Baseline.txt')
for tag, pre_tag, name in [
    ('a',  'a-precut',  'Baseline A'),
    ('b',  'b_precut',  'Baseline B'),
]:
    if tag in bl and pre_tag in bl:
        ki, Fi = peak_metrics(bl[tag][:, 0], bl[tag][:, 1])
        _,  Fp = peak_metrics(bl[pre_tag][:, 0], bl[pre_tag][:, 1])
        records.append(dict(family='Baseline', label=name, k=ki, F_max=Fi,
                            eta=eta_ratio(Fi, Fp)))


sym_int = pd.read_csv(BASE / 'Symmetric/f-d-all.txt',    sep='\t', header=0)
sym_pre = pd.read_csv(BASE / 'Symmetric/f-d-precut.txt', sep='\t', header=0)

designs = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']
for i, des in enumerate(designs):
    ci = i * 3

    def _col(df, c):
        v = pd.to_numeric(df.iloc[:, c], errors='coerce').dropna().values
        return v

    d_i, f_i = _col(sym_int, ci), _col(sym_int, ci + 1)
    d_p, f_p = _col(sym_pre, ci), _col(sym_pre, ci + 1)
    n = min(len(d_i), len(f_i)); d_i, f_i = d_i[:n], f_i[:n]
    n = min(len(d_p), len(f_p)); d_p, f_p = d_p[:n], f_p[:n]

    if len(d_i) > 3 and len(d_p) > 3:
        ki, Fi = peak_metrics(d_i, f_i)
        _,  Fp = peak_metrics(d_p, f_p)
        records.append(dict(family='Symmetric', label=f'Sym {des}', k=ki,
                            F_max=Fi, eta=eta_ratio(Fi, Fp)))


def load_g1_rate1(path):
    df = pd.read_csv(path, sep='\t', header=0)
    d = pd.to_numeric(df.iloc[:, 6], errors='coerce').dropna().values
    f = pd.to_numeric(df.iloc[:, 7], errors='coerce')
    f = f.dropna().values
    n = min(len(d), len(f))
    return d[:n], f[:n]


def load_g1_precut(path):
    df = pd.read_csv(path, sep='\t', skiprows=[1], header=0)                       
    d = pd.to_numeric(df.iloc[:, 0], errors='coerce').dropna().values
    f = pd.to_numeric(df.iloc[:, 1], errors='coerce')
    f = f.dropna().values
    n = min(len(d), len(f))
    return d[:n], f[:n]


for v, name in [('a', 'G1-A'), ('b', 'G1-B'), ('c', 'G1-C')]:
    try:
        d_i, f_i = load_g1_rate1(BASE / f'G1/g1-{v}.txt')
        d_p, f_p = load_g1_precut(BASE / f'G1/{v}-precut.txt')
        if len(d_i) > 3 and len(d_p) > 3:
            ki, Fi = peak_metrics(d_i, f_i)
            _,  Fp = peak_metrics(d_p, f_p)
            records.append(dict(family='G1', label=name, k=ki, F_max=Fi,
                                eta=eta_ratio(Fi, Fp)))
    except Exception as exc:
        print(f"G1-{v}: {exc}")


for v in ['A', 'B']:
    secs = parse_sections(BASE / f'G2/G2-{v}.txt')
    intact = {k: v2 for k, v2 in secs.items() if 'precut' not in k}
    precut = {k: v2 for k, v2 in secs.items() if 'precut' in k}
    if intact and precut:
        di, fi = next(iter(intact.values()))[:, 0], next(iter(intact.values()))[:, 1]
        dp, fp = next(iter(precut.values()))[:, 0], next(iter(precut.values()))[:, 1]
        ki, Fi = peak_metrics(di, fi)
        _,  Fp = peak_metrics(dp, fp)
        records.append(dict(family='G2', label=f'G2-{v}', k=ki, F_max=Fi,
                            eta=eta_ratio(Fi, Fp)))


g3 = parse_sections(BASE / 'G3/g3-f-d.txt')
g3_int = {k: v for k, v in g3.items() if 'precut' not in k}
g3_pre = {k: v for k, v in g3.items() if 'precut' in k}
if g3_int and g3_pre:
    di, fi = next(iter(g3_int.values()))[:, 0], next(iter(g3_int.values()))[:, 1]
    dp, fp = next(iter(g3_pre.values()))[:, 0], next(iter(g3_pre.values()))[:, 1]
    ki, Fi = peak_metrics(di, fi)
    _,  Fp = peak_metrics(dp, fp)
    records.append(dict(family='G3', label='G3', k=ki, F_max=Fi,
                        eta=eta_ratio(Fi, Fp)))


df = pd.DataFrame(records).dropna()
print("\n=== Computed metrics ===")
print(df.to_string(index=False))

                                                                                
def is_pareto(df, col_y):
    pareto = []
    for _, row in df.iterrows():
        dominated = any(
            (other['eta'] <= row['eta']) and (other[col_y] >= row[col_y])
            and not (other['eta'] == row['eta'] and other[col_y] == row[col_y])
            for _, other in df.iterrows()
        )
        pareto.append(not dominated)
    return np.array(pareto)


FAMILY = {
    'Baseline':  dict(color='#4575b4', marker='D', zorder=4),
    'Symmetric': dict(color='#f46d43', marker='o', zorder=4),
    'G1':        dict(color='#74add1', marker='s', zorder=4),
    'G2':        dict(color='#1a9641', marker='^', zorder=5),
    'G3':        dict(color='#9970ab', marker='P', zorder=4),
}

                                                                   
ANNOTATE = {
    'Baseline A', 'Baseline B',
    'Sym 0.4', 'Sym 0.5', 'Sym 0.9',
    'G1-A', 'G1-B', 'G1-C',
    'G2-A', 'G2-B',
    'G3',
}

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 9,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for ax, col_y, ylabel, title in [
    (axes[0], 'F_max', 'Peak force  $F_{\\mathrm{max}}$ (kN)',
                       '(a) Intact strength vs. notch sensitivity'),
    (axes[1], 'k',     'Initial stiffness  $k$ (kN/mm)',
                       '(b) Stiffness vs. notch sensitivity'),
]:
    pareto_mask = is_pareto(df, col_y)

                             
    for family, grp in df.groupby('family'):
        st = FAMILY[family]
        ax.scatter(grp['eta'], grp[col_y],
                   c=st['color'], marker=st['marker'],
                   s=70, zorder=st['zorder'],
                   edgecolors='white', linewidths=0.6,
                   label=family)

                            
    pareto_pts = df[pareto_mask]
    ax.scatter(pareto_pts['eta'], pareto_pts[col_y],
               facecolors='none', edgecolors='black',
               linewidths=1.2, s=120, zorder=6, label='Pareto-optimal')

                                                                       
    ax.axvline(1.0, color='#666', lw=0.8, ls=':', alpha=0.8)
    ax.axvline(1.2, color='#999', lw=0.8, ls='--', alpha=0.6)

                       
    LABEL_OFFSET = {
        'G2-B':       ( 42,   4, 'left'),
        'Sym 0.5':    (  5,   4, 'left'),
        'Sym 0.4':    (  5,   4, 'left'),
        'Sym 0.9':    ( -6,   4, 'right'),
        'Baseline A': (  5, -11, 'left'),
        'Baseline B': ( -6,   4, 'right'),
        'G1-A':       (  5,   4, 'left'),
        'G1-B':       (  5,   4, 'left'),
        'G1-C':       (  5,   4, 'left'),
        'G2-A':       (  5,  -11, 'left'),
        'G3':         (  5,   4, 'left'),
    }
    for _, row in df.iterrows():
        if row['label'] not in ANNOTATE:
            continue
        color = FAMILY[row['family']]['color']
        xo, yo, ha = LABEL_OFFSET.get(row['label'], (5, 4, 'left'))
        ax.annotate(row['label'],
                    xy=(row['eta'], row[col_y]),
                    xytext=(xo, yo),
                    textcoords='offset points',
                    fontsize=7.2, color=color, ha=ha,
                    arrowprops=dict(arrowstyle='-', color=color,
                                   lw=0.4, alpha=0.6))

                                                              
    ylo, yhi = ax.get_ylim()
    y_text = ylo + 0.96 * (yhi - ylo)
    ax.text(1.005, y_text, r'$\eta=1.0$', fontsize=7.5,
            color='#555', va='top', ha='left')
    ax.text(1.205, y_text, r'$\eta=1.2$', fontsize=7.5,
            color='#888', va='top', ha='left')

    y_arrow = ylo + 0.15 * (yhi - ylo)
    ax.annotate('', xy=(1.02, y_arrow), xytext=(1.52, y_arrow),
                arrowprops=dict(arrowstyle='->', color='#bbb', lw=0.9))
    ax.text(1.27, y_arrow + 0.005 * (yhi - ylo), 'better tolerance',
            fontsize=7, color='#aaa', ha='center', va='bottom')

    ax.set_xlabel('Notch sensitivity  $\\eta = F_{\\mathrm{intact}} / F_{\\mathrm{precut}}$',
                  fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=6)
    ax.legend(fontsize=7.5, loc='upper left',
              framealpha=0.5, edgecolor='none', ncol=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.set_xlim(left=0.95)


fig.suptitle('Design performance trade-off: strength and stiffness vs. fracture tolerance',
             fontsize=11, y=1.01)
fig.tight_layout()

out = BASE / 'discussion_pareto_scatter.png'
fig.savefig(out, dpi=200, bbox_inches='tight')
print(f"\nSaved {out}")
plt.show()
