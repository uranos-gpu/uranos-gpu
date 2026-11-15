#!/usr/bin/env python3
# Paper figure fig:tgv (Fig. 10, Table 5) from the 168 runs of this folder:
# reads out/<tag>/TGV_MONITOR.txt for every case listed in campaign.json.
# (a,b) enstrophy against time on the Cartesian and wavy grids, insets on the peak;
# (c) deficit of the peak enstrophy against the DNS; (d) change of the mean total energy,
# scaled by its initial value, against N_g.  Color family: viscous scheme; shade: resolution;
# line style and marker: grid.  The four KEP fluxes and the two gas models are drawn alike.
import os, re, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = HERE
plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'lines.linewidth': 0.8, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})

VISC = [('laplacian', 'Laplacian', '#374a95'), ('staggered', 'staggered', '#ab3013'),
        ('conservative', 'conservative', '#3db81d')]
RES = ['N0064', 'N0128', 'N0256', 'N0512']
NG = {r: int(r[1:]) for r in RES}
MIX = {'N0064': 0.60, 'N0128': 0.40, 'N0256': 0.20, 'N0512': 0.0}     # light (coarse) to dark (fine)
GRID = {'cart': ('-', 'o', 'Cartesian'), 'wavy': ('--', 's', 'wavy')}
ZV = {'laplacian': 1, 'staggered': 2, 'conservative': 3}

def shade(c, mix):
    r = np.array(to_rgb(c)); return tuple(r + (1 - r) * mix)

def load(f):
    U0 = float(re.search(r'U0\s*=\s*([0-9.Ee+-]+)', open(f).read()).group(1))
    d = np.loadtxt(f)
    return d[:, 2], d[:, 5] / U0**2, d[:, 10]          # t*, enstrophy, <rho E>

runs = {}
for c in json.load(open(os.path.join(DATA, 'campaign.json')))['cases']:
    f = os.path.join(DATA, 'out', c['tag'], 'TGV_MONITOR.txt')
    if os.path.isfile(f):
        runs[(c['grid'], c['gas'], c['conv'], c['visc'], f"N{c['N']:04d}")] = load(f)
print(len(runs), 'runs')
dns = np.loadtxt(os.path.join(DATA, 'reference', 'dns_vanrees2011_Re1600_512.txt'))
Epk_dns = dns[:, 3].max()

# conventions of the viscous verification figures: color family = viscous scheme
# (Laplacian blue, staggered red, conservative green), shade = resolution (light coarse,
# dark fine), solid with circles = Cartesian, dashed with squares = wavy.  The KEP fluxes
# overlap; markers (on the KEEP runs, staggered per scheme) reveal the overlapping schemes.
VCOL = {v: c for v, _, c in VISC}
OFFV = {'laplacian': 0.0, 'staggered': 1/3, 'conservative': 2/3}
fig = plt.figure(figsize=(7.0, 2.9))
gs = fig.add_gridspec(1, 2, wspace=0.22)
GAS = {'CP': 'calorically perfect gas', 'TP': 'thermally perfect gas'}
for col, gasp in enumerate(['CP', 'TP']):
    ax = fig.add_subplot(gs[0, col])
    a = ax.inset_axes([0.645, 0.53, 0.235, 0.41])
    for (grid, gas, conv, visc, res), (t, E, _) in sorted(runs.items(), key=lambda kv: (RES.index(kv[0][4]), ZV[kv[0][3]])):
        if gas != gasp: continue
        c = shade(VCOL[visc], MIX[res]); ls, mk, _ = GRID[grid]
        z = 2 + 3 * RES.index(res) + ZV[visc]
        kw = {}
        if conv == 'KEEP':
            step = 5.0; t0 = 1.5 + step * OFFV[visc] + (0.8 if grid == 'wavy' else 0.0)
            idx = [int(np.argmin(abs(t - tm))) for tm in np.arange(t0, 20, step)]
            kw = dict(marker=mk, ms=2.4, mfc='white', mew=0.5, markevery=idx)
        ax.plot(t, E, color=c, ls=ls, lw=0.7, zorder=z, **kw)
        a.plot(t, E, color=c, ls=ls, lw=0.7, zorder=z)
    ax.plot(dns[::100, 0], dns[::100, 3], ls='none', marker='o', ms=3.0, mfc='none', mec='k', mew=0.6, zorder=40)
    a.plot(dns[::8, 0], dns[::8, 3], ls='none', marker='o', ms=2.6, mfc='none', mec='k', mew=0.5, zorder=40)
    ax.set_xlim(0, 20); ax.set_ylim(0, 12.5)
    # resolution trend on the decaying branch: from the coarse (light) to the fine (dark) curves
    ax.annotate('', xy=(17.0, 5.0), xytext=(11.2, 3.6),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=0.6, mutation_scale=7), zorder=50)
    ax.text(11.0, 3.6, 'finer', ha='right', va='center', fontsize=7, zorder=50)
    ax.xaxis.set_major_locator(plt.MultipleLocator(5)); ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(2)); ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.tick_params(which='minor', length=1.5, width=0.5)
    ax.set_xlabel(r'$t\,u_\infty/\delta_0$')
    if col == 0: ax.set_ylabel(r'$\mathcal{E}\,\delta_0^2/u_\infty^2$')
    a.set_xlim(8.4, 9.6); a.set_ylim(9.65, 10.4)
    a.xaxis.set_major_locator(plt.MultipleLocator(0.5)); a.yaxis.set_major_locator(plt.MultipleLocator(0.25))
    a.tick_params(labelsize=6, length=2, width=0.5); a.yaxis.tick_right()
    for s_ in a.spines.values(): s_.set_linewidth(0.5)
    mark_inset(ax, a, loc1=2, loc2=3, fc='none', ec='0.55', lw=0.4)
    ax.text(0.03, 0.95, f"({'ab'[col]}) {GAS[gasp]}", transform=ax.transAxes, va='top')

# legend: one row, grouped (viscous scheme | resolution | grid | reference)
hv = [plt.Line2D([], [], color=c, lw=1.2, label=l) for _, l, c in VISC]
hr = [plt.Line2D([], [], color=shade('0.15', MIX[r]), lw=1.2, label=rf'${NG[r]}^3$') for r in RES]
hg = [plt.Line2D([], [], color='0.25', ls=GRID[g][0], marker=GRID[g][1], ms=2.6, mfc='white', mew=0.5,
                 label=GRID[g][2]) for g in ('cart', 'wavy')]
hd = [plt.Line2D([], [], ls='none', marker='o', ms=3.0, mfc='none', mec='k', mew=0.6, label='reference spectral DNS')]
blank = plt.Line2D([], [], color='none', label=' ')
# four columns: viscous scheme | resolution | grid | reference (column-major, 4 rows)
cols = [hv, hr, hg, hd]
h = [x for c in cols for x in (c + [blank] * (4 - len(c)))]
leg = fig.legend(handles=h, loc='upper center', ncol=4, frameon=True, fancybox=False, edgecolor='k',
                 bbox_to_anchor=(0.5, 0.0), handlelength=2.2, handletextpad=0.5, columnspacing=3.0,
                 borderpad=0.6, labelspacing=0.35)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'tgv_validation.pdf'); fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
