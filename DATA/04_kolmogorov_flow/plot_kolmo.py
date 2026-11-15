#!/usr/bin/env python3
# Paper figure: steady Kolmogorov flow, the three viscous schemes against the exact solution.
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import kolmo_lib as rs
from run_kolmo import N, SCHEMES, MESHES

plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'lines.linewidth': 0.8,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})

# sixth order shades of the viscous error scaling figure; line style and marker = grid
CASES = [('laplacian', 'Laplacian', '#374a95'), ('staggered', 'staggered', '#ab3013'),
         ('conservative', 'conservative', '#3db81d')]
GRID = {'real_cart': ('-', 'o', 'Cartesian grid'), 'wavy': ('--', 's', 'wavy grid')}
OFF = {'laplacian': 0, 'staggered': 1, 'conservative': 2}
U = rs.A_AN                                    # amplitude of the exact solution, in units of u_0

def mk(key, mesh, start, step):
    return dict(marker=GRID[mesh][1], ms=2.8, mfc='white', mew=0.6,
                markevery=(start + OFF[key]*step/3, step))

def line(diff, mesh):
    nx, ny, nz = rs.case_dims(mesh, N)
    XX, YY = rs.build_mesh(mesh, nx, ny)
    fn = rs.final_dump(os.path.join(HERE, f'{diff}_{mesh}', 'BINARY'))
    t, _, _, _, u = rs.read_binary_u(fn)
    i = nx//2
    y = YY[rs.NG + i, rs.NG:rs.NG + ny]            # grid line of constant xi
    uu = u[i, :, nz//2]
    yy = YY[rs.NG:rs.NG+nx, rs.NG:rs.NG+ny][:, :, None] + 0*u
    detJ = rs.detJ2D_interior(XX, YY, nx, ny)[:, :, None] + 0*u
    L2 = np.sqrt(np.sum((u - U*np.sin(rs.KWAVE*yy))**2*detJ)/np.sum(detJ))/U
    return y, uu/U, L2

data = {(d, m): line(d, m) for m in MESHES for d, *_ in CASES}
for k, v in data.items(): print(k, 'L2/U = %.3e' % v[2])

fig = plt.figure(figsize=(7.0, 4.4))
gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.24)
labels = iter('abcd')
yf = np.linspace(rs.YMIN, rs.YMAX, 2001)
for col, mesh in enumerate(MESHES):
    ls, m, title = GRID[mesh]
    ax = fig.add_subplot(gs[0, col])
    ax.plot(yf, np.sin(rs.KWAVE*yf), color='k', lw=0.8)
    for key, lab, c in CASES:
        y, u, _ = data[(key, mesh)]
        ax.plot(y, u, color=c, ls=ls, **mk(key, mesh, 0.02, 0.18))
    ax.set_xlim(rs.YMIN, rs.YMAX); ax.set_ylim(-1.2, 1.45)
    ax.set_xlabel(r'$y/\delta_0$'); ax.set_ylabel(r'$u/U_4$')
    ax.text(0.03, 0.90, f'({next(labels)}) {title}', transform=ax.transAxes)
for col, mesh in enumerate(MESHES):
    ls, m, title = GRID[mesh]
    ax = fig.add_subplot(gs[1, col])
    for key, lab, c in CASES:
        y, u, _ = data[(key, mesh)]
        ax.semilogy(y, np.abs(u - np.sin(rs.KWAVE*y)), color=c, ls=ls, **mk(key, mesh, 0.02, 0.18))
    ax.set_xlim(rs.YMIN, rs.YMAX); ax.set_ylim(1e-8, 1e-2)
    ax.set_xlabel(r'$y/\delta_0$'); ax.set_ylabel(r'$|u-u_{\mathrm{exact}}|/U_4$')
    ax.text(0.03, 0.90, f'({next(labels)}) {title}', transform=ax.transAxes)
handles = [plt.Line2D([], [], color='k', lw=0.8, label='Exact')] + \
          [plt.Line2D([], [], color=c, lw=0.8, label=l) for _, l, c in CASES] + \
          [plt.Line2D([], [], color='0.3', ls=GRID[g][0], marker=GRID[g][1], ms=2.8, mfc='white', mew=0.6,
                      label=GRID[g][2].split()[0]) for g in MESHES]
leg = fig.legend(handles=handles, loc='upper center', ncol=6, frameon=True, fancybox=True, edgecolor='k',
                 bbox_to_anchor=(0.5, 0.035), handlelength=2.8, columnspacing=1.6)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'kolmogorov_profiles.pdf'); fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
