#!/usr/bin/env python3
# Paper figure: temporal order of the three stage Runge--Kutta scheme (self convergence in time).
import os, sys, json, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import rk_lib as sw
import kolmo_lib as rs
plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'lines.linewidth': 0.8, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
S = json.load(open(os.path.join(HERE, 'setup.json')))

def last(name):
    f = sorted(glob.glob(os.path.join(HERE, name, 'BINARY', '*.bin')))[-1]
    return sw.read_bin(f)

def weights(grid, nx, ny, nz):
    if grid == 'wavy':
        XX, YY = rs.build_mesh('wavy', nx, ny); w = rs.detJ2D_interior(XX, YY, nx, ny)[:, :, None]
        return np.broadcast_to(w, (nx, ny, nz))
    return np.ones((nx, ny, nz))

def spectral_radius(grid):
    # sum over the in plane directions of (|u.grad xi_j| + c |grad xi_j|), xi_j in index units
    nx = ny = 32; c = np.sqrt(1.4 / 0.8); u = 0.5   # largest sound speed (rho = 0.8), wave speed
    XX, YY = rs.build_mesh(grid if grid == 'wavy' else 'real_cart', nx, ny)
    g = rs.NG; sl = (slice(g, g + nx), slice(g, g + ny))
    xe, xn = rs.d_idx(XX, 0)[sl], rs.d_idx(XX, 1)[sl]; ye, yn = rs.d_idx(YY, 0)[sl], rs.d_idx(YY, 1)[sl]
    J = xe * yn - xn * ye
    gxi = np.stack([yn / J, -xn / J]); geta = np.stack([-ye / J, xe / J])
    return float(np.max(abs(u * gxi[0]) + c * np.hypot(*gxi) + abs(u * geta[0]) + c * np.hypot(*geta)))

def sweep(prefix, grid):
    ref = last(f'{prefix}_r{S["ref"]:03d}' if 'dwave' in prefix else f'{prefix}_r064')
    nx, ny, nz = ref['nx'], ref['ny'], ref['nz']; w = weights(grid, nx, ny, nz)
    refs = S['refs'] if 'dwave' in prefix else [1, 2, 4, 8]
    key = 'dwave' if 'dwave' in prefix else 'tgv'
    dt, err = [], []
    for r in refs:
        d = last(f'{prefix}_r{r:03d}')
        e = d['phi'][..., 0] - ref['phi'][..., 0]
        dt.append(S['dt0'][key] / r); err.append(np.sqrt(np.sum(e**2 * w) / np.sum(w)))
    return np.array(dt), np.array(err)

COL = {'KEEP': '#374a95', 'PEP': '#ab3013'}; LS = {'KEEP': '-', 'PEP': '--'}
MK = {'cart': 'o', 'wavy': 's'}
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw=dict(wspace=0.3))
rows = []
for grid in ('cart', 'wavy'):
    Sr = spectral_radius(grid)
    for ep in ('KEEP', 'PEP'):
        dt, err = sweep(f'dwave_{grid}_{ep}', grid)
        cfl = dt * Sr
        ax.loglog(cfl, err, color=COL[ep], ls=LS[ep], marker=MK[grid], ms=3.2, mfc='white', mew=0.6)
        p = np.log(err[:-1] / err[1:]) / np.log(dt[:-1] / dt[1:]); rows.append((f'dwave {grid} {ep}', cfl, err, p))
x = np.array([0.06, 1.8]); ax.loglog(x, 0.15 * 4.76e-7 / 1.21**3 * x**3, 'k:', lw=0.6)
ax.set_xlabel(r'$\mathrm{CFL}$'); ax.set_ylabel(r'$\mathcal{L}_2$')
ax.text(0.03, 0.95, '(a) density wave', transform=ax.transAxes, va='top')
# TGV: spectral radius from the initial field, sum over the three directions
U0 = 1.1832159566199232; c = np.sqrt(1.4 * 100.0); dx = 2 * np.pi / 32
Stg = (2 * (U0 + c) + c) / dx
dt, err = sweep('tgv_cart_KEEP', 'cart'); cfl = dt * Stg
bx.loglog(cfl, err, color=COL['KEEP'], ls='-', marker='o', ms=3.2, mfc='white', mew=0.6)
p = np.log(err[:-1] / err[1:]) / np.log(dt[:-1] / dt[1:]); rows.append(('tgv cart KEEP', cfl, err, p))
x = np.array([0.08, 1.0]); bx.loglog(x, err[0] / cfl[0]**3 * 0.3 * x**3, 'k:', lw=0.6)
bx.set_xlabel(r'$\mathrm{CFL}$'); bx.set_ylabel(r'$\mathcal{L}_2$')
bx.text(0.03, 0.95, r'(b) Taylor--Green vortex, $32^3$', transform=bx.transAxes, va='top')
for r in rows:
    print(r[0], 'CFL', np.round(r[1], 3), 'err', ['%.2e' % e for e in r[2]], 'order', np.round(r[3], 2))
h = [plt.Line2D([], [], color=COL['KEEP'], ls='-', label=r'\texttt{KEEP}, calorically perfect'),
     plt.Line2D([], [], color=COL['PEP'], ls='--', label=r'\texttt{PEP}, thermally perfect'),
     plt.Line2D([], [], color='0.3', ls='none', marker='o', ms=3.2, mfc='white', mew=0.6, label='Cartesian'),
     plt.Line2D([], [], color='0.3', ls='none', marker='s', ms=3.2, mfc='white', mew=0.6, label='wavy'),
     plt.Line2D([], [], color='k', ls=':', lw=0.6, label=r'$\mathrm{CFL}^3$ reference')]
leg = fig.legend(handles=h, loc='upper center', ncol=5, frameon=True, fancybox=False, edgecolor='k',
                 bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.4)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'rk_order.pdf'); fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
# stability scan: final dump at T = 20, finite and with the initial amplitude 0.2
import re
for ep in ('KEEP', 'PEP'):
    for d in sorted(glob.glob(os.path.join(HERE, f'stab_{ep}_cfl*/'))):
        n = os.path.basename(d.rstrip('/')); cfl = float(n.split('cfl')[1])
        fs = sorted(glob.glob(os.path.join(d, 'BINARY', '*.bin')))
        if not fs: print(n, 'no output'); continue
        b = sw.read_bin(fs[-1]); r = b['phi'][..., 0]
        amp = float(np.nanmax(np.abs(r - 1))) if np.all(np.isfinite(r)) else float('nan')
        print(f"{ep} CFL={cfl:.2f}  t={b['time']:.3f}  max|rho-1|={amp:.4f}")
