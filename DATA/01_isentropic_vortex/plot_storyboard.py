#!/usr/bin/env python3
# Paper figure: isentropic vortex storyboard on the wavy grid, sixth order JP, calorically perfect gas.
# Advection along x reuses the run jp_wavy (run_vortex.py); advection along y is run here as jp_wavy_y.
import os, sys, glob, subprocess
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import vortex_lib as rv
import make_gif as mg

plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})

N, NZ = rv.N_INPLANE, rv.NZ
MESH = os.path.join(HERE, 'mesh_wavy.txt')
REL = os.path.relpath(HERE, os.path.join(REPO, 'DATA'))

def run_y():
    name = 'jp_wavy_y'
    if glob.glob(os.path.join(HERE, name, 'BINARY', '*.bin')): return name
    dat = rv.DAT_TEMPLATE.format(
        xmin=rv.XMIN, xmax=rv.XMAX, ymin=rv.YMIN, ymax=rv.YMAX, zmin=rv.ZMIN, zmax=rv.ZMAX,
        mesh_file=os.path.relpath(MESH, REPO), tmax=f'{rv.TMAX:.6f}', itmax=rv.ITMAX,
        nx=N, ny=N, nz=NZ, cfl=rv.CFL, mach=rv.MACH, itout=rv.ITOUT, dir='y',
        ofn=name, data_dir=f'{REL}/{name}')
    p = os.path.join(HERE, f'{name}.dat'); open(p, 'w').write(dat)
    os.makedirs(os.path.join(HERE, name, 'BINARY'), exist_ok=True)
    subprocess.run(['mpirun', '--oversubscribe', '-np', '1', rv.EXE, p], cwd=REPO, check=True,
                   stdout=open(p.replace('.dat', '.log'), 'w'), stderr=subprocess.STDOUT)
    return name

def frames(name, metr):
    files = sorted(glob.glob(os.path.join(HERE, name, 'BINARY', '*.bin')))
    stamped = sorted((mg.read_binary(f)[0], f) for f in files)
    out = []
    for frac in (0.0, 1/3, 2/3, 1.0):                         # t u_inf/L_x = 0, 1/3, 2/3, 1
        t, f = min(stamped, key=lambda q: abs(q[0]/rv.TMAX - frac))
        phi = mg.read_binary(f)[4]
        out.append((frac, mg.omega_z_slice(phi, metr, phi.shape[3]//2)))
    return out

nx, ny, nz, xx, yy = mg.read_mesh(MESH)
metr = mg.mesh_metrics(xx, yy, nx, ny)
X = xx[mg.NG:mg.NG+nx, mg.NG:mg.NG+ny]; Y = yy[mg.NG:mg.NG+nx, mg.NG:mg.NG+ny]
rows = [frames('jp_wavy', metr), frames(run_y(), metr)]
vmax = float(np.abs(rows[0][0][1]).max())

fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.6), sharex=True, sharey=True,
                         gridspec_kw=dict(wspace=0.06, hspace=0.10))
labels = iter('abcdefgh')
tlab = {0.0: '0', 1/3: '1/3', 2/3: '2/3', 1.0: '1'}
for r, row in enumerate(rows):
    for c, (frac, w) in enumerate(row):
        ax = axes[r, c]
        pc = ax.pcolormesh(X, Y, w, cmap='RdBu_r', shading='gouraud', vmin=-vmax, vmax=vmax, rasterized=True)
        for i in range(0, nx, 3): ax.plot(X[i, :], Y[i, :], color='0.35', lw=0.2)
        for j in range(0, ny, 3): ax.plot(X[:, j], Y[:, j], color='0.35', lw=0.2)
        ax.set_aspect('equal')
        ax.text(0.03, 0.97, f'({next(labels)})', transform=ax.transAxes, va='top')
        ax.text(0.97, 0.97, rf'$t\,u_\infty/L_{"xy"[r]}={tlab[frac]}$', transform=ax.transAxes, va='top', ha='right',
                bbox=dict(fc='white', ec='none', pad=0.8, alpha=0.8))
        ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4])
        if r == 1: ax.set_xlabel(r'$x/\delta_0$')
        if c == 0: ax.set_ylabel(r'$y/\delta_0$')
cb = fig.colorbar(pc, ax=axes, fraction=0.025, pad=0.02)
cb.set_label(r'$\omega_z\,\delta_0/u_0$'); cb.ax.tick_params(width=0.6)
out = os.path.join(HERE, 'vortex_storyboard.pdf'); fig.savefig(out, dpi=300); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
print('vmax', vmax)
