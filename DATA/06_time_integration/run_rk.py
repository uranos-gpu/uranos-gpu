#!/usr/bin/env python3
# Temporal order of accuracy of the three stage Runge--Kutta scheme inside URANOS, for the paper.
# Self convergence in time: grid and spatial operators fixed, only Dt refined, error measured
# against a reference run with Dt_ref = Dt0/REF on the same grid (the spatial error cancels).
#   dwave : harmonic density wave, 2D periodic, inviscid; Cartesian and wavy grid; KEEP (CP), PEP (TP)
#   tgv   : Taylor--Green vortex 32^3, viscous, nonlinear; KEEP (CP)
#   stab  : density wave on the Cartesian grid at increasing CFL, to find the stability limit
import os, sys, shutil, subprocess, json
import numpy as np
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import rk_lib as sw          # deck template COMMON, read_bin
import kolmo_lib as rs              # wavy mesh builder (domain [-1,1]^2, a_w = 0.1)
EXE = os.path.join(REPO, 'Uranos.exe')
REL = os.path.relpath(HERE, os.path.join(REPO, 'DATA'))

DW = dict(ic='periodic_euler_x', dims=2, nx=32, ny=32, nz=0, xmin=-1.0, xmax=1.0, ymin=-1.0, ymax=1.0,
          zmin=0.0, zmax=0.0, Re=1000.0, Ma=0.5, viscous='.false.')   # rho = 1 + 0.2 sin(pi x), u = 0.5
TG = dict(ic='TGV', dims=3, nx=32, ny=32, nz=32, xmin=0.0, xmax=2*np.pi, ymin=0.0, ymax=2*np.pi,
          zmin=0.0, zmax=2*np.pi, Re=160.0, Ma=0.2, viscous='.true.')
PER = {f'bc{i}': 'periodic' for i in range(1, 7)}
DT0 = {'dwave': 0.024, 'tgv': 4.0e-3}
N0 = {'dwave': 40, 'tgv': 100}
REFS = [1, 2, 4, 8, 16]
REF = 128

def mesh_wavy():
    m = os.path.join(HERE, 'mesh_wavy_32.txt')
    if not os.path.isfile(m):
        XX, YY = rs.build_mesh('wavy', 32, 32); rs.write_mesh_file(m, XX, YY, 32, 32, 1)
    return os.path.relpath(m, REPO)

def deck(name, cfg, dt, itmax, ep, tp, mesh=''):
    d = sw.COMMON.format(dt=dt, itmax=itmax, data_dir=f'{REL}/{name}', **PER, **cfg)
    d = d.replace("ep_type    = 'KEEP'", f"ep_type    = '{ep}'")
    d = d.replace('thermally_perfect_gas = 0', f'thermally_perfect_gas = {tp}')
    if cfg['ic'] == 'TGV': d = d.replace('vis_flag = 0', 'vis_flag = 2')
    if mesh: d = d.replace("gridpoint_z = 'uniform'", f"gridpoint_z = 'uniform'\n    mesh_file = '{mesh}'")
    return d

def run(job):
    name, cfg, dt, itmax, ep, tp, mesh, np_ = job
    out = os.path.join(HERE, name)
    if os.path.isfile(os.path.join(out, 'done')): return name, 0
    p = os.path.join(HERE, f'{name}.dat'); open(p, 'w').write(deck(name, cfg, dt, itmax, ep, tp, mesh))
    shutil.rmtree(out, ignore_errors=True); os.makedirs(os.path.join(out, 'BINARY'))
    r = subprocess.run(['mpirun', '--oversubscribe', '-np', str(np_), EXE, p], cwd=REPO,
                       stdout=open(p.replace('.dat', '.log'), 'w'), stderr=subprocess.STDOUT)
    if r.returncode == 0: open(os.path.join(out, 'done'), 'w').close()
    return name, r.returncode

if __name__ == '__main__':
    jobs = []
    mw = mesh_wavy()
    for grid in ('cart', 'wavy'):
        for ep, tp in (('KEEP', 0), ('PEP', 1)):
            for r in REFS + [REF]:
                jobs.append((f'dwave_{grid}_{ep}_r{r:03d}', DW, DT0['dwave'] / r, N0['dwave'] * r, ep, tp,
                             mw if grid == 'wavy' else '', 1))
    for r in [1, 2, 4, 8, 64]:
        jobs.append((f'tgv_cart_KEEP_r{r:03d}', TG, DT0['tgv'] / r, N0['tgv'] * r, 'KEEP', 0, '', 4))
    # stability scan at fixed physical time T = 20 (five crossings), spectral radius of the
    # Cartesian 32^2 grid summed over the two directions, largest sound speed (rho = 0.8)
    c = np.sqrt(1.4 / 0.8); u = 0.5; dx = 2.0 / 32
    S = (u + c) / dx + c / dx
    for cfl in np.arange(0.2, 1.61, 0.1):
        dt = cfl / S; n = int(np.ceil(20.0 / dt)); dt = 20.0 / n
        for ep, tp in (('KEEP', 0), ('PEP', 1)):
            jobs.append((f'stab_{ep}_cfl{cfl:.2f}', DW, dt, n, ep, tp, '', 1))
    json.dump(dict(S=S, dt0=DT0, n0=N0, refs=REFS, ref=REF), open(os.path.join(HERE, 'setup.json'), 'w'))
    with ThreadPoolExecutor(max_workers=6) as ex:
        for name, rc in ex.map(run, jobs): print(name, 'rc', rc, flush=True)
