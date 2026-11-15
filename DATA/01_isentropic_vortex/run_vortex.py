#!/usr/bin/env python3
# Isentropic vortex runs for the URANOS-3.0 paper: all central (KEP) schemes,
# real Cartesian grid and wavy grid. Everything lives in DATA/01_isentropic_vortex.
import os, sys, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import vortex_lib as rv                       # mesh builder, deck template and constants

N, NZ = rv.N_INPLANE, rv.NZ
CASES = [('JP', 0), ('KEEP', 0), ('AEC', 0), ('PEP', 1)]    # (ep_type, thermally perfect)
MESH = os.path.join(HERE, 'mesh_wavy.txt')
REL = os.path.relpath(HERE, os.path.join(REPO, 'DATA'))

NCROSS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
SUF = '' if NCROSS == 1 else f'_x{NCROSS}'
ITOUT = rv.ITOUT if NCROSS == 1 else 100

def run(case):
    ep, tp, grid = case
    name = f'{ep.lower()}_{grid}{SUF}'
    dat = rv.DAT_TEMPLATE.format(
        xmin=rv.XMIN, xmax=rv.XMAX, ymin=rv.YMIN, ymax=rv.YMAX, zmin=rv.ZMIN, zmax=rv.ZMAX,
        mesh_file=(os.path.relpath(MESH, REPO) if grid == 'wavy' else ''),
        tmax=f'{NCROSS*rv.TMAX:.6f}', itmax=rv.ITMAX, nx=N, ny=N, nz=NZ, cfl=rv.CFL, mach=rv.MACH,
        itout=ITOUT, dir='x', ofn=name, data_dir=f'{REL}/{name}')
    dat = dat.replace("ep_type     = 'JP'", f"ep_type     = '{ep}'")
    dat = dat.replace('&end', f'\tthermally_perfect_gas = {tp}\n\tTref = 298.15\n&end')
    p = os.path.join(HERE, f'{name}.dat'); open(p, 'w').write(dat)
    shutil.rmtree(os.path.join(HERE, name), ignore_errors=True)
    os.makedirs(os.path.join(HERE, name, 'BINARY'), exist_ok=True)
    r = subprocess.run(['mpirun', '--oversubscribe', '-np', '1', rv.EXE, p], cwd=REPO,
                       stdout=open(p.replace('.dat', '.log'), 'w'), stderr=subprocess.STDOUT)
    return name, r.returncode

if __name__ == '__main__':
    XX, YY = rv.build_wavy(N, N); rv.write_mesh_file(MESH, XX, YY, N, N, NZ)
    cases = [(ep, tp, g) for g in ('cart', 'wavy') for ep, tp in CASES]
    with ThreadPoolExecutor(max_workers=8) as ex:
        for name, rc in ex.map(run, cases): print(name, 'exit', rc, flush=True)
