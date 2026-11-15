#!/usr/bin/env python3
# Steady Kolmogorov flow for the URANOS-3.0 paper: the three viscous schemes at sixth order,
# Cartesian and wavy grids, one resolution. Everything lives in DATA/04_kolmogorov_flow.
import os, sys, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import kolmo_lib as rs                      # deck template, mesh builder and constants

N, FD = 48, 6
SCHEMES = ['laplacian', 'staggered', 'conservative']
MESHES = ['real_cart', 'wavy']
REL = os.path.relpath(HERE, os.path.join(REPO, 'DATA'))

def run(case):
    diff, mesh = case
    name = f'{diff}_{mesh}'
    nx, ny, nz = rs.case_dims(mesh, N)
    mesh_file = ''
    if mesh != 'real_cart':
        XX, YY = rs.build_mesh(mesh, nx, ny)
        mesh_file = os.path.join(HERE, f'mesh_{mesh}.txt'); rs.write_mesh_file(mesh_file, XX, YY, nx, ny, nz)
        mesh_file = os.path.relpath(mesh_file, REPO)
    dat = rs.DAT_TEMPLATE.format(
        xmin=rs.XMIN, xmax=rs.XMAX, ymin=rs.YMIN, ymax=rs.YMAX, zmin=rs.ZMIN, zmax=rs.ZMAX,
        mesh_file=mesh_file, itmax=rs.itmax_for(mesh, N), nx=nx, ny=ny, nz=nz, cfl=rs.CFL,
        re=rs.RE, mach=rs.MACH, forcing=f'{rs.FORCING:.6E}', fd_order=FD, ofn=name,
        data_dir=f'{REL}/{name}', diffusion=diff)
    p = os.path.join(HERE, f'{name}.dat'); open(p, 'w').write(dat)
    shutil.rmtree(os.path.join(HERE, name), ignore_errors=True)
    os.makedirs(os.path.join(HERE, name, 'BINARY'), exist_ok=True)
    r = subprocess.run(['mpirun', '--oversubscribe', '-np', '1', rs.EXE, p], cwd=REPO,
                       stdout=open(p.replace('.dat', '.log'), 'w'), stderr=subprocess.STDOUT)
    return name, r.returncode

if __name__ == '__main__':
    cases = [(d, m) for m in MESHES for d in SCHEMES]
    with ThreadPoolExecutor(max_workers=6) as ex:
        for name, rc in ex.map(run, cases): print(name, 'exit', rc, flush=True)
