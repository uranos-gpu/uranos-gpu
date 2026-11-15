#!/usr/bin/env python3
# Sod shock tube runs for the URANOS-3.0 paper figure (WENO/TENO comparison).
# Everything (mesh, decks, logs, outputs) lives in DATA/03_sod_shock_tube, so the
# whole folder can be deleted in one go.  Run from anywhere:
#   python3 DATA/03_sod_shock_tube/run_sod.py
import os, re, subprocess, shutil
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE  = os.path.join(REPO, 'Uranos.exe')
BASE = os.path.join(HERE, 'shock_tube_base.dat')
REL  = os.path.relpath(HERE, os.path.join(REPO, 'DATA'))      # 03_sod_shock_tube
NX, NY, NZ, TEND, DT = 200, 5, 5, 1.4, 1.0e-3
MESH = os.path.join(HERE, 'mesh_fake_cartesian.txt')

CASES = [('weno3', 'hybrid_wenoEP', 3), ('weno5', 'hybrid_wenoEP', 5),
         ('weno7', 'hybrid_wenoEP', 7), ('teno5', 'hybrid_tenoEP', 5),
         ('teno7', 'hybrid_tenoEP', 7), ('tenoa5', 'hybrid_tenoaEP', 5)]

def make_mesh():
    subprocess.run(['python3', os.path.join(HERE, 'gen_fake_cartesian_mesh.py'),
                    str(NX), str(NY), str(NZ), '-5', '5', '-1', '1', MESH],
                   check=True, stdout=subprocess.DEVNULL)

def run(case, gas=0, suffix='', tref=298.15):
    name, scheme, order = case
    name = name + suffix
    txt = open(BASE).read()
    nit = int(round(TEND / DT))
    ov = (f"\tmesh_file = '{os.path.relpath(MESH, REPO)}'\n\tdata_dir = \"{REL}/{name}\"\n"
          f"\tnx={NX}\n\tny={NY}\n\tnz={NZ}\n\tthermally_perfect_gas={gas}\n\tTref={tref}\n"
          f"\tscheme='{scheme}'\n\tfd_order=6\n\tweno_order={order}\n\tsensor='density_jump'\n"
          f"\tlogical_CFL=.false.\n\tDt={DT:.1E}\n\titmax={nit}\n\titout={nit}\n"
          f"\tcart_dims(1)=1\n\tcart_dims(2)=1\n\tcart_dims(3)=1\n")
    txt = re.sub(r'&end', ov + '&end', txt, count=1)
    dat = os.path.join(HERE, f'{name}.dat'); open(dat, 'w').write(txt)
    shutil.rmtree(os.path.join(HERE, name), ignore_errors=True)
    r = subprocess.run(['mpirun', '--oversubscribe', '-np', '1', EXE, dat], cwd=REPO,
                       stdout=open(dat.replace('.dat', '.log'), 'w'), stderr=subprocess.STDOUT)
    return name, r.returncode

if __name__ == '__main__':
    make_mesh()
    import sys
    gas = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    tref = float(sys.argv[2]) if len(sys.argv) > 2 else 298.15
    suffix = sys.argv[3] if len(sys.argv) > 3 else ('' if gas == 0 else '_tp')
    with ThreadPoolExecutor(max_workers=6) as ex:
        for name, rc in ex.map(lambda c: run(c, gas, suffix, tref), CASES):
            print(name, 'exit', rc, flush=True)
