#!/usr/bin/env python3
# ===========================================================================
#  Multi-GPU scaling of the URANOS-3.0 paper (Sec. 4.2, Figs. 26-27, Table 9):
#  triperiodic flow at rest ('ambient' IC, unit density, pressure and
#  temperature), (2 pi)^3 box with 1024 x 1024 x 512 nodes, extended along z
#  for the weak scaling, no output inside the time loop.
#
#    strong (S): 1024 x 1024 x 512 nodes on 1 ... 256 nodes (4 GPUs each)
#    weak   (W): 1024 x 1024 x (512 n) on n nodes
#    configurations C1 ... C6 (CONFIGS below), 100 timed iterations; the paper
#    uses the mean of three runs after a warm-up of 10 iterations
#
#  Writes decks/scal_<S|W>_<C>_n<nodes>.dat (and the 10-iteration warm-up
#  decks *_w.dat) and grids/wavy_1024.txt in this folder.
#  usage:  python3 DATA/08_multi_gpu_scaling/make_case.py
# ===========================================================================
import os, subprocess, sys
TWO_PI = 6.283185307179586
MESH = "DATA/08_multi_gpu_scaling/grids/wavy_1024.txt"

CONFIGS = {   # name: (scheme, ep_type, weno_order, sensor, viscous, tpg, mesh)
    'C1': ('energy_preserving', 'JP',  5, 'ducros', 'staggered',    0, ''),
    'C2': ('hybrid_wenoEP',     'JP',  5, 'weno',   'staggered',    0, ''),
    'C3': ('hybrid_wenoEP',     'JP',  5, 'weno',   'staggered',    1, ''),
    'C4': ('hybrid_wenoEP',     'JP',  5, 'weno',   'staggered',    0, MESH),
    'C5': ('hybrid_tenoEP',     'AEC', 7, 'weno',   'conservative', 1, MESH),
    'C6': ('energy_preserving', 'JP',  5, 'ducros', 'staggered',    0, MESH),
}

def deck(cfg, nz, niter, data_dir, nx=1024, ny=1024):
    scheme, ep, worder, sensor, visc, tpg, mesh = CONFIGS[cfg]
    lz = TWO_PI * nz / 512
    return f"""&input_list
    ! scaling campaign {cfg}: triperiodic flow at rest, {nx}x{ny}x{nz}
    xmin = 0.
    xmax = {TWO_PI}
    ymin = 0.
    ymax = {TWO_PI}
    zmin = 0.
    zmax = {lz}
    gridpoint_x = 'uniform'
    gridpoint_y = 'uniform'
    gridpoint_z = 'uniform'
    mesh_file   = '{mesh}'
    dims = 3
    nx = {nx}
    ny = {ny}
    nz = {nz}
    cart_dims(1) = 0
    cart_dims(2) = 0
    cart_dims(3) = 0
    mpi_opt_level = 2
    cuda_aware    = .true.
    tmax  = 1.0E+12
    itmax = {niter}
    logical_CFL = .true.
    CFL = 0.5
    Dt  = 1.E-03
    n_step = 3
    bc(1) = 'periodic'
    bc(2) = 'periodic'
    bc(3) = 'periodic'
    bc(4) = 'periodic'
    bc(5) = 'periodic'
    bc(6) = 'periodic'
    Trat = 1.0
    inflow_profile = 'none'
    smooth_inflow  = .false.
    turb_inflow    = .false.
    Reynolds = 1000.0
    Mach     = 0.1
    Prandtl  = 0.72
    viscous  = .true.
    thermally_perfect_gas = {tpg}
    scheme           = '{scheme}'
    ep_type          = '{ep}'
    weno_order       = {worder}
    fd_order         = 6
    sensor           = '{sensor}'
    ducrosToll       = 0.05
    diffusion_scheme = '{visc}'
    enable_limiter   = .false.
    LES       = .false.
    sgs_model = 'WALE'
    ic  = 'ambient'
    printres    = .false.
    restartStat = .false.
    itout = 1000000
    StFlg       = .false.
    rsOut = 1000000
    data_dir         = '08_multi_gpu_scaling/runs/{data_dir}'
    output_file_name = 'scal'
&end
"""


HERE  = os.path.dirname(os.path.abspath(__file__))
NODES = [1, 2, 4, 8, 16, 32, 64, 128, 256]

if __name__ == '__main__':
    mesh = os.path.join(HERE, 'grids', 'wavy_1024.txt')
    if not os.path.isfile(mesh):
        os.makedirs(os.path.dirname(mesh), exist_ok=True)
        subprocess.run([sys.executable, os.path.join(HERE, 'gen_curved_mesh.py'), '1024', '1024', '512',
                        '0', str(TWO_PI), '0', str(TWO_PI), '0.2', mesh], check=True)
    os.makedirs(os.path.join(HERE, 'decks'), exist_ok=True)
    n_decks = 0
    for c in CONFIGS:
        for n in NODES:
            for s in 'SW':
                nz = 512 if s == 'S' else 512 * n
                for suf, it in (('', 100), ('_w', 10)):
                    with open(os.path.join(HERE, 'decks', f'scal_{s}_{c}_n{n}{suf}.dat'), 'w') as f:
                        f.write(deck(c, nz, it, f'scal_{s}_{c}_n{n}'))
                    n_decks += 1
    print(n_decks, 'decks in', os.path.join(HERE, 'decks'))
