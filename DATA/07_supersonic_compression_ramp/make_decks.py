#!/usr/bin/env python3
# ===========================================================================
#  Decks of the 24-deg compression ramp at Ma = 2.91 (URANOS-3.0 paper,
#  Sec. 3.5, Table 7, DNS rows), on the body-fitted DNS mesh
#  mesh/dns/MeshXY_ramp.txt (2560 x 384 x 192, -29 <= x/delta0 <= 12,
#  wall-normal extent 6 delta0, span 2.2 delta0).
#
#    KEEP-CP  hybrid WENO-Z 5, Laplacian viscous scheme
#    AEC-TP   hybrid WENO-Z 7, conservative viscous scheme
#    PEP-TP   hybrid WENO-Z 7, conservative viscous scheme
#
#  Common setup: Ducros sensor (threshold 0.2), sixth-order central fluxes,
#  RK3 at CFL = 0.5, digital-filter turbulent inflow, NSCBC outflow,
#  supersonic (extrapolation) top boundary, isothermal wall Tw/Tr = 1.14
#  (Tw/T_inf = 2.87), T_inf = 107.1 K, spanwise periodicity.
#
#  Each scheme has two decks, chained by a restart:
#    <case>_init.dat   transient, statistics off
#    <case>_stats.dat  restarts from the last field of the transient and
#                      accumulates the statistics (StOut)
#
#  usage:  python3 DATA/07_supersonic_compression_ramp/make_decks.py
# ===========================================================================
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.basename(HERE)
MESH = f'DATA/{CASE}/mesh/dns/MeshXY_ramp.txt'

GRID = dict(nx=2560, ny=384, nz=192, xmin=-29.0, xmax=12.0, ymin=0.0, ymax=6.0,
            zmin=-1.1, zmax=1.1, cart=(8, 1, 4))      # 32 GPUs; cart_dims(2) must be 1
INFLOW_RETAU = 380.0     # digital-filter inflow Re_tau (deck 'Reynolds' for ic='turbulent_BL')
ITERS = {'init': 400000, 'stats': 2000000}

CASES = {   # name: (ep_type, thermally_perfect_gas, weno_order, diffusion_scheme)
    'dns_keep_cp': ('KEEP', 0, 5, 'laplacian'),
    'dns_aec_tp':  ('AEC',  1, 7, 'conservative'),
    'dns_pep_tp':  ('PEP',  1, 7, 'conservative'),
}

DECK = """&input_list

! Mach 2.91 / 24-deg compression ramp (URANOS-3.0 paper, Sec. 3.5): {name}, {phase}

! --- body-fitted mesh (x-y from mesh_file, z uniform) ---
mesh_file = '{mesh}'
xmin = {xmin}
xmax = {xmax}
ymin = {ymin}
ymax = {ymax}
zmin = {zmin}
zmax = {zmax}
gridpoint_x = 'uniform'
gridpoint_y = 'uniform'
gridpoint_z = 'uniform'
dims = 3
nx = {nx}
ny = {ny}
nz = {nz}
! cart_dims(2) = 1 is required by the turbulent-inflow wall helper
cart_dims(1) = {c1}
cart_dims(2) = {c2}
cart_dims(3) = {c3}

CFL  = 0.5
Dt   = 1.E-04
logical_CFL = .true.
n_step = 3

bc(1) = 'nscbc_inflow_generalized'
bc(2) = 'nscbc_outflow'
bc(3) = 'diabatic_wall'
bc(4) = 'extrapolation'
bc(5) = 'periodic'
bc(6) = 'periodic'

! isothermal wall: Tw = Trat*T_recovery = 2.87 T_inf
Trat = 1.14

! digital-filter turbulent inflow
inflow_profile = 'TurbulentBoundaryLayerInflow'
smooth_inflow  = .false.
turb_inflow    = .true.
turbulent_intensity = 0.03

! for ic = 'turbulent_BL' Reynolds is the INFLOW Re_tau
Reynolds = {retau}
Mach     = 2.91
Prandtl  = 0.72

thermally_perfect_gas = {tp}
Tref     = 107.1
vis_flag = 1

scheme           = 'hybrid_wenoEP'
ep_type          = '{ep}'
weno_order       = {weno}
diffusion_scheme = '{diff}'
sensor           = 'ducros'
ducrosToll       = 0.2
fd_order         = 6

enable_limiter     = .true.
limiter_rho_floor  = 0.001
limiter_e_floor    = 0.0001
limiter_abort_frac = 0.0001

LES       = .false.
sgs_model = 'WALE'
viscous   = .true.

tmax  = 1.E9
itmax = {itmax}
StFlg = {stflg}
restartStat = {stflg}
{itprb}itout = 5000
StOut = 2000
rsOut = 1000

ic = 'turbulent_BL'
output_file_name = 'ramp'
data_dir = '{data_dir}'



&end
"""

if __name__ == '__main__':
    g = GRID
    for name, (ep, tp, weno, diff) in CASES.items():
        for phase in ('init', 'stats'):
            stats = phase == 'stats'
            txt = DECK.format(name=name, phase=phase, mesh=MESH, xmin=g['xmin'], xmax=g['xmax'],
                              ymin=g['ymin'], ymax=g['ymax'], zmin=g['zmin'], zmax=g['zmax'],
                              nx=g['nx'], ny=g['ny'], nz=g['nz'], c1=g['cart'][0], c2=g['cart'][1],
                              c3=g['cart'][2], retau=INFLOW_RETAU, tp=tp, ep=ep, weno=weno, diff=diff,
                              itmax=ITERS[phase], stflg='.true.' if stats else '.false.',
                              itprb='itPrb = 0\n' if stats else '', data_dir=f'{CASE}/{name}',
                              rey='.true.' if stats else '.false.')
            with open(os.path.join(HERE, f'{name}_{phase}.dat'), 'w') as f:
                f.write(txt)
            print('wrote', f'{name}_{phase}.dat')
