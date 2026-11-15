#!/usr/bin/env python3
# ===========================================================================
#  Taylor-Green vortex campaign of the URANOS-3.0 paper (Sec. 3.2.2, Fig. 10,
#  Table 5): every convective scheme x viscous scheme x resolution x grid.
#
#    convective + gas : JP, KEEP, AEC with the calorically (CP) and thermally
#                       perfect (TP) gas, PEP with TP only            -> 7
#    viscous          : laplacian, staggered, conservative           -> 3
#    resolution       : 64^3, 128^3, 256^3, 512^3                     -> 4
#    grid             : Cartesian (analytic), wavy (mesh file)        -> 2
#                                                              total  168
#
#  Canonical compressible TGV (van Rees 2011): Re_flow = 1600, Ma = 0.1,
#  Pr = 0.71, triply periodic (2 pi)^3, up to t* = 20.  The IC uses p_z = 100,
#  so U0 = 10 u_inf and Reynolds(namelist) = 160; the viscosity is constant
#  (vis_flag = 2).  Wavy grid: gen_curved_mesh.py, amp = 0.2 (separable warp).
#
#  Writes grids/wavy_tgv_<N>.txt, decks/<tag>.dat and campaign.json in this
#  folder (already shipped: rerun only to regenerate them).
#  usage:  python3 DATA/05_taylor_green_vortex/make_case.py
# ===========================================================================
import os, json, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
NAME = os.path.basename(HERE)            # 05_taylor_green_vortex
CAMP = HERE
GRID = os.path.join(HERE, 'grids')

TWO_PI = 6.283185307179586
U0 = 1.1832159566199232          # Mach*sqrt(gamma*p_z), t* = t*U0
TSTAR_END = 20.0
AMP = 0.2

CONV = [('JP', 0), ('KEEP', 0), ('AEC', 0),
        ('JP', 1), ('KEEP', 1), ('AEC', 1), ('PEP', 1)]
VISC = {'laplacian': 'lap', 'staggered': 'stag', 'conservative': 'cons'}
RES = [64, 128, 256, 512]
GRIDS = ['cart', 'wavy']

ITOUT = {64: 200000, 128: 200000, 256: 20000, 512: 10000}
STOUT = {64: 20, 128: 20, 256: 50, 512: 50}

DAT = """&input_list
    xmin = 0.
    xmax = {L}
    ymin = 0.
    ymax = {L}
    zmin = 0.
    zmax = {L}
    gridpoint_x = 'uniform'
    gridpoint_y = 'uniform'
    gridpoint_z = 'uniform'
    mesh_file = '{mesh}'
    dims = 3
    nx = {n}
    ny = {n}
    nz = {n}
    cart_dims(1) = 0
    cart_dims(2) = 0
    cart_dims(3) = 0
    tmax  = {tmax:.6f}
    itmax = 100000000
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
    Reynolds = 160          ! Re_flow = 10 x Reynolds = 1600
    Mach     = 0.1
    Prandtl  = 0.71
    vis_flag = 2            ! constant mu (mandatory for the TGV IC)
    thermally_perfect_gas = {tp}
    scheme           = 'energy_preserving'
    ep_type          = '{ep}'
    diffusion_scheme = '{diff}'
    sensor           = 'density_jump'
    fd_order         = 6
    weno_order       = 5
    LES              = .false.
    sgs_model        = 'WALE'
    viscous = .true.
    itout = {itout}
    StOut = {stout}
    StFlg = .true.
    rsOut = 200
    ic               = 'TGV'
    output_file_name = 'TGV'
    data_dir         = '{data_dir}'
    mpi_opt_level = 2
    cuda_aware    = .true.
&end
"""

def tag_of(grid, ep, tp, diff, n):
    return f"{grid}_{'TP' if tp else 'CP'}_{ep}_{VISC[diff]}_N{n:04d}"


def main():
    for d in ('decks',):
        os.makedirs(os.path.join(CAMP, d), exist_ok=True)
    os.makedirs(GRID, exist_ok=True)
    gen = os.path.join(HERE, 'gen_curved_mesh.py')
    meshes = {}
    for n in RES:
        m = os.path.join(GRID, f'wavy_tgv_{n}.txt')
        if not os.path.isfile(m):
            subprocess.run([sys.executable, gen, str(n), str(n), str(n), '0', str(TWO_PI),
                            '0', str(TWO_PI), str(AMP), m], check=True, cwd=REPO)
        meshes[n] = os.path.relpath(m, REPO)

    cases = []
    for grid in GRIDS:
        for ep, tp in CONV:
            for n in RES:
                for diff in VISC:
                    tag = tag_of(grid, ep, tp, diff, n)
                    with open(os.path.join(CAMP, 'decks', f'{tag}.dat'), 'w') as f:
                        f.write(DAT.format(L=TWO_PI, n=n, tmax=TSTAR_END / U0,
                                           mesh=meshes[n] if grid == 'wavy' else '',
                                           tp=tp, ep=ep, diff=diff, itout=ITOUT[n],
                                           stout=STOUT[n], data_dir=f'{NAME}/out/{tag}'))
                    cases.append(dict(tag=tag, grid=grid, gas='TP' if tp else 'CP', conv=ep,
                                      visc=diff, N=n))
    with open(os.path.join(CAMP, 'campaign.json'), 'w') as f:
        json.dump(dict(U0=U0, tstar_end=TSTAR_END, amp=AMP, cases=cases), f, indent=1)
    print(f'{len(cases)} decks in {CAMP}/decks')


if __name__ == '__main__':
    main()
