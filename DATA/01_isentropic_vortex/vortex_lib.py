#!/usr/bin/env python3
# ===========================================================================
#  Isentropic-vortex advection on a WAVY (curvilinear) mesh — energy-preserving
#  (EP) scheme, ep_type = 'JP' (ep_id = 1).
#
#  A 2D isentropic vortex (Kim 2005) is superimposed on a uniform free stream
#  u_inf = sqrt(gamma)*Mach and let translate across a triply periodic box.
#  With ic = 'isentropic_vortex_x' the stream is along x; with '..._y' along y.
#  The vortex is seeded at the PHYSICAL node positions xx(i,j)/yy(i,j) (the IC
#  is curvilinear), so on the wavy mesh the core sits at the true geometric
#  centre and is convected through the curved cells.
#
#  The mesh is the Soldati periodic warp on the in-plane (x,y) = (xi,eta) plane,
#      x = xmin + Lx*xi  + a*sin(2*pi*eta)
#      y = ymin + Ly*eta + b*sin(4*pi*xi)
#  with a = b = AMP chosen so the Jacobian Lx*Ly - 8*pi^2*a*b*cos*cos > 0 every-
#  where (a*b = 0.36 << Lx*Ly/(8 pi^2) = 1.27).  z is the 1D Cartesian extrusion;
#  the vortex is z-invariant so nz is kept small.
#
#  The run stops after exactly one domain crossing (tmax = Lx/u_inf), so the
#  vortex returns to its starting point and make_gif.py produces a seamless loop.
#
#  Helper module, imported by the run_*.py / plot_*.py scripts of this folder.
# ===========================================================================
import os, sys, shutil, subprocess, argparse
import numpy as np

# --------------------------------------------------------------------------- #
#  fixed problem / mesh constants
# --------------------------------------------------------------------------- #
XMIN, XMAX = -5.0, 5.0
YMIN, YMAX = -5.0, 5.0
ZMIN, ZMAX =  0.0, 1.0
LX, LY     = XMAX-XMIN, YMAX-YMIN
NG    = 7                 # ghost layers carried by a curvilinear mesh file (GN+mesh_gh)
AMP   = 0.6               # wavy warp amplitude (a*b = 0.36 < Lx*Ly/(8 pi^2) = 1.27 -> J>0)
N_INPLANE = 100           # in-plane resolution (NxN), the vortex spans ~10 cells
NZ    = 4                 # z extrusion (vortex is z-invariant -> few cells suffice)
MACH  = 0.7
GAMMA = 1.4               # gamma_Tref for the default calorically perfect gas
U_INF = np.sqrt(GAMMA)*MACH
CFL   = 0.5
TMAX  = LX/U_INF          # exactly one domain crossing -> seamless periodic loop
ITOUT = 10                # field dump every ITOUT iterations
ITMAX = 200000            # safety ceiling; the run actually stops on tmax

HERE   = os.path.dirname(os.path.abspath(__file__))
REPO   = os.path.abspath(os.path.join(HERE, '..', '..'))     # DATA/<case> -> DATA -> repo
EXE    = os.path.join(REPO, 'Uranos.exe')
GRIDIR = HERE
DATDIR = HERE


# --------------------------------------------------------------------------- #
#  wavy mesh (ghost-padded, exactly as the solver expects)
# --------------------------------------------------------------------------- #
def build_wavy(nx, ny):
    """ghost-padded XX,YY of shape (nx+2NG, ny+2NG), cell-centred node convention
       s(i)=smin+L*(i-0.5)/n (matching uniformGrid_f0) plus the Soldati warp."""
    ii  = np.arange(1-NG, nx+NG+1)
    jj  = np.arange(1-NG, ny+NG+1)
    csi = (ii - 0.5)/nx
    eta = (jj - 0.5)/ny
    CSI, ETA = np.meshgrid(csi, eta, indexing='ij')
    XX = XMIN + LX*CSI + AMP*np.sin(2.0*np.pi*ETA)
    YY = YMIN + LY*ETA + AMP*np.sin(4.0*np.pi*CSI)
    return XX, YY

def write_mesh_file(path, XX, YY, nx, ny, nz):
    """URANOS curvilinear format: header 'nx ny nz', x-block then y-block,
       i fastest (Fortran order), one value per line."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'{nx} {ny} {nz}\n')
        np.savetxt(f, XX.flatten(order='F')[:, None], fmt='%25.16e')
        np.savetxt(f, YY.flatten(order='F')[:, None], fmt='%25.16e')


# --------------------------------------------------------------------------- #
#  input-file (namelist) template
# --------------------------------------------------------------------------- #
DAT_TEMPLATE = """&input_list
	xmin = {xmin}
	xmax = {xmax}
	ymin = {ymin}
	ymax = {ymax}
	zmin = {zmin}
	zmax = {zmax}
	gridpoint_x = 'uniform'
	gridpoint_y = 'uniform'
	gridpoint_z = 'uniform'
	mesh_file   = '{mesh_file}'
	dims = 3
	tmax = {tmax}
	itmax = {itmax}
	nx = {nx}
	ny = {ny}
	nz = {nz}
	cart_dims(1) = 1
	cart_dims(2) = 1
	cart_dims(3) = 1
	CFL  = {cfl}
	Dt   = 1.0E-02
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
	Reynolds = 10E+10
	Mach = {mach}
	Prandtl = 0.71
	logical_CFL = .true.
	scheme      = 'energy_preserving'
	ep_type     = 'JP'
	sensor      = 'ducros'
	ducrosToll  = 1.0
	fd_order    = 6
	LES         = .false.
	sgs_model   = 'WALE'
	viscous = .false.
	n_step = 3
	itout = {itout}
	StOut = 100000
	StFlg = .false.
	ic = 'isentropic_vortex_{dir}'
	output_file_name = '{ofn}'
	data_dir = '{data_dir}'
&end
"""


# --------------------------------------------------------------------------- #
def run(direction, np_ranks=1):
    if direction not in ('x', 'y'):
        sys.exit(f" --dir must be 'x' or 'y' (got {direction!r})")
    if not os.path.isfile(EXE):
        sys.exit(f' {EXE} not found — build first (make comp=gnu mode=cpu)')

    tag       = f'vortex_wavy_{direction}'
    mesh_file = os.path.join(GRIDIR, f'{tag}.txt')
    XX, YY    = build_wavy(N_INPLANE, N_INPLANE)
    write_mesh_file(mesh_file, XX, YY, N_INPLANE, N_INPLANE, NZ)

    os.makedirs(DATDIR, exist_ok=True)
    data_dir = tag
    outroot  = os.path.join(REPO, 'DATA', data_dir)
    if os.path.isdir(outroot):
        shutil.rmtree(outroot)
    os.makedirs(os.path.join(outroot, 'BINARY'), exist_ok=True)

    dat = DAT_TEMPLATE.format(
        xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX, zmin=ZMIN, zmax=ZMAX,
        mesh_file=mesh_file, tmax=f'{TMAX:.6f}', itmax=ITMAX,
        nx=N_INPLANE, ny=N_INPLANE, nz=NZ, cfl=CFL, mach=MACH, itout=ITOUT,
        dir=direction, ofn=tag, data_dir=data_dir)
    datpath = os.path.join(DATDIR, f'{tag}.dat')
    with open(datpath, 'w') as f:
        f.write(dat)

    print(f' [{tag}] wavy mesh {N_INPLANE}x{N_INPLANE}x{NZ}, AMP={AMP}, '
          f'Mach={MACH}, u_inf={U_INF:.4f}, tmax={TMAX:.3f} (one crossing)')
    print(f' [{tag}] launching solver (np={np_ranks}) ...')
    log = subprocess.run(['mpirun', '-np', str(np_ranks), EXE, datpath],
                         cwd=REPO, capture_output=True, text=True)
    if log.returncode != 0:
        print(log.stdout[-1500:]); print(log.stderr[-1500:])
        sys.exit(f' [{tag}] solver failed (rc={log.returncode})')

    import glob
    ndump = len(glob.glob(os.path.join(outroot, 'BINARY', '*.bin')))
    print(f' [{tag}] done — {ndump} BINARY dumps in DATA/{data_dir}/BINARY')
    return outroot, mesh_file


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, choices=['x', 'y'],
                    help='free-stream direction of the vortex advection')
    ap.add_argument('--np', type=int, default=1, help='MPI ranks')
    args = ap.parse_args()
    run(args.dir, args.np)


if __name__ == '__main__':
    main()
