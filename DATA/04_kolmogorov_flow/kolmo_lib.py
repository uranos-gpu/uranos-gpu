#!/usr/bin/env python3
# ===========================================================================
#  Error-scaling harness for the VISCOUS (Laplacian-form) operator on
#  curvilinear grids — steady Kolmogorov flow.
#
#  Body force (rhs_module::kolmogorov_flow_forcing), curvilinear-aware:
#       f_x(y) = 0.1 u_inf sin(8 pi yy / Ly) ,   yy = PHYSICAL y
#  Pure-viscous problem (scheme='none' zeroes the convective RHS,
#  diffusion_scheme='laplacian' selects the ported chain-rule operator).
#  The steady state is the EXACT solution of the (solenoidal, unidirectional)
#  balance  mu_inf*VIS*d2u/dy2 + f_x = 0 :
#       u*(yy) = A sin(8 pi yy / Ly) ,   A = 0.1 Re /(VIS k^2) = 1/(16 pi^2)
#  with k = 8 pi / Ly = 4 pi (Ly=2), Re=10 (=> 0.1*Re=1), VIS->1 at T=Tref.
#  Because u depends only on the PHYSICAL y, du/dx_phys = 0 (metric identity),
#  v=w=0 are preserved, so u* is exact on ANY mesh (cartesian / fake / wavy);
#  the steady solver error is purely the viscous-operator truncation error.
#
#  Norms (paper De Vanna et al. CPC 287 (2023) eq.20, volume-weighted):
#       Lp = ( sum_i |u_i-u*_i|^p dV_i / sum_i dV_i )^(1/p) ,  Linf = max|.|
#  dV = detJ2D/(nx*ny) * dz  (physical cell volume); on a uniform mesh -> (1/N).
#
#  fd_order = {2,4,6} sets the ORDER of the viscous scheme: storage_module
#  fills central_1 (1st deriv) and central_2 (2nd deriv) with the 2/4/6-order
#  stencils, so fd_order is the formal order of the whole Laplacian operator.
#
#  Three mesh families, same as the advection (periodic_euler) study:
#    * real_cart : separable 1D-tensor Cartesian grid     (mesh_file = '')
#    * fake_cart : 2D curvilinear mesh file, Cartesian geometry (uniform)
#                  -> exercises the curvilinear metric path; must match real_cart
#    * wavy      : genuinely curved periodic map (Soldati), u* remapped on yy
#
#  Steady state is reached deterministically by stopping on ITMAX ~ C*N^2
#  (dt is viscous-limited, dt ~ dy^2, so steps-to-steady ~ N^2); the final
#  dump (itout = itmax) is the steady field.
#
#  usage (from repo root):
#    python3 .../run_scaling.py --smoke
#    python3 .../run_scaling.py --jobs 8 --orders 2,4,6 --meshes real_cart,fake_cart,wavy
# ===========================================================================
import os, sys, glob, json, shutil, subprocess, argparse, math
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

# --------------------------------------------------------------------------- #
#  fixed problem / discretisation constants
# --------------------------------------------------------------------------- #
XMIN, XMAX = -1.0, 1.0
YMIN, YMAX = -1.0, 1.0
ZMIN, ZMAX = -1.0, 1.0
LX, LY, LZ = XMAX-XMIN, YMAX-YMIN, ZMAX-ZMIN
NG    = 7                       # ghost layers in a curvilinear mesh file (GN+mesh_gh)
RE    = 10.0
MACH  = 0.05
FORCING = 1.0e-3               # Kolmogorov body-force amplitude. LINEAR REGIME: the
                               # compressible variable-mu(T) coupling generates a steady
                               # 3rd harmonic + amplitude correction ~ forcing^2, so a
                               # small forcing keeps u* a single sine and that physical
                               # correction sits far below the truncation error.
                               # (production default is 0.1; verified analytically that
                               #  the order is forcing-independent.)
KWAVE = 8.0*np.pi/LY           # forcing wavenumber = 4 pi
A_AN  = FORCING*RE/(KWAVE*KWAVE)   # = forcing/(16 pi^2 / 10)  (VIS=1)
CFL   = 0.5
BASE  = 8                       # transverse points on a cartesian run
WAVY_AX = 0.10                  # Soldati warp amplitudes (J>0 guaranteed)
WAVY_AY = 0.10
STEADY_C_CART = 3.0             # itmax = C * N^2  (steady ~2.2 N^2, margin)
STEADY_C_WAVY = 4.0             # wavy cells are compressed -> smaller dt

C1 = np.array([-1.0/60, 3.0/20, -3.0/4, 0.0, 3.0/4, -3.0/20, 1.0/60])  # order-6 d/dxi

HERE   = os.path.dirname(os.path.abspath(__file__))
REPO   = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE    = os.path.join(REPO, 'Uranos.exe')
GRIDIR = os.path.join(HERE, 'grids_scaling')
DATDIR = os.path.join(HERE, 'decks_scaling')
RESDIR = os.path.join(HERE, 'results')


# --------------------------------------------------------------------------- #
#  mesh construction (ghost-padded, cell-centred — matches uniformGrid_f0)
# --------------------------------------------------------------------------- #
def grid_indices(n):
    return np.arange(1-NG, n+NG+1)

def build_mesh(meshtype, nx, ny):
    ii, jj = grid_indices(nx), grid_indices(ny)
    csi = (ii - 0.5)/nx
    eta = (jj - 0.5)/ny
    CSI, ETA = np.meshgrid(csi, eta, indexing='ij')
    if meshtype == 'wavy':
        XX = XMIN + LX*CSI + WAVY_AX*np.sin(2.0*np.pi*ETA)
        YY = YMIN + LY*ETA + WAVY_AY*np.sin(4.0*np.pi*CSI)
    else:                                    # real_cart / fake_cart geometry
        XX = XMIN + LX*CSI
        YY = YMIN + LY*ETA
    return XX, YY

def write_mesh_file(path, XX, YY, nx, ny, nz):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'{nx} {ny} {nz}\n')
        np.savetxt(f, XX.flatten(order='F')[:, None], fmt='%25.16e')
        np.savetxt(f, YY.flatten(order='F')[:, None], fmt='%25.16e')


# --------------------------------------------------------------------------- #
#  metrics / volume element (same order-6 central FD as the solver)
# --------------------------------------------------------------------------- #
def d_idx(a, axis):
    out = np.zeros_like(a)
    for l in range(-3, 4):
        out += C1[l+3]*np.roll(a, -l, axis=axis)
    return out

def detJ2D_interior(XX, YY, nx, ny):
    x_xi, x_eta = d_idx(XX, 0), d_idx(XX, 1)
    y_xi, y_eta = d_idx(YY, 0), d_idx(YY, 1)
    sl = (slice(NG, NG+nx), slice(NG, NG+ny))
    return x_xi[sl]*y_eta[sl] - x_eta[sl]*y_xi[sl]      # (nx,ny), per-index Jacobian


# --------------------------------------------------------------------------- #
#  URANOS BINARY reader  (phi(1)=rho, phi(2)=rho*u, ... -> u)
# --------------------------------------------------------------------------- #
def read_binary_u(fn):
    with open(fn, 'rb') as f:
        nx, ny, nz, it = np.fromfile(f, np.int32,  4)
        time, dt       = np.fromfile(f, np.float64, 2)
        n   = int(nx)*int(ny)*int(nz)
        flat = np.fromfile(f, np.float64, 5*n)
    rho  = flat[0:n   ].reshape((nx, ny, nz), order='F')
    rhou = flat[n:2*n ].reshape((nx, ny, nz), order='F')
    return float(time), int(nx), int(ny), int(nz), rhou/rho

def final_dump(bindir):
    files = glob.glob(os.path.join(bindir, '*.bin'))
    if not files:
        raise RuntimeError(f'no .bin in {bindir}')
    best, bestit = None, -1
    for fn in files:
        with open(fn, 'rb') as f:
            it = int(np.fromfile(f, np.int32, 4)[3])
        if it > bestit:
            best, bestit = fn, it
    return best


# --------------------------------------------------------------------------- #
#  input-file template (pure-viscous Kolmogorov, Laplacian diffusion)
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
	tmax = 1.0E+09
	itmax = {itmax}
	nx = {nx}
	ny = {ny}
	nz = {nz}
	cart_dims(1) = 1
	cart_dims(2) = 1
	cart_dims(3) = 1
	CFL  = {cfl}
	Dt   = 1.0E-03
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
	Reynolds = {re}
	Mach = {mach}
	Prandtl = 0.71
	logical_CFL = .true.
	scheme      = 'none'
	diffusion_scheme = '{diffusion}'
	kolmogorov_forcing = {forcing}
	sensor      = 'density'
	fd_order    = {fd_order}
	LES         = .false.
	sgs_model   = 'WALE'
	viscous = .true.
	n_step = 3
	itout = {itmax}
	StOut = 100000000
	StFlg = .false.
	ic = 'KolmogorovFlow'
	output_file_name = '{ofn}'
	data_dir = '{data_dir}'
&end
"""


# --------------------------------------------------------------------------- #
def case_dims(meshtype, N):
    """cartesian/fake: the steady solution is 1D in y -> refine ny only.
       wavy: yy(x,y) varies in both indices -> refine in-plane isotropically."""
    if meshtype == 'wavy':
        return N, N, BASE
    return BASE, N, BASE          # real_cart / fake_cart : refine ny

def itmax_for(meshtype, N):
    C = STEADY_C_WAVY if meshtype == 'wavy' else STEADY_C_CART
    return max(3000, int(math.ceil(C*N*N)))


def run_one(meshtype, N, fd_order=6, verbose=False, diffusion='laplacian'):
    tag = f'{diffusion[:4]}_{meshtype}_O{fd_order}_N{N}'
    nx, ny, nz = case_dims(meshtype, N)
    itmax = itmax_for(meshtype, N)

    XX, YY = build_mesh(meshtype, nx, ny)
    mesh_file = '' if meshtype == 'real_cart' else os.path.join(GRIDIR, f'{tag}.txt')
    if mesh_file:
        write_mesh_file(mesh_file, XX, YY, nx, ny, nz)

    os.makedirs(DATDIR, exist_ok=True)
    data_dir = f'visc_scaling_{tag}'
    outroot  = os.path.join(REPO, 'DATA', data_dir)
    shutil.rmtree(outroot, ignore_errors=True)
    os.makedirs(os.path.join(outroot, 'BINARY'), exist_ok=True)

    dat = DAT_TEMPLATE.format(
        xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX, zmin=ZMIN, zmax=ZMAX,
        mesh_file=mesh_file, itmax=itmax, nx=nx, ny=ny, nz=nz, cfl=CFL,
        re=RE, mach=MACH, forcing=f'{FORCING:.6E}', fd_order=fd_order, ofn=tag, data_dir=data_dir,
        diffusion=diffusion)
    datpath = os.path.join(DATDIR, f'{tag}.dat')
    with open(datpath, 'w') as f:
        f.write(dat)

    log = subprocess.run(['mpirun', '-np', '1', EXE, datpath],
                         cwd=REPO, capture_output=True, text=True)
    if log.returncode != 0:
        return {'tag': tag, 'error': f'solver rc={log.returncode}',
                'stderr': log.stderr[-800:], 'stdout': log.stdout[-600:]}

    fn = final_dump(os.path.join(outroot, 'BINARY'))
    t, bnx, bny, bnz, u = read_binary_u(fn)

    # exact steady solution u*(yy) = A sin(k yy) at the physical nodes
    yy = YY[NG:NG+nx, NG:NG+ny][:, :, None]
    u_exact = A_AN*np.sin(KWAVE*yy) + np.zeros_like(u)

    # physical cell-volume weight
    detJ = detJ2D_interior(XX, YY, nx, ny)
    dV   = (detJ[:, :, None]/(nx*ny))*(LZ/nz) + np.zeros_like(u)

    err  = np.abs(u - u_exact)
    Vtot = np.sum(dV)
    L1   = np.sum(err*dV)/Vtot
    L2   = np.sqrt(np.sum(err*err*dV)/Vtot)
    Linf = float(np.max(err))
    # amplitude diagnostic (best-fit to sin(k yy))
    s = np.sin(KWAVE*yy) + np.zeros_like(u)
    A_fit = float(np.sum(u*s*dV)/np.sum(s*s*dV))
    h    = LY/N

    if verbose:
        print(f'  {tag:22s} it_steady={itmax:6d} t={t:8.2f}  '
              f'L1={L1:.3e} L2={L2:.3e} Linf={Linf:.3e}  A/A_an={A_fit/A_AN:.5f}')
    return {'tag': tag, 'meshtype': meshtype, 'fd_order': fd_order, 'N': N,
            'nx': nx, 'ny': ny, 'nz': nz, 't': t, 'h': h,
            'L1': L1, 'L2': L2, 'Linf': Linf, 'A_ratio': A_fit/A_AN}


def fit_order(hs, errs, floor=1e-12):
    hs, errs = np.asarray(hs, float), np.asarray(errs, float)
    m = errs > floor
    if m.sum() < 2:
        m = np.ones_like(errs, bool)
    return np.polyfit(np.log(hs[m]), np.log(errs[m]), 1)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 4)//2))
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--meshes', default='real_cart,fake_cart,wavy')
    ap.add_argument('--orders', default='2,4,6')
    ap.add_argument('--diffusion', default='laplacian', choices=['laplacian', 'staggered'])
    args = ap.parse_args()

    if not os.path.isfile(EXE):
        sys.exit(f' {EXE} not found — build first (make comp=gnu mode=cpu)')
    os.makedirs(RESDIR, exist_ok=True)

    CART_LADDER = [24, 32, 48, 64, 96, 128]
    WAVY_LADDER = [24, 32, 48, 64, 96]

    if args.smoke:
        for N in (32, 64):
            print(json.dumps(run_one('real_cart', N, 6, verbose=True, diffusion=args.diffusion)))
        return

    meshes = args.meshes.split(',')
    orders = [int(o) for o in args.orders.split(',')]
    jobs = []
    for fd in orders:
        for meshtype in meshes:
            ladder = WAVY_LADDER if meshtype == 'wavy' else CART_LADDER
            for N in ladder:
                jobs.append((meshtype, N, fd))

    print(f' launching {len(jobs)} cases on {args.jobs} workers ...')
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_one, m, N, fd, True, args.diffusion): (m, N, fd) for m, N, fd in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            if 'error' in r:
                print(f'  !! {r["tag"]}: {r["error"]}\n{r.get("stderr","")}')
            results.append(r)

    results = [r for r in results if 'error' not in r]
    sfx = '' if args.diffusion == 'laplacian' else f'_{args.diffusion}'
    out = os.path.join(RESDIR, f'scaling_results{sfx}.json')
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\n wrote {out}  ({len(results)} successful cases)')

    print('\n ===== observed order of accuracy (slope of log-error vs log-h) =====')
    print(f'  {"fd":>2} {"mesh":>10} | {"order L1":>9} {"order L2":>9} {"order Linf":>10}')
    for fd in orders:
        for meshtype in meshes:
            sub = sorted([r for r in results if r['fd_order'] == fd and r['meshtype'] == meshtype],
                         key=lambda r: r['h'])
            if len(sub) < 2:
                continue
            hs = [r['h'] for r in sub]
            o1 = fit_order(hs, [r['L1']   for r in sub])
            o2 = fit_order(hs, [r['L2']   for r in sub])
            oi = fit_order(hs, [r['Linf'] for r in sub])
            print(f'  {fd:>2} {meshtype:>10} | {o1:9.2f} {o2:9.2f} {oi:10.2f}')


if __name__ == '__main__':
    main()
