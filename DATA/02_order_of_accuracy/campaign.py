#!/usr/bin/env python3
# ===========================================================================
#  URANOS 3.0 -- error-scaling campaign over every PURE scheme.
#
#  Three blocks, 324 curves, 2349 runs.  One curve = one ladder in N = one
#  output directory = one `scaling.txt` ready to plot.
#
#    A_ep       pure energy-preserving split forms.  Periodic density wave
#               (De Vanna et al., CPC 287 (2023) 108717, A.1):
#                   (rho,u,v,w,p)_0 = (1 + 0.2 sin(pi s), 0.05, 0, 0, 1)
#               exact:  rho(s,t) = 1 + 0.2 sin(pi (s - 0.05 t)).
#               Valid for a thermally-perfect gas too: with p and u uniform
#               rho*e is advected passively and is a function of rho alone.
#               ep_type in {JP, KEEP, AEC, PEP} x gas {CP, TP} x fd_order
#               {2,4,6}.  AEC auto-branches on thermally_perfect_gas (AEC2-CP
#               / AEC-TP); PEP is thermally-perfect only, the calorically-
#               perfect path is rejected by input_module -- so CP x PEP is
#               absent by design, not by omission.
#
#    B_weno     pure WENO/TENO.  Same density wave with sensor='weno', which
#               forces every cell into the shock branch, so the EP half never
#               runs and the reconstruction is measured alone.
#               WENO {3,5,7}, TENO {5,7}, TENO-A {5} x gas {CP, TP}.
#
#    C_viscous  pure viscous operators.  Steady Kolmogorov flow with
#               scheme='none':  u*(y) = (forcing*Re/k^2) sin(k y), k = 8 pi/Ly,
#               exact on any mesh because it depends on the physical ordinate
#               only.  Steady state is reached by stopping on itmax ~ C N^2, so
#               the time integrator contributes no error at all.
#               diffusion_scheme in {laplacian, staggered, conservative} x
#               fd_order {2,4,6}.  Sutherland mu(T) throughout (vis_flag=1).
#
#  Every block runs on three mesh families:
#    real_cart : separable tensor grid, no mesh file
#    fake_cart : curvilinear mesh file whose geometry is Cartesian -- exercises
#                the metric path and must reproduce real_cart to round-off
#    wavy      : genuinely curved periodic (Soldati-type) map
#
#  Errors are volume-weighted Lp norms with the physical cell volume as weight,
#  which collapses to the paper's (1/N) sum on a uniform mesh.
#
#  BUILD REQUIREMENT -- OpenACC, never CUDA:
#      make comp=nvhpc mode=gpu gpu_cc=cc80
#  With -DUSE_CUDA the advection backend ignores ep_id and always runs the
#  hard-wired JP / calorically-perfect / Cartesian split, so blocks A and B
#  would silently produce identical curves; and diffusion_scheme='conservative'
#  has no CUDA backend at all.
#
#  Everything lives in DATA/02_order_of_accuracy: decks/, grids/ (shipped),
#  runs/, logs/, results/ (created when the runs are made).  Run from anywhere.
#
#  usage:
#    campaign.py list                       # the 324 curve names
#    campaign.py plan                       # run/curve counts and totals
#    campaign.py prepare                    # write every deck and mesh file (no run)
#    campaign.py run --curve <name>         # one curve (what a SLURM task does)
#    campaign.py run --index <i> --block A  # same, addressed by index
#    campaign.py gate                       # fake_cart == real_cart check
# ===========================================================================
import os, sys, glob, json, math, shutil, argparse, subprocess, datetime
import numpy as np

# --------------------------------------------------------------------------- #
#  paths
# --------------------------------------------------------------------------- #
HERE      = os.path.dirname(os.path.abspath(__file__))
REPO      = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE       = os.path.join(REPO, 'Uranos.exe')
CAMPAIGN  = os.path.basename(HERE)                 # 02_order_of_accuracy (DATA/<CAMPAIGN>)
DATA_ROOT = os.path.join(HERE, 'runs')             # solver output, one folder per run
GRID_ROOT = os.path.join(HERE, 'grids')            # mesh files, one per geometry
INP_ROOT  = os.path.join(HERE, 'decks')            # one deck per run
LOG_ROOT  = os.path.join(HERE, 'logs')
RES_ROOT  = os.path.join(HERE, 'results')

BLOCK_DIR = {'A': 'A_ep', 'B': 'B_weno', 'C': 'C_viscous'}

# --------------------------------------------------------------------------- #
#  geometry (shared by all three blocks)
# --------------------------------------------------------------------------- #
XMIN, XMAX = -1.0, 1.0
YMIN, YMAX = -1.0, 1.0
ZMIN, ZMAX = -1.0, 1.0
LX, LY, LZ = XMAX-XMIN, YMAX-YMIN, ZMAX-ZMIN
NG = 7                     # ghost layers carried by a curvilinear mesh file (GN+mesh_gh)
WAVY_AX = 0.10             # Soldati warp amplitudes; a*b < Lx*Ly/(8 pi^2) keeps J > 0
WAVY_AY = 0.10

# order-6 central first derivative, unit index spacing (the solver's own stencil)
C1 = np.array([-1.0/60, 3.0/20, -3.0/4, 0.0, 3.0/4, -3.0/20, 1.0/60])

# --------------------------------------------------------------------------- #
#  blocks A/B -- advected density wave
# --------------------------------------------------------------------------- #
A_WAVE   = 0.2
U_INF    = 0.05            # = Mach
DT_ADV   = 1.0e-3          # fixed Dt: the RK3 error ~1e-13 stays under the spatial one
ITMAX_ADV = 1000           # t_final = 1.0, stop on itmax (deterministic)
BASE_ADV = 5               # transverse points where the solution is uniform

# --------------------------------------------------------------------------- #
#  block C -- steady Kolmogorov flow
# --------------------------------------------------------------------------- #
RE_VIS   = 10.0
MACH_VIS = 0.05
FORCING  = 1.0e-3          # linear regime: the mu(T) coupling makes a 3rd harmonic
                           # ~forcing^2, which must stay below the truncation error
KWAVE    = 8.0*np.pi/LY
A_AN     = FORCING*RE_VIS/(KWAVE*KWAVE)
CFL_VIS  = 0.5
BASE_VIS = 8
STEADY_C_CART = 3.0        # itmax = C N^2 (steady at ~2.2 N^2, margin on top)
STEADY_C_WAVY = 4.0        # wavy cells are compressed -> smaller dt -> more steps

# --------------------------------------------------------------------------- #
#  ladders
# --------------------------------------------------------------------------- #
LADDER = {
    'A': {'cart': [16, 24, 32, 48, 64, 96, 128], 'wavy': [16, 24, 32, 48, 64, 96]},
    'B': {'cart': [8, 12, 16, 24, 32, 48, 64, 96, 128], 'wavy': [8, 12, 16, 24, 32, 48, 64, 96]},
    'C': {'cart': [24, 32, 48, 64, 96, 128], 'wavy': [24, 32, 48, 64, 96]},
}

MESHES = ['real_cart', 'fake_cart', 'wavy']
DIRS   = ['x', 'y', 'z']
ORDERS = [2, 4, 6]

# block A: (gas, ep_type).  CP x PEP is rejected by input_module -- see banner.
A_KERNELS = [('CP', 'JP'), ('CP', 'KEEP'), ('CP', 'AEC'),
             ('TP', 'JP'), ('TP', 'KEEP'), ('TP', 'AEC'), ('TP', 'PEP')]

# block B: (scheme, weno_order).  TENO has no order 3; TENO-A implements only
# the weno_num=3 stencil, hence order 5 alone (input_module.f90:161).
B_SCHEMES = ([('hybrid_wenoEP',  o) for o in (3, 5, 7)] +
             [('hybrid_tenoEP',  o) for o in (5, 7)] +
             [('hybrid_tenoaEP', o) for o in (5,)])
B_ABBR = {'hybrid_wenoEP': 'weno', 'hybrid_tenoEP': 'teno', 'hybrid_tenoaEP': 'tenoa'}

C_DIFFUSION = ['laplacian', 'staggered', 'conservative']


# =========================================================================== #
#  curve enumeration -- the single source of truth for the whole campaign
# =========================================================================== #
def all_curves():
    curves = []
    for gas, ep in A_KERNELS:
        for fd in ORDERS:
            for d in DIRS:
                for m in MESHES:
                    curves.append({
                        'block': 'A', 'name': f'A_{ep}_{gas}_O{fd}_{d}_{m}',
                        'gas': gas, 'ep_type': ep, 'fd_order': fd,
                        'dir': d, 'mesh': m,
                        'ladder': LADDER['A']['wavy' if m == 'wavy' else 'cart'],
                    })
    for sch, wo in B_SCHEMES:
        for gas in ('CP', 'TP'):
            for d in DIRS:
                for m in MESHES:
                    curves.append({
                        'block': 'B', 'name': f'B_{B_ABBR[sch]}{wo}_{gas}_{d}_{m}',
                        'gas': gas, 'scheme': sch, 'weno_order': wo,
                        'dir': d, 'mesh': m,
                        'ladder': LADDER['B']['wavy' if m == 'wavy' else 'cart'],
                    })
    for diff in C_DIFFUSION:
        for fd in ORDERS:
            for m in MESHES:
                curves.append({
                    'block': 'C', 'name': f'C_{diff}_O{fd}_{m}',
                    'diffusion': diff, 'fd_order': fd, 'mesh': m,
                    'ladder': LADDER['C']['wavy' if m == 'wavy' else 'cart'],
                })
    return curves


def curves_of(block):
    return [c for c in all_curves() if c['block'] == block]


def find_curve(name):
    for c in all_curves():
        if c['name'] == name:
            return c
    raise SystemExit(f'unknown curve: {name}')


# =========================================================================== #
#  mesh
# =========================================================================== #
def grid_indices(n):
    """integer node indices 1-NG .. n+NG (ghost-inclusive)."""
    return np.arange(1-NG, n+NG+1)


def build_mesh(meshtype, nx, ny):
    """ghost-padded XX,YY of shape (nx+2NG, ny+2NG).  The Cartesian node
       convention is cell-centred, s(i) = smin + L*(i-0.5)/n, matching
       uniformGrid_f0; wavy adds the periodic Soldati warp on top."""
    ii, jj = grid_indices(nx), grid_indices(ny)
    csi = (ii - 0.5)/nx
    eta = (jj - 0.5)/ny
    CSI, ETA = np.meshgrid(csi, eta, indexing='ij')
    if meshtype == 'wavy':
        XX = XMIN + LX*CSI + WAVY_AX*np.sin(2.0*np.pi*ETA)
        YY = YMIN + LY*ETA + WAVY_AY*np.sin(4.0*np.pi*CSI)
    else:                                   # real_cart / fake_cart geometry
        XX = XMIN + LX*CSI
        YY = YMIN + LY*ETA
    return XX, YY


def write_mesh_file(path, XX, YY, nx, ny, nz):
    """URANOS curvilinear format: header 'nx ny nz', x-block then y-block,
       i fastest (Fortran order), one value per line."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'{nx} {ny} {nz}\n')
        np.savetxt(f, XX.flatten(order='F')[:, None], fmt='%25.16e')
        np.savetxt(f, YY.flatten(order='F')[:, None], fmt='%25.16e')


def d_idx(a, axis):
    out = np.zeros_like(a)
    for l in range(-3, 4):
        out += C1[l+3]*np.roll(a, -l, axis=axis)
    return out


def detJ2D_interior(XX, YY, nx, ny):
    """in-plane Jacobian d(x,y)/d(csi,eta) on the interior, by the same order-6
       central FD the solver uses.  d_idx assumes unit index spacing, so this is
       the physical area per unit index^2; the caller divides by nx*ny."""
    x_xi, x_eta = d_idx(XX, 0), d_idx(XX, 1)
    y_xi, y_eta = d_idx(YY, 0), d_idx(YY, 1)
    sl = (slice(NG, NG+nx), slice(NG, NG+ny))
    return x_xi[sl]*y_eta[sl] - x_eta[sl]*y_xi[sl]


# =========================================================================== #
#  BINARY reader
# =========================================================================== #
def read_binary(fn):
    """returns (time, nx, ny, nz, rho, u) -- phi(1)=rho, phi(2)=rho*u."""
    with open(fn, 'rb') as f:
        nx, ny, nz, it = np.fromfile(f, np.int32, 4)
        time, dt       = np.fromfile(f, np.float64, 2)
        n = int(nx)*int(ny)*int(nz)
        flat = np.fromfile(f, np.float64, 5*n)
    rho  = flat[0:n].reshape((nx, ny, nz), order='F')
    rhou = flat[n:2*n].reshape((nx, ny, nz), order='F')
    return float(time), int(nx), int(ny), int(nz), rho, rhou/rho


def final_dump(bindir):
    files = glob.glob(os.path.join(bindir, '*.bin'))
    if not files:
        raise RuntimeError(f'no .bin in {bindir}')
    best, bestt = None, -1.0
    for fn in files:
        t = read_binary(fn)[0]
        if t > bestt:
            best, bestt = fn, t
    return best


# =========================================================================== #
#  namelist template -- one template, every block
# =========================================================================== #
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
	Fourier = 0.1
	Dt   = {dt}
	bc(1) = 'periodic'
	bc(2) = 'periodic'
	bc(3) = 'periodic'
	bc(4) = 'periodic'
	bc(5) = 'periodic'
	bc(6) = 'periodic'
	Trat = 1.0
	Tref = 298.15
	species_name = "air"
	vis_flag = 1
	thermally_perfect_gas = {tpg}
	inflow_profile = 'none'
	smooth_inflow  = .false.
	turb_inflow    = .false.
	Reynolds = {re}
	Mach = {mach}
	Prandtl = 0.71
	logical_CFL = {logical_cfl}
	scheme      = '{scheme}'
	ep_type     = '{ep_type}'
	diffusion_scheme = '{diffusion}'
	kolmogorov_forcing = {forcing}
	mu_jump_ratio = 1.0
	sensor      = '{sensor}'
	ducrosToll  = 1.0
	enable_limiter = .false.
	fd_order    = {fd_order}
	weno_order  = {weno_order}
	LES         = .false.
	sgs_model   = 'WALE'
	yoshizawa   = .false.
	viscous = {viscous}
	n_step = 3
	printres = .true.
	itout = {itmax}
	StOut = 100000000
	rsOut = 100000000
	StFlg = .false.
	ic = '{ic}'
	output_file_name = '{ofn}'
	data_dir = '{data_dir}'
&end
"""


def case_dims(curve, N):
    """which index the ladder refines.

       A/B: the solution is 1D along the advection direction.  On a Cartesian
       mesh only that index is refined; on wavy the in-plane metric varies in
       both i and j, so x and y refine isotropically (z stays separable).
       C: the steady solution is 1D in y -- same argument."""
    if curve['block'] == 'C':
        base = BASE_VIS
        return (N, N, base) if curve['mesh'] == 'wavy' else (base, N, base)
    base = BASE_ADV
    d, wavy = curve['dir'], curve['mesh'] == 'wavy'
    if d == 'z':
        return base, base, N                    # z is separable on every mesh
    if wavy:
        return N, N, base
    return (N, base, base) if d == 'x' else (base, N, base)


def itmax_for(curve, N):
    if curve['block'] != 'C':
        return ITMAX_ADV
    C = STEADY_C_WAVY if curve['mesh'] == 'wavy' else STEADY_C_CART
    return max(3000, int(math.ceil(C*N*N)))


def namelist_for(curve, N, nx, ny, nz, mesh_file, data_dir, tag):
    b = curve['block']
    common = dict(
        xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX, zmin=ZMIN, zmax=ZMAX,
        mesh_file=mesh_file, itmax=itmax_for(curve, N), nx=nx, ny=ny, nz=nz,
        ofn=tag, data_dir=data_dir, fd_order=curve.get('fd_order', 6),
    )
    if b == 'C':
        common.update(
            cfl=CFL_VIS, dt='1.0E-03', logical_cfl='.true.',
            tpg=0, re=RE_VIS, mach=MACH_VIS,
            scheme='none', ep_type='JP', diffusion=curve['diffusion'],
            forcing=f'{FORCING:.6E}', sensor='density', weno_order=5,
            viscous='.true.', ic='KolmogorovFlow',
        )
    else:
        common.update(
            cfl=0.1, dt=f'{DT_ADV:.6E}', logical_cfl='.false.',
            tpg=1 if curve['gas'] == 'TP' else 0,
            re='10E+10', mach=U_INF, diffusion='staggered',
            forcing='0.0', viscous='.false.',
            ic=f"periodic_euler_{curve['dir']}",
        )
        if b == 'A':
            common.update(scheme='energy_preserving', ep_type=curve['ep_type'],
                          sensor='ducros', weno_order=5)
        else:
            # sensor='weno' forces compute_pure_weno_flag: every cell takes the
            # shock branch, so the EP half is never reached and ep_type is inert.
            common.update(scheme=curve['scheme'], ep_type='JP',
                          sensor='weno', weno_order=curve['weno_order'])
    return DAT_TEMPLATE.format(**common)


# =========================================================================== #
#  one run
# =========================================================================== #
def write_case(curve, N):
    """Mesh file (one per geometry, shared by every curve that uses it) and
       deck of one run; returns (datpath, data_dir, XX, YY, nx, ny, nz).
       Paths in the deck are relative to the repository root, where URANOS runs."""
    name = curve['name']
    tag  = f'N{N:04d}'
    nx, ny, nz = case_dims(curve, N)

    XX, YY = build_mesh(curve['mesh'], nx, ny)
    if curve['mesh'] == 'real_cart':
        mesh_file = ''
    else:
        path = os.path.join(GRID_ROOT, f"{curve['mesh']}_{nx}x{ny}x{nz}.txt")
        if not os.path.isfile(path):
            write_mesh_file(path, XX, YY, nx, ny, nz)
        mesh_file = os.path.relpath(path, REPO)

    data_dir = f"{CAMPAIGN}/runs/{BLOCK_DIR[curve['block']]}/{name}/{tag}"
    inpdir = os.path.join(INP_ROOT, BLOCK_DIR[curve['block']], name)
    os.makedirs(inpdir, exist_ok=True)
    datpath = os.path.join(inpdir, f'{tag}.dat')
    with open(datpath, 'w') as f:
        f.write(namelist_for(curve, N, nx, ny, nz, mesh_file, data_dir, f'{name}_{tag}'))
    return datpath, data_dir, XX, YY, nx, ny, nz


def run_one(curve, N, launcher, keep_binary=True):
    name = curve['name']
    tag  = f'N{N:04d}'
    datpath, data_dir, XX, YY, nx, ny, nz = write_case(curve, N)

    outroot  = os.path.join(REPO, 'DATA', data_dir)
    shutil.rmtree(outroot, ignore_errors=True)
    os.makedirs(os.path.join(outroot, 'BINARY'), exist_ok=True)

    logdir = os.path.join(LOG_ROOT, name)
    os.makedirs(logdir, exist_ok=True)

    cmd = launcher.split() + [EXE, datpath]
    t0  = datetime.datetime.now()
    log = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    wall = (datetime.datetime.now() - t0).total_seconds()
    with open(os.path.join(logdir, f'{tag}.log'), 'w') as f:
        f.write(log.stdout)
        if log.stderr:
            f.write('\n===== stderr =====\n' + log.stderr)

    rec = {'curve': name, 'N': N, 'nx': nx, 'ny': ny, 'nz': nz,
           'itmax': itmax_for(curve, N), 'wall_s': round(wall, 2)}
    if log.returncode != 0:
        rec['error'] = f'solver rc={log.returncode}'
        write_error_txt(curve, rec, outroot)
        return rec

    try:
        fn = final_dump(os.path.join(outroot, 'BINARY'))
        t, bnx, bny, bnz, rho, u = read_binary(fn)
    except Exception as exc:                     # noqa: BLE001 -- recorded, not raised
        rec['error'] = f'read: {exc}'
        write_error_txt(curve, rec, outroot)
        return rec

    # ---- exact solution at the physical node positions -------------------- #
    if curve['block'] == 'C':
        yy = YY[NG:NG+nx, NG:NG+ny][:, :, None]
        num, exact = u, A_AN*np.sin(KWAVE*yy) + np.zeros_like(u)
        h = LY/N
    else:
        d = curve['dir']
        if d == 'x':
            s = XX[NG:NG+nx, NG:NG+ny][:, :, None]
        elif d == 'y':
            s = YY[NG:NG+nx, NG:NG+ny][:, :, None]
        else:
            kk = np.arange(1, nz+1)
            s = (ZMIN + LZ*(kk-0.5)/nz)[None, None, :]
        num   = rho
        exact = 1.0 + A_WAVE*np.sin(np.pi*(s - U_INF*t)) + np.zeros_like(rho)
        h = LX/N

    # ---- volume-weighted norms ------------------------------------------- #
    detJ = detJ2D_interior(XX, YY, nx, ny)
    dV   = (detJ[:, :, None]/(nx*ny))*(LZ/nz) + np.zeros_like(num)
    err  = np.abs(num - exact)
    Vtot = np.sum(dV)
    rec.update({'t': t, 'h': h,
                'L1': float(np.sum(err*dV)/Vtot),
                'L2': float(np.sqrt(np.sum(err*err*dV)/Vtot)),
                'Linf': float(np.max(err))})
    if curve['block'] == 'C':
        sn = np.sin(KWAVE*yy) + np.zeros_like(u)
        rec['A_ratio'] = float(np.sum(u*sn*dV)/np.sum(sn*sn*dV)/A_AN)
    if not keep_binary:
        shutil.rmtree(os.path.join(outroot, 'BINARY'), ignore_errors=True)
    write_error_txt(curve, rec, outroot)
    return rec


# =========================================================================== #
#  the .txt files
# =========================================================================== #
def describe(curve):
    """the header block shared by error.txt and scaling.txt."""
    b = curve['block']
    L = [f"# curve      = {curve['name']}",
         f"# block      = {b}_{ {'A':'ep','B':'weno','C':'viscous'}[b] }",
         f"# mesh       = {curve['mesh']}"]
    if b == 'A':
        L += [f"# scheme     = energy_preserving",
              f"# ep_type    = {curve['ep_type']}",
              f"# gas        = {curve['gas']}  (thermally_perfect_gas = {1 if curve['gas']=='TP' else 0})",
              f"# fd_order   = {curve['fd_order']}",
              f"# direction  = {curve['dir']}",
              f"# case       = periodic density wave, rho = 1 + {A_WAVE} sin(pi s), u_inf = {U_INF}, t_fin = 1",
              f"# error on   = density",
              f"# expected   = {curve['fd_order']}"]
    elif b == 'B':
        L += [f"# scheme     = {curve['scheme']}  (sensor = 'weno': pure reconstruction, no EP blending)",
              f"# weno_order = {curve['weno_order']}",
              f"# gas        = {curve['gas']}  (thermally_perfect_gas = {1 if curve['gas']=='TP' else 0})",
              f"# direction  = {curve['dir']}",
              f"# case       = periodic density wave, rho = 1 + {A_WAVE} sin(pi s), u_inf = {U_INF}, t_fin = 1",
              f"# error on   = density",
              f"# expected   = {curve['weno_order']}  (WENO3 drops to ~2 at critical points, by construction)"]
    else:
        L += [f"# scheme     = none  (pure viscous)",
              f"# diffusion  = {curve['diffusion']}",
              f"# fd_order   = {curve['fd_order']}",
              f"# case       = steady Kolmogorov flow, u* = {A_AN:.6e} sin({KWAVE:.6f} y), Re = {RE_VIS}, forcing = {FORCING}",
              f"# error on   = streamwise velocity",
              f"# expected   = {curve['fd_order']}  (staggered saturates at 2 on wavy: variable face coefficient)"]
    L += ["# norm       = volume-weighted Lp with the physical cell volume as weight",
          f"# generated  = {datetime.datetime.now().isoformat(timespec='seconds')}"]
    return L


def write_error_txt(curve, rec, outroot):
    """self-contained record of a single run, so the point can be recomputed or
       relaunched without replaying the whole ladder."""
    lines = describe(curve)
    lines += [f"# nx,ny,nz   = {rec.get('nx')},{rec.get('ny')},{rec.get('nz')}",
              f"# itmax      = {rec.get('itmax')}",
              f"# wall_s     = {rec.get('wall_s')}"]
    if 'error' in rec:
        lines += ['#', f"# FAILED: {rec['error']}", '#',
                  '#    N             h               L1               L2             Linf',
                  f"{rec['N']:6d}" + '           nan              nan              nan              nan']
    else:
        lines += ['#',
                  '#    N             h               L1               L2             Linf',
                  f"{rec['N']:6d}  {rec['h']:.8e}  {rec['L1']:.8e}  {rec['L2']:.8e}  {rec['Linf']:.8e}"]
    with open(os.path.join(outroot, 'error.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')


def fit_order(hs, errs, floor=1e-12):
    """least squares slope over the points above the round-off floor.  Points
       that have hit double-precision saturation would otherwise flatten the fit
       and read as a loss of order."""
    hs, errs = np.asarray(hs, float), np.asarray(errs, float)
    m = np.isfinite(errs) & (errs > floor)
    if m.sum() < 2:
        m = np.isfinite(errs)
    if m.sum() < 2:
        return float('nan')
    return float(np.polyfit(np.log(hs[m]), np.log(errs[m]), 1)[0])


def write_scaling_txt(curve, recs, curvedir):
    """the plottable file: one row per N, local order against the previous row,
       least-squares order in the header.

          gnuplot> set log xy
          gnuplot> plot 'scaling.txt' u 2:4 w lp    # L2 vs h
    """
    ok = [r for r in recs if 'error' not in r]
    hs = [r['h'] for r in ok]
    fits = {k: fit_order(hs, [r[k] for r in ok]) for k in ('L1', 'L2', 'Linf')}
    lines = describe(curve)
    lines += ['#',
              f"# fit (least squares, points above 1e-12):"
              f"  L1 = {fits['L1']:6.3f}   L2 = {fits['L2']:6.3f}   Linf = {fits['Linf']:6.3f}",
              f"# runs       = {len(ok)}/{len(recs)} completed",
              '#',
              '#    N             h               L1               L2             Linf'
              '      pL1     pL2   pLinf']
    prev = None
    for r in recs:
        if 'error' in r:
            lines.append(f"{r['N']:6d}" + '           nan              nan              nan'
                                          '              nan      nan     nan     nan')
            continue
        cols = f"{r['N']:6d}  {r['h']:.8e}  {r['L1']:.8e}  {r['L2']:.8e}  {r['Linf']:.8e}"
        if prev is None:
            cols += '      nan     nan     nan'
        else:
            lr = math.log(prev['h']/r['h'])
            for k in ('L1', 'L2', 'Linf'):
                p = math.log(prev[k]/r[k])/lr if (r[k] > 0 and prev[k] > 0) else float('nan')
                cols += f'  {p:7.3f}'
        lines.append(cols)
        prev = r
    with open(os.path.join(curvedir, 'scaling.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return fits


# =========================================================================== #
#  one curve = one SLURM array task
# =========================================================================== #
def run_curve(curve, launcher, keep_binary=True):
    name = curve['name']
    curvedir = os.path.join(DATA_ROOT, BLOCK_DIR[curve['block']], name)
    os.makedirs(curvedir, exist_ok=True)
    recs = []
    for N in curve['ladder']:
        rec = run_one(curve, N, launcher, keep_binary)
        recs.append(rec)
        if 'error' in rec:
            print(f"  {name:38s} N={N:4d}  FAILED  {rec['error']}", flush=True)
        else:
            print(f"  {name:38s} N={N:4d}  L1={rec['L1']:.4e} L2={rec['L2']:.4e} "
                  f"Linf={rec['Linf']:.4e}  {rec['wall_s']:6.1f}s", flush=True)
    fits = write_scaling_txt(curve, recs, curvedir)
    os.makedirs(RES_ROOT, exist_ok=True)
    with open(os.path.join(RES_ROOT, f'{name}.json'), 'w') as f:
        json.dump({'curve': {k: v for k, v in curve.items()},
                   'fit': fits, 'runs': recs}, f, indent=1)
    print(f"  {name:38s} FIT  L1={fits['L1']:.3f} L2={fits['L2']:.3f} "
          f"Linf={fits['Linf']:.3f}", flush=True)
    return recs, fits


# =========================================================================== #
def cmd_list(args):
    for c in all_curves():
        if args.block and c['block'] != args.block:
            continue
        print(c['name'])


def cmd_plan(args):
    tot_runs = 0
    for b in ('A', 'B', 'C'):
        cs = curves_of(b)
        n = sum(len(c['ladder']) for c in cs)
        tot_runs += n
        print(f'  block {b} ({BLOCK_DIR[b]:10s}): {len(cs):4d} curves  {n:5d} runs')
    print(f'  {"total":24s}: {len(all_curves()):4d} curves  {tot_runs:5d} runs')


def cmd_prepare(args):
    n = 0
    for c in all_curves():
        for N in c['ladder']:
            write_case(c, N); n += 1
    print(f'  {n} decks in {INP_ROOT}, {len(os.listdir(GRID_ROOT))} mesh files in {GRID_ROOT}')


def cmd_run(args):
    if args.curve:
        curve = find_curve(args.curve)
    else:
        cs = curves_of(args.block)
        if not 0 <= args.index < len(cs):
            raise SystemExit(f'index {args.index} out of range for block {args.block} (0..{len(cs)-1})')
        curve = cs[args.index]
    if not os.path.isfile(EXE):
        raise SystemExit(f'missing executable: {EXE}')
    print(f"[{datetime.datetime.now().isoformat(timespec='seconds')}] curve {curve['name']}", flush=True)
    run_curve(curve, args.launcher, keep_binary=not args.drop_binary)


def cmd_gate(args):
    """fake_cart must reproduce real_cart to round-off: the curvilinear metric
       path on a Cartesian geometry is the degenerate case of the general one.
       Run this on the same executable that will run the campaign -- GPU and CPU
       differ at ~8e-15, so the comparison is only meaningful within a backend."""
    stems = ['A_JP_TP_O6_x', 'A_PEP_TP_O6_x', 'B_weno5_CP_x', 'B_teno7_TP_x',
             'C_laplacian_O6', 'C_staggered_O6', 'C_conservative_O6']
    N, bad = args.N, 0
    for stem in stems:
        vals = {m: run_one(find_curve(f'{stem}_{m}'), N, args.launcher)
                for m in ('real_cart', 'fake_cart')}
        a, b = vals['real_cart'], vals['fake_cart']
        if 'error' in a or 'error' in b:
            print(f'  {stem:22s} FAILED to run'); bad += 1; continue
        da = max(abs(a[k]-b[k]) for k in ('L1', 'L2', 'Linf'))
        dr = max(abs(a[k]-b[k])/max(abs(a[k]), 1e-300) for k in ('L1', 'L2', 'Linf'))
        # The two paths are the same discretization evaluated in a different
        # order, so they agree to floating-point associativity, not bitwise.
        # The field is O(1), hence an absolute drift at round-off is the right
        # criterion; a relative one divides by an error that is itself ~1e-11
        # on the finer schemes and turns 1e-15 of noise into a 1e-3 "mismatch".
        ok = (da < args.atol) or (dr < 1e-10)
        bad += 0 if ok else 1
        print(f'  {stem:22s} N={N}  L2 real={a["L2"]:.6e} fake={b["L2"]:.6e}  '
              f'abs={da:.2e} rel={dr:.2e}   {"OK" if ok else "MISMATCH"}')
    print(f'  gate: {len(stems)-bad}/{len(stems)} passed')
    raise SystemExit(0 if bad == 0 else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('list');  p.add_argument('--block', choices=list(BLOCK_DIR))
    p.set_defaults(func=cmd_list)

    p = sub.add_parser('plan');  p.set_defaults(func=cmd_plan)

    p = sub.add_parser('prepare'); p.set_defaults(func=cmd_prepare)

    p = sub.add_parser('run')
    p.add_argument('--curve')
    p.add_argument('--block', choices=list(BLOCK_DIR))
    p.add_argument('--index', type=int)
    p.add_argument('--launcher', default='mpirun -np 1')
    p.add_argument('--drop-binary', action='store_true', help='delete BINARY after reducing')
    p.set_defaults(func=cmd_run)

    p = sub.add_parser('gate')
    p.add_argument('--N', type=int, default=32)
    p.add_argument('--atol', type=float, default=1e-13,
                   help='absolute agreement required between the real_cart and fake_cart norms')
    p.add_argument('--launcher', default='mpirun -np 1')
    p.set_defaults(func=cmd_gate)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
