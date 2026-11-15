#!/usr/bin/env python3
"""
RK order verification, part 2 -- the scheme INSIDE URANOS.

Temporal self-convergence: the grid, the scheme and every spatial operator are
held FIXED and only Dt is refined, so the spatial discretisation error cancels
exactly and what is left is purely the Runge-Kutta error.

    err(Dt) = || phi(T; Dt) - phi(T; Dt_ref) ||_2 ,   Dt_ref = Dt_coarse / 32

All runs integrate to the SAME final time T = N*Dt with logical_CFL = .false.
(fixed Dt), so the comparison is at identical physical time.

Cases
-----
  dwave   passive density wave, 2D, periodic, INVISCID   -> autonomous, linear
  tgv     Taylor-Green vortex, 3D, periodic, VISCOUS     -> autonomous, nonlinear
  stokes  Stokes 2nd problem, oscillating wall           -> NON-AUTONOMOUS
          (bc_module.f90:1464  vel = u_inf*cos(0.25*pi*time))

Expected: order 3 on the two autonomous cases; the non-autonomous case exposes
the stage-time bug (see rk_tableau_and_ode.py).

Usage:  python3 uranos_dt_sweep.py [--exe ./Uranos.exe] [--tag baseline] [case ...]
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))   # the URANOS repo root

# ---------------------------------------------------------------------------
# case definitions:  a namelist template + the Dt sweep
# ---------------------------------------------------------------------------
COMMON = """&input_list
    xmin = {xmin}
    xmax = {xmax}
    ymin = {ymin}
    ymax = {ymax}
    zmin = {zmin}
    zmax = {zmax}
    gridpoint_x = 'uniform'
    gridpoint_y = 'uniform'
    gridpoint_z = 'uniform'
    dims = {dims}
    nx = {nx}
    ny = {ny}
    nz = {nz}
    cart_dims(1) = 0
    cart_dims(2) = 0
    cart_dims(3) = 0

    ! ---- fixed time step: this is what the sweep refines ----
    logical_CFL = .false.
    Dt   = {dt:.16e}
    CFL  = 0.5
    tmax = 1.e9
    itmax = {itmax}
    n_step = 3

    bc(1) = '{bc1}'
    bc(2) = '{bc2}'
    bc(3) = '{bc3}'
    bc(4) = '{bc4}'
    bc(5) = '{bc5}'
    bc(6) = '{bc6}'
    Trat = 1.0

    inflow_profile = 'none'
    smooth_inflow  = .false.
    turb_inflow    = .false.

    Reynolds = {Re}
    Mach     = {Ma}
    Prandtl  = 0.71
    viscous  = {viscous}
    vis_flag = 0
    thermally_perfect_gas = 0

    scheme     = 'energy_preserving'
    ep_type    = 'KEEP'
    sensor     = 'density'
    fd_order   = 6
    weno_order = 5
    diffusion_scheme = 'laplacian'
    LES        = .false.
    sgs_model  = 'WALE'
    yoshizawa  = .false.

    ic = '{ic}'
    itout  = {itmax}
    StOut  = 1000000
    StFlg  = .false.
    rsOut  = 1000000
    printres = .true.
    output_file_name = 'fld'
    data_dir = '{data_dir}'
    mpi_opt_level = 2
    cuda_aware = .false.

&end
"""

CASES = {
    # ---- autonomous, linear: a passive density wave on a periodic box ------
    "dwave": dict(
        ic="density_wave", dims=2, nx=32, ny=32, nz=0,
        xmin=-1.0, xmax=1.0, ymin=-1.0, ymax=1.0, zmin=0.0, zmax=0.0,
        bc1="periodic", bc2="periodic", bc3="periodic",
        bc4="periodic", bc5="periodic", bc6="periodic",
        Re=1000.0, Ma=0.2, viscous=".false.",
        dt0=8.0e-3, n0=20, np=1,
        blurb="passive density wave, 2D periodic, INVISCID (autonomous, linear RHS)",
    ),
    # ---- autonomous, nonlinear + viscous: Taylor-Green ---------------------
    "tgv": dict(
        ic="TGV", dims=3, nx=32, ny=32, nz=32,
        xmin=0.0, xmax=6.283185307179586, ymin=0.0, ymax=6.283185307179586,
        zmin=0.0, zmax=6.283185307179586,
        bc1="periodic", bc2="periodic", bc3="periodic",
        bc4="periodic", bc5="periodic", bc6="periodic",
        Re=160.0, Ma=0.2, viscous=".true.",
        dt0=4.0e-3, n0=20, np=4,
        blurb="Taylor-Green vortex, 3D periodic, VISCOUS (autonomous, fully nonlinear RHS)",
    ),
    # ---- NON-autonomous: oscillating wall (Stokes 2nd problem) -------------
    "stokes": dict(
        ic="II_stokes_problem_x", dims=2, nx=16, ny=96, nz=0,
        xmin=0.0, xmax=1.0, ymin=0.0, ymax=2.0, zmin=0.0, zmax=0.0,
        bc1="periodic", bc2="periodic", bc3="isothermal_wall",
        bc4="neumann", bc5="periodic", bc6="periodic",
        Re=1.0, Ma=0.01, viscous=".true.",
        dt0=2.5e-3, n0=800, np=1,
        blurb="Stokes 2nd problem, wall velocity u_w = u_inf*cos(0.25*pi*t)  "
              "(NON-AUTONOMOUS RHS: bc_module.f90:1464)",
    ),
}

REFINEMENTS = [1, 2, 4, 8]     # Dt0/1, Dt0/2, Dt0/4, Dt0/8
REF_REFINE = 32                # the reference run: Dt0/32


def read_bin(path):
    """[nx ny nz it]_i4 + [time dt]_f8 + phi(1:5) as 5 Fortran-ordered blocks."""
    with open(path, "rb") as f:
        nx, ny, nz, it = np.fromfile(f, dtype=np.int32, count=4)
        time, dt = np.fromfile(f, dtype=np.float64, count=2)
        n = int(nx) * int(ny) * int(nz)
        raw = np.fromfile(f, dtype=np.float64, count=5 * n)
    if raw.size != 5 * n:
        raise RuntimeError(f"{path}: short read ({raw.size} != {5*n})")
    phi = raw.reshape((int(nx), int(ny), int(nz), 5), order="F")
    return dict(nx=int(nx), ny=int(ny), nz=int(nz), it=int(it),
                time=float(time), dt=float(dt), phi=phi)


def run_one(case, name, exe, refine, tag):
    cfg = CASES[case]
    dt = cfg["dt0"] / refine
    itmax = cfg["n0"] * refine
    data_dir = f"rkv_{tag}_{case}_r{refine}"
    deck = os.path.join(HERE, "decks", f"{data_dir}.dat")
    os.makedirs(os.path.dirname(deck), exist_ok=True)
    with open(deck, "w") as f:
        f.write(COMMON.format(dt=dt, itmax=itmax, data_dir=data_dir, **{
            k: v for k, v in cfg.items() if k not in ("dt0", "n0", "np", "blurb")
        }))
    out_dir = os.path.join(ROOT, "DATA", data_dir)
    shutil.rmtree(out_dir, ignore_errors=True)

    cmd = ["mpirun", "-np", str(cfg["np"]), exe, deck]
    log = os.path.join(HERE, "runs", f"{data_dir}.log")
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "w") as lf:
        r = subprocess.run(cmd, cwd=ROOT, stdout=lf, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f"run failed (rc={r.returncode}): {' '.join(cmd)}  -- see {log}")

    fbin = os.path.join(out_dir, "BINARY", f"fld{itmax:07d}.bin")
    if not os.path.exists(fbin):
        raise RuntimeError(f"missing output {fbin} -- see {log}")
    d = read_bin(fbin)
    got_t, want_t = d["time"], dt * itmax
    if abs(got_t - want_t) > 1e-10 * max(1.0, want_t):
        raise RuntimeError(f"final time mismatch: {got_t} != {want_t}")
    return d, dt


def sweep(case, exe, tag):
    cfg = CASES[case]
    print(f"\n{'='*78}\n{case.upper()}  --  {cfg['blurb']}\n{'='*78}")
    print(f"  grid {cfg['nx']}x{cfg['ny']}x{max(cfg['nz'],1)}   "
          f"final time T = {cfg['dt0']*cfg['n0']:.6g}   exe = {exe}")

    ref, dtref = run_one(case, case, exe, REF_REFINE, tag)
    print(f"  reference: Dt = {dtref:.4e}  ({cfg['n0']*REF_REFINE} steps)")

    rows = []
    for r in REFINEMENTS:
        d, dt = run_one(case, case, exe, r, tag)
        diff = d["phi"] - ref["phi"]
        # normalised L2 over the 5 conservative variables
        err = np.sqrt(np.mean(diff**2))
        scale = np.sqrt(np.mean(ref["phi"]**2))
        rows.append((dt, err, err / scale))

    print(f"\n  {'Dt':>12} {'L2 err':>14} {'rel L2':>12} {'order':>8}")
    orders = []
    for i, (dt, err, rel) in enumerate(rows):
        if i == 0:
            print(f"  {dt:>12.4e} {err:>14.6e} {rel:>12.4e} {'--':>8}")
        else:
            p = np.log(rows[i-1][1] / err) / np.log(rows[i-1][0] / dt)
            orders.append(p)
            print(f"  {dt:>12.4e} {err:>14.6e} {rel:>12.4e} {p:>8.3f}")
    obs = float(np.mean(orders[-2:]))
    print(f"\n  ==> observed temporal order = {obs:.3f}")
    return dict(case=case, rows=rows, order=obs, blurb=cfg["blurb"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="./Uranos.exe")
    ap.add_argument("--tag", default="baseline")
    ap.add_argument("cases", nargs="*", default=list(CASES))
    a = ap.parse_args()

    exe = a.exe if os.path.isabs(a.exe) else os.path.join(ROOT, a.exe)
    if not os.path.exists(exe):
        sys.exit(f"executable not found: {exe}")

    out = {}
    for c in a.cases:
        out[c] = sweep(c, exe, a.tag)

    print(f"\n{'='*78}\nSUMMARY  (tag = {a.tag})\n{'='*78}")
    for c, r in out.items():
        print(f"  {c:<8} observed temporal order = {r['order']:.3f}   -- {r['blurb']}")

    np.save(os.path.join(HERE, "runs", f"orders_{a.tag}.npy"),
            {c: dict(rows=r["rows"], order=r["order"]) for c, r in out.items()},
            allow_pickle=True)
