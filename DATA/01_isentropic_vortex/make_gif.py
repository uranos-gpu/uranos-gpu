#!/usr/bin/env python3
# ===========================================================================
#  Animate the isentropic vortex advecting on the WAVY curvilinear mesh.
#
#  Reads the URANOS BINARY dumps written by run_vortex.py and the curvilinear
#  mesh file, builds the z-vorticity by the SAME order-6 central FD + curvilinear
#  chain rule the solver uses, and renders it ON the physical (wavy) grid:
#    * vortex_wavy_<dir>.gif   animated z-vorticity over one domain crossing
#    * vortex_wavy_<dir>.png   4-frame storyboard (t = 0, 1/3, 2/3, end)
#
#  The colour scale is fixed from the initial frame: an energy-preserving scheme
#  convects the vortex without numerical dissipation, so the core amplitude stays
#  put and the blob simply translates through the curved cells and loops back.
#
#  Helper module, imported by the run_*.py / plot_*.py scripts of this folder.
# ===========================================================================
import os, sys, glob, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# order-6 central first-derivative stencil (l = -3..3), unit index spacing
C1 = np.array([-1.0/60, 3.0/20, -3.0/4, 0.0, 3.0/4, -3.0/20, 1.0/60])
NG = 7                      # ghost layers carried by the mesh file (GN+mesh_gh)

HERE   = os.path.dirname(os.path.abspath(__file__))
REPO   = os.path.abspath(os.path.join(HERE, '..', '..'))
OUTDIR = os.path.join(HERE, 'figs')


# --------------------------------------------------------------------------- #
def read_binary(fn):
    """URANOS BINARY: header [nx,ny,nz,it]_int32 + [time,dt]_f64, then
       phi(1:5)=(rho,rho u,rho v,rho w,rho E) global interior, Fortran order."""
    with open(fn, 'rb') as f:
        nx, ny, nz, it = np.fromfile(f, np.int32,  4)
        time, dt       = np.fromfile(f, np.float64, 2)
        n   = int(nx)*int(ny)*int(nz)
        flat = np.fromfile(f, np.float64, 5*n)
    phi = np.empty((5, nx, ny, nz))
    for v in range(5):
        phi[v] = flat[v*n:(v+1)*n].reshape((nx, ny, nz), order='F')
    return float(time), int(nx), int(ny), int(nz), phi

def read_mesh(fn):
    """Curvilinear mesh: header 'nx ny nz', x-block then y-block (i fastest,
       NG ghosts per side).  Returns padded xx,yy of shape (nx+2NG, ny+2NG)."""
    with open(fn) as f:
        nx, ny, nz = (int(t) for t in f.readline().split())
        vals = np.array([float(t) for t in f.read().split()])
    npad = (nx+2*NG)*(ny+2*NG)
    xx = vals[:npad].reshape((nx+2*NG, ny+2*NG), order='F')
    yy = vals[npad:2*npad].reshape((nx+2*NG, ny+2*NG), order='F')
    return nx, ny, nz, xx, yy

def d_idx(a, axis):
    """order-6 central derivative w.r.t. the (unit) index along `axis` (np.roll
       wrap — valid for the ghost-padded mesh and the periodic interior field)."""
    out = np.zeros_like(a)
    for l in range(-3, 4):
        out += C1[l+3]*np.roll(a, -l, axis=axis)
    return out

def mesh_metrics(xx, yy, nx, ny):
    """interior (nx,ny) inverse in-plane metrics by central FD of xx,yy."""
    x_xi, x_eta = d_idx(xx, 0), d_idx(xx, 1)
    y_xi, y_eta = d_idx(yy, 0), d_idx(yy, 1)
    sl = (slice(NG, NG+nx), slice(NG, NG+ny))
    x_xi, x_eta = x_xi[sl], x_eta[sl]
    y_xi, y_eta = y_xi[sl], y_eta[sl]
    detJ = x_xi*y_eta - x_eta*y_xi
    return (y_eta/detJ, -x_eta/detJ, -y_xi/detJ, x_xi/detJ)   # dxi_dx,dxi_dy,deta_dx,deta_dy

def omega_z_slice(phi, metr, k):
    """z-vorticity on the k-th z-slice via the curvilinear chain rule (periodic)."""
    dxi_dx, dxi_dy, deta_dx, deta_dy = metr
    rho = phi[0, :, :, k]; u = phi[1, :, :, k]/rho; v = phi[2, :, :, k]/rho
    u_xi, u_eta = d_idx(u, 0), d_idx(u, 1)
    v_xi, v_eta = d_idx(v, 0), d_idx(v, 1)
    vx = v_xi*dxi_dx + v_eta*deta_dx
    uy = u_xi*dxi_dy + u_eta*deta_dy
    return vx - uy


# --------------------------------------------------------------------------- #
def load_frames(bindir, metr, nx, ny, nz):
    """all dumps, sorted by physical time -> (times, list of omega_z fields)."""
    files = glob.glob(os.path.join(bindir, '*.bin'))
    if not files:
        sys.exit(f' no .bin files in {bindir} — run run_vortex.py first')
    stamped = sorted((read_binary(fn)[0], fn) for fn in files)
    ksl = nz//2
    times, fields = [], []
    for tm, fn in stamped:
        _, _, _, _, phi = read_binary(fn)
        times.append(tm)
        fields.append(omega_z_slice(phi, metr, ksl))
    return np.array(times), fields


def draw(ax, xi, yi, field, vmax, nx, ny, title):
    ax.clear()
    pc = ax.pcolormesh(xi, yi, field, cmap='RdBu_r', shading='gouraud',
                       vmin=-vmax, vmax=vmax)
    step = max(1, nx//40)
    for i in range(0, nx, step): ax.plot(xi[i, :], yi[i, :], 'k-', lw=0.15, alpha=0.30)
    for j in range(0, ny, step): ax.plot(xi[:, j], yi[:, j], 'k-', lw=0.15, alpha=0.30)
    ax.set_aspect('equal'); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_title(title)
    return pc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, choices=['x', 'y'])
    ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--frames', type=int, default=80, help='max frames in the GIF')
    args = ap.parse_args()

    tag     = f'vortex_wavy_{args.dir}'
    bindir  = os.path.join(REPO, 'DATA', tag, 'BINARY')
    meshf   = os.path.join(HERE, 'mesh_wavy.txt')
    os.makedirs(OUTDIR, exist_ok=True)

    nx, ny, nz, xx, yy = read_mesh(meshf)
    metr = mesh_metrics(xx, yy, nx, ny)
    xi = xx[NG:NG+nx, NG:NG+ny]; yi = yy[NG:NG+nx, NG:NG+ny]

    times, fields = load_frames(bindir, metr, nx, ny, nz)
    vmax = float(np.percentile(np.abs(fields[0]), 99.8))      # fixed scale from t=0
    print(f' [{tag}] {len(fields)} dumps, t in [{times[0]:.2f}, {times[-1]:.2f}], '
          f'omega_z scale = +/-{vmax:.2f}')

    plt.rcParams.update({'font.size': 11})
    arrow = {'x': r'$u_\infty\!\rightarrow$', 'y': r'$v_\infty\!\uparrow$'}[args.dir]

    # ---- storyboard png (t = 0, 1/3, 2/3, end) ---------------------------- #
    picks = [int(round(f*(len(fields)-1))) for f in (0.0, 1/3, 2/3, 1.0)]
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.2))
    for ax, idx in zip(axes, picks):
        pc = draw(ax, xi, yi, fields[idx], vmax, nx, ny, f't = {times[idx]:.2f}')
        fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, label=r'$\omega_z$')
    fig.suptitle(f'Isentropic vortex on the WAVY curvilinear grid — EP scheme '
                 f'(ep_id=1/JP), advection along {args.dir}  {arrow}')
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, f'{tag}.png'), dpi=120)
    plt.close(fig)
    print(f' [{tag}] wrote {tag}.png')

    # ---- animated gif ----------------------------------------------------- #
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    idxs = np.unique(np.linspace(0, len(fields)-1,
                                 min(args.frames, len(fields))).astype(int))
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    sm = ScalarMappable(norm=Normalize(-vmax, vmax), cmap='RdBu_r')
    fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label=r'$\omega_z$')
    def upd(kk):
        k = idxs[kk]
        draw(ax, xi, yi, fields[k], vmax, nx, ny,
             f'isentropic vortex on wavy grid ({args.dir} {arrow})   t = {times[k]:.2f}')
    anim = FuncAnimation(fig, upd, frames=len(idxs), interval=1000.0/args.fps)
    gif = os.path.join(OUTDIR, f'{tag}.gif')
    anim.save(gif, writer=PillowWriter(fps=args.fps))
    plt.close(fig)
    print(f' [{tag}] wrote {tag}.gif  ({len(idxs)} frames @ {args.fps} fps)')


if __name__ == '__main__':
    main()
