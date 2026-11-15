#!/usr/bin/env python3
"""
view_mesh.py : visualize + sanity-check a GENERATED URANOS curvilinear mesh .txt
               (the file written by ramp_mesh_gen.exe / gen_martin_mesh.py).

Unlike martin_mesh_explorer.py (which rebuilds the mesh from parameters), this reads
the ACTUAL file on disk -- so it shows exactly what URANOS will read, ghosts included.

USAGE
-----
    python3 view_mesh.py                                   # default 6delta grid
    python3 view_mesh.py ../../../../grids/les_RAMP/mesh_ramp_martin_6delta.txt
    python3 view_mesh.py <file.txt> --ghosts               # also draw the GN ghost layers
    python3 view_mesh.py <file.txt> --H 6                  # check the wall-normal gap vs H

Reports (and prints in the title):
  * inlet / outlet wall->top column length  (should both equal H; a mismatch or a
    value != H means the wall-normal grid overshoot -- e.g. a hardcoded NyWr; use NyWr=0)
  * min cell Jacobian + folded-cell count
  * TOP-edge inflow fraction for a horizontal free-stream (un<0 -> nscbc_outflow fails)
  * first wall dy and dy+ (needs --retau)
"""
import sys, os
import numpy as np

GN = 7
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "../../../../grids/les_RAMP/mesh_ramp_martin_12delta.txt")


def read_mesh(path):
    """Read URANOS curvilinear PLOT3D block: header 'nx ny 1' then X block then Y
    block, each (nx+2*GN)*(ny+2*GN) reals, i fastest then j.  Returns full arrays
    (with ghosts) shaped [i,j] and the interior counts nx,ny."""
    tok = open(path).read().split()
    nx, ny = int(tok[0]), int(tok[1])
    nxt, nyt = nx + 2 * GN, ny + 2 * GN
    v = np.asarray(tok[3:3 + 2 * nxt * nyt], dtype=float)
    X = v[:nxt * nyt].reshape(nyt, nxt).T
    Y = v[nxt * nyt:2 * nxt * nyt].reshape(nyt, nxt).T
    return X, Y, nx, ny


def col_len(cx, cy):
    return float(np.sum(np.hypot(np.diff(cx), np.diff(cy))))


def min_jacobian(X, Y):
    xxi = 0.5 * (X[2:, 1:-1] - X[:-2, 1:-1]); yxi = 0.5 * (Y[2:, 1:-1] - Y[:-2, 1:-1])
    xet = 0.5 * (X[1:-1, 2:] - X[1:-1, :-2]); yet = 0.5 * (Y[1:-1, 2:] - Y[1:-1, :-2])
    J = xxi * yet - xet * yxi
    return float(J.min()), int((J <= 0).sum())


def top_inflow_frac(xt, yt):
    """% of the top edge where un<0 vs a horizontal free-stream u=(1,0)."""
    dx, dy = np.gradient(xt), np.gradient(yt)
    un = -dy / np.maximum(np.hypot(dx, dy), 1e-30)
    return 100.0 * float(np.mean(un < -1e-3)), un


def arg(name, default, cast=str):
    a = sys.argv[1:]
    return cast(a[a.index(name) + 1]) if name in a and a.index(name) + 1 < len(a) else default


def main():
    a = sys.argv[1:]
    valopts = ('--H', '--retau', '--save', '--at')
    path = next((s for s in a if not s.startswith('--') and
                 (a.index(s) == 0 or a[a.index(s) - 1] not in valopts)), None)
    path = path or DEFAULT
    if not os.path.exists(path):
        sys.exit(f"  mesh file not found: {path}")
    show_ghosts = '--ghosts' in a
    H = float(arg('--H', 'nan', str))
    retau = float(arg('--retau', '180', str))
    save = arg('--save', None)

    X, Y, nx, ny = read_mesh(path)
    Xi, Yi = X[GN:GN + nx, GN:GN + ny], Y[GN:GN + nx, GN:GN + ny]

    gi = col_len(Xi[0], Yi[0])                  # inlet wall->top column length
    go = col_len(Xi[-1], Yi[-1])                # outlet wall->top column length
    jmin, nbad = min_jacobian(X, Y)             # full field incl. ghosts
    frac, un = top_inflow_frac(Xi[:, -1], Yi[:, -1])
    dy1 = np.hypot(Xi[0, 1] - Xi[0, 0], Yi[0, 1] - Yi[0, 0])

    print(f"  file        : {path}")
    print(f"  interior    : {nx} x {ny}   (block {nx+2*GN} x {ny+2*GN} incl. {GN} ghosts/side)")
    print(f"  x range     : [{Xi.min():.3f}, {Xi.max():.3f}]   y range: [{Yi.min():.3f}, {Yi.max():.3f}]")
    print(f"  INLET  gap  : {gi:.4f}  (wall {Yi[0,0]:.3f} -> top {Yi[0,-1]:.3f})")
    print(f"  OUTLET gap  : {go:.4f}  (wall {Yi[-1,0]:.3f} -> top {Yi[-1,-1]:.3f})")
    if not np.isnan(H):
        f_i, f_o = gi / H - 1, go / H - 1
        flag = "  <-- OVERSHOOT (use NyWr=0 auto-fit!)" if max(abs(f_i), abs(f_o)) > 0.02 else "  OK"
        print(f"  vs H={H:g}    : inlet {f_i*100:+.1f}%  outlet {f_o*100:+.1f}%{flag}")
    print(f"  in/out skew : {abs(gi-go)/max(gi,go)*100:.2f}%  (should be ~0)")
    print(f"  min Jacobian: {jmin:.3e}   folded cells: {nbad}   {'FOLDED!' if nbad else 'valid'}")
    print(f"  1st wall dy : {dy1:.5f}   (dy+ ~ {dy1*retau:.2f} at ReTau={retau:g})")
    # dy+ readout along the inlet column (eta grid) at requested heights
    ats = arg('--at', None)
    if ats:
        cy = Yi[0, :]                                   # inlet column y = eta (wall at 0)
        ymid = 0.5 * (cy[1:] + cy[:-1])
        dyp = np.diff(cy) * retau
        print("  dy+ at height (inlet column):")
        for Y in [float(s) for s in ats.split(',')]:
            print(f"     y = {Y:5.2f}  ->  dy+ = {np.interp(Y, ymid, dyp):6.2f}")
        for thr in (5.0, 10.0):
            ok = ymid[dyp <= thr]
            print(f"     dy+ <= {thr:4.0f} up to y = {(ok.max() if ok.size else 0.0):.2f}")
    print(f"  TOP inflow  : {frac:.0f}% of top edge has un<0 (horizontal free-stream)")

    import matplotlib
    matplotlib.use('Agg' if save else 'TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    fig, ax = plt.subplots(figsize=(15, 6))
    Xp, Yp = (X, Y) if show_ghosts else (Xi, Yi)
    step = max(1, Xp.shape[0] // 170)
    stey = max(1, Xp.shape[1] // 90)
    seg = [np.column_stack([Xp[i, :], Yp[i, :]]) for i in range(0, Xp.shape[0], step)]
    seg += [np.column_stack([Xp[:, j], Yp[:, j]]) for j in range(0, Xp.shape[1], stey)]
    ax.add_collection(LineCollection(seg, colors='k', lw=0.2, alpha=0.5))
    ax.plot(Xi[:, 0], Yi[:, 0], 'r-', lw=1.6, label='wall')
    ax.plot(Xi[0, :], Yi[0, :], 'g-', lw=1.6, label='inlet')
    ax.plot(Xi[-1, :], Yi[-1, :], 'm-', lw=1.6, label='outlet')
    ax.plot(Xi[:, -1], Yi[:, -1], 'b-', lw=1.3, label='top (un>=0)')
    infl = un < -1e-3
    if infl.any():
        ax.plot(Xi[infl, -1], Yi[infl, -1], 'r.', ms=4, label='top un<0 (nscbc fails)')
    if show_ghosts:
        ax.plot(X[:, GN + ny - 1], Y[:, GN + ny - 1], color='orange', lw=0.6, alpha=0.5)
    ax.set_aspect('equal'); ax.legend(loc='upper left', fontsize=8)
    ax.set_xlabel(r'$x/\delta$'); ax.set_ylabel(r'$y/\delta$')
    ax.set_title(f"{os.path.basename(path)}  {nx}x{ny}   inlet_gap={gi:.3f} outlet_gap={go:.3f}   "
                 f"minJ={jmin:.1e} folded={nbad}   TOPinflow={frac:.0f}%",
                 color=('red' if (nbad or frac > 0) else 'black'), fontsize=10)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=140); print(f"  saved {save}")
    else:
        png = path[:-4] + "_view.png"
        fig.savefig(png, dpi=140); print(f"  saved {png}")
        try:
            plt.show()
        except Exception:
            pass


if __name__ == '__main__':
    main()
