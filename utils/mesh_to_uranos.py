#!/usr/bin/env python3
# -----------------------------------------------------------------------------
#  Structured O-grid  ->  URANOS curvilinear mesh_file.
#
#  Reads a 2D structured O-grid (PLOT3D or gmsh) and emits the block read by
#  read_curvilinear_mesh (mesh_module.f90:457): header "nx ny nz", then all x,
#  then all y, one value per line, i fastest, over 1-ng..nx+ng and 1-ng..ny+ng
#  with ng = GN + mesh_gh = 4 + 3 = 7.
#
#    *.x    PLOT3D, single block, 2D.  The (i,j) map is explicit in the header;
#           the seam column is repeated (i = ni equals i = 1) and gets dropped.
#           PREFER THIS when both formats are available: nothing is reconstructed.
#    *.msh  gmsh 2.2 from Pointwise.  Exported as an UNstructured quad mesh -- no
#           (i,j) map -- so the map is rebuilt from the quad incidence graph alone:
#           boundary loops, then a radial march.  Works, but it is the fragile path.
#
#  Four conversions, none of them optional:
#
#   * wrap line   - see above; keeping the duplicate would put two coincident
#                   columns in the csi halo.  A .msh carries unique nodes, so
#                   there the drop is a no-op; it is detected, not assumed.
#   * handedness  - URANOS aborts on detJ2D <= 1e-12 (mesh_module.f90:604), so the
#                   wall loop must run CLOCKWISE.  The verse is DETECTED from the
#                   signed area, never assumed: of the meshes in this campaign the
#                   blunt ones are counter-clockwise and the sharp ones are already
#                   clockwise, so a hardcoded flip folds the grid silently.
#   * staggering  - URANOS is CELL-CENTRED with the physical boundary on the half
#                   node (uniformGrid_f0: x(1) = xmin + dx/2; the wall BC mirrors
#                   ghost j=0 onto node j=1, bc_module.f90:253-254).  The source
#                   NODE lines are therefore cell FACES: the grid points are their
#                   eta midpoints, and ny = nj - 1.
#   * ghosts      - xx/yy never go through mpi_share, so all 7 layers are built
#                   here: exact periodic wrap in csi, odd mirror through the wall
#                   for the first NREF rows and a smoothed ray beyond, geometric
#                   continuation at the far field.
#
#  The angle of attack is baked in by a rigid rotation about the leading edge:
#  URANOS has no AoA namelist key, and both uniform ICs put the stream along +x.
#
#  ⚠️ A SHARP trailing edge is NOT convertible on an O-grid, and it is not a
#     limitation of this script.  The first ghost row must be the exact odd
#     reflection (the wall BC mirrors onto it), but at a cusp the body is thinner
#     than the first cell -- measured 0.000245 c against 0.00073 c one cell from
#     the TE -- so that row lands past the opposite surface and detJ2D goes
#     negative.  NREF in {1,2,4,7} x NSM in {0..1000}, bounded inward geometric
#     continuation and csi position smoothing were all tried and all fail.
#     A TE about 1 % of chord thick would fit.
#
#  Usage:  python3 utils/mesh_to_uranos.py <in.x|in.msh> <out.txt> [aoa_deg]
#  (aoa_deg defaults to 0: the geometry is left as delivered)
# -----------------------------------------------------------------------------
import numpy as np, os, sys, re
from collections import defaultdict

NG   = 7     # GN (mpi_module.f90:39) + mesh_gh (mesh_module.f90:68)
GN   = 4     # halo over which URANOS evaluates detJ2D
NREF = 2     # wall-ghost rows built by EXACT reflection (the rest ride a ray)
NSM  = 20    # smoothing passes on the ghost launch angle; 0 folds a sharp corner
AOA  = 0.0  # deg, baked into the mesh; 0 leaves the geometry as delivered

def read_mesh(path):
    """Dispatch on the extension; both branches return node arrays (ni, nj)."""
    e = os.path.splitext(path)[1].lower()
    if e == '.x':
        return read_plot3d(path)
    if e == '.msh':
        return structured_map(*read_msh(path))
    raise SystemExit('%s: unknown mesh format (expected .x or .msh)' % path)


def read_msh(path):
    """Nodes and elements of a gmsh 2.2 ASCII file (quads = type 3, bars = type 1)."""
    L = open(path).read().splitlines()
    i = L.index('$Nodes'); n = int(L[i + 1])
    xy = np.zeros((n + 1, 2))
    for k in range(n):
        p = L[i + 2 + k].split(); xy[int(p[0])] = (float(p[1]), float(p[2]))
    j = L.index('$Elements'); ne = int(L[j + 1])
    quads, bars = [], []
    for k in range(ne):
        p = L[j + 2 + k].split(); t = int(p[1]); nt = int(p[2])
        c = [int(v) for v in p[3 + nt:]]
        (quads if t == 3 else bars).append(c)
    return xy, np.array(quads), bars


def structured_map(xy, quads, bars):
    """Rebuild the (i,j) O-grid map from topology: boundary loops + radial march."""
    adj = defaultdict(list)
    for a, b in bars: adj[a].append(b); adj[b].append(a)
    seen, loops = set(), []
    for s in adj:
        if s in seen: continue
        loop, cur, prev = [s], s, None
        seen.add(s)
        while True:
            nxt = [v for v in adj[cur] if v != prev and v not in seen]
            if not nxt: break
            prev, cur = cur, nxt[0]; seen.add(cur); loop.append(cur)
        loops.append(loop)
    if len(loops) != 2: raise SystemExit('expected 2 boundary loops, got %d' % len(loops))
    # the wall is the loop of smaller radial extent
    ext = [np.hypot(xy[l][:, 0], xy[l][:, 1]).max() for l in loops]
    wall = loops[int(np.argmin(ext))]
    NI = len(wall); NJ = len(quads) // NI + 1
    if NI * (NJ - 1) != len(quads): raise SystemExit('not a structured O-grid')

    q_of = defaultdict(list)
    for qi, q in enumerate(quads):
        for v in q: q_of[v].append(qi)

    def outer(a, b, banned):
        """The quad on edge (a,b) not containing the previous ring, and its far pair."""
        cand = [quads[x] for x in set(q_of[a]) & set(q_of[b])]
        cand = [q for q in cand if not banned & set(q)]
        if len(cand) != 1: raise SystemExit('ambiguous march at %d-%d' % (a, b))
        q = list(cand[0]); ia = q.index(a)
        if q[(ia + 1) % 4] == b: return q[(ia + 3) % 4], q[(ia + 2) % 4]
        return q[(ia + 1) % 4], q[(ia + 2) % 4]

    g = np.zeros((NI, NJ), dtype=int); g[:, 0] = wall
    for jj in range(NJ - 1):
        banned = set() if jj == 0 else set(g[:, jj - 1])
        for i in range(NI):
            a, b = g[i, jj], g[(i + 1) % NI, jj]
            oa, ob = outer(a, b, banned)
            g[i, jj + 1] = oa; g[(i + 1) % NI, jj + 1] = ob
    if len(set(g.ravel())) != NI * NJ: raise SystemExit('march is not a bijection')
    return xy[g][:, :, 0], xy[g][:, :, 1]


def _p3d_dims(path, hdr, nval):
    """(ni, nj) from the integer header, whichever PLOT3D flavour wrote it:
       [ni nj], [ni nj nk], [nblk ni nj] or [nblk ni nj nk].  The value count
       (x,y or x,y,z) decides, so the flavour is never guessed."""
    cand = set()
    for nb, d in ((1, hdr), (hdr[0] if hdr else 0, hdr[1:])):
        if nb != 1 or len(d) not in (2, 3): continue
        ni, nj = d[0], d[1]; nk = d[2] if len(d) == 3 else 1
        if nk == 1 and ni > 1 and nj > 1 and nval in (2 * ni * nj, 3 * ni * nj):
            cand.add((ni, nj))
    if hdr and hdr[0] > 1 and len(hdr) >= 3 and not cand:
        raise SystemExit('%s: header %s -- multi-block or 3D file, this tool wants ONE 2D block' % (path, hdr))
    if len(cand) != 1:
        raise SystemExit('%s: header %s does not match %d coordinate values' % (path, hdr, nval))
    return cand.pop()


def _fortran_records(raw):
    """Sequential unformatted records (4-byte little-endian markers), or None."""
    recs, p = [], 0
    while p < len(raw):
        n = int.from_bytes(raw[p:p + 4], 'little')
        if p + 8 + n > len(raw) or int.from_bytes(raw[p + 4 + n:p + 8 + n], 'little') != n:
            return None
        recs.append(raw[p + 4:p + 4 + n]); p += 8 + n
    return recs


def read_plot3d(path):
    """Single-block 2D PLOT3D -> node arrays (ni, nj), i fastest on file.
       ASCII or Fortran unformatted (double or single), with or without the
       block count and the nk = 1 entry."""
    raw = open(path, 'rb').read()
    if b'\x00' in raw[:256]:
        recs = _fortran_records(raw)
        if not recs or len(recs) < 2:
            raise SystemExit('%s: binary but not Fortran sequential records; re-export as ASCII' % path)
        hdr = [int(k) for r in recs[:-1] for k in np.frombuffer(r, '<i4')]
        for dt in ('<f8', '<f4'):
            v = np.frombuffer(recs[-1], dt).astype(float)
            try:
                ni, nj = _p3d_dims(path, hdr, v.size); break
            except SystemExit:
                if dt == '<f4': raise
    else:
        t = raw.decode().split()
        nint = next(k for k, s in enumerate(t) if not re.fullmatch(r'[+-]?\d+', s))
        hdr, v = [int(s) for s in t[:nint]], np.array(t[nint:], dtype=float)
        ni, nj = _p3d_dims(path, hdr, v.size)
    nb = ni * nj
    X = v[:nb].reshape(nj, ni).T          # -> (ni, nj)
    Y = v[nb:2 * nb].reshape(nj, ni).T
    return X, Y


def drop_wrap(X, Y):
    """An O-grid .x repeats the seam column; drop it and return the periodic set."""
    if np.allclose(X[0], X[-1]) and np.allclose(Y[0], Y[-1]):
        return X[:-1], Y[:-1], True
    return X, Y, False


def rotate(X, Y, aoa_deg):
    """Rigid rotation by -aoa about the leading edge, so a +x stream sees +aoa."""
    t = -np.radians(aoa_deg); c, s = np.cos(t), np.sin(t)
    return c * X - s * Y, s * X + c * Y


def to_clockwise(X, Y):
    """Reverse the wrap only if the wall loop is counter-clockwise. Returns (X,Y,flipped)."""
    xw, yw = X[:, 0], Y[:, 0]
    area = 0.5 * np.sum(xw * np.roll(yw, -1) - np.roll(xw, -1) * yw)
    if area < 0:
        return X, Y, False
    idx = (-np.arange(X.shape[0])) % X.shape[0]
    return X[idx], Y[idx], True


def build(Xn, Yn):
    """Node lines -> ghost-inclusive URANOS cell-centre arrays, indexed [i,j] from 1-NG."""
    nx, ny = Xn.shape[0], Xn.shape[1] - 1
    xw, yw = Xn[:, 0], Yn[:, 0]
    area = 0.5 * np.sum(xw * np.roll(yw, -1) - np.roll(xw, -1) * yw)
    if area >= 0:
        raise SystemExit('the wall loop is still counter-clockwise: detJ2D would be negative')

    Xc = 0.5 * (Xn[:, :-1] + Xn[:, 1:]); Yc = 0.5 * (Yn[:, :-1] + Yn[:, 1:])
    Xw, Yw = xw, yw
    XX = np.zeros((nx + 2 * NG, ny + 2 * NG)); YY = np.zeros_like(XX)
    XX[NG:NG + nx, NG:NG + ny] = Xc; YY[NG:NG + nx, NG:NG + ny] = Yc

    seg = np.hypot(np.diff(Xc, axis=1), np.diff(Yc, axis=1))
    s = np.zeros((nx, ny)); s[:, 0] = np.hypot(Xc[:, 0] - Xw, Yc[:, 0] - Yw)
    for k in range(1, ny): s[:, k] = s[:, k - 1] + seg[:, k - 1]

    # launch angle of the ghost ray, smoothed periodically: the raw angle steps
    # by tens of degrees at a blunt-base corner or a sharp trailing edge and
    # folds the deep ghost rows
    ang = np.arctan2(Yc[:, 0] - Yw, Xc[:, 0] - Xw)
    for _ in range(NSM):
        p, n = np.roll(ang, 1), np.roll(ang, -1)
        p += 2 * np.pi * np.rint((ang - p) / (2 * np.pi))
        n += 2 * np.pi * np.rint((ang - n) / (2 * np.pi))
        ang = 0.25 * p + 0.5 * ang + 0.25 * n
    ux, uy = np.cos(ang), np.sin(ang)

    for k in range(1, NG + 1):                                # wall ghosts j = 0..-6
        if k <= NREF:
            XX[NG:NG + nx, NG - k] = 2 * Xw - Xc[:, k - 1]
            YY[NG:NG + nx, NG - k] = 2 * Yw - Yc[:, k - 1]
        else:
            XX[NG:NG + nx, NG - k] = Xw - s[:, k - 1] * ux
            YY[NG:NG + nx, NG - k] = Yw - s[:, k - 1] * uy

    dX, dY = Xc[:, -1] - Xc[:, -2], Yc[:, -1] - Yc[:, -2]     # far-field ghosts
    r = np.clip(np.hypot(dX, dY) / np.hypot(Xc[:, -2] - Xc[:, -3], Yc[:, -2] - Yc[:, -3]), 1.0, 1.3)
    px, py, sx, sy = Xc[:, -1].copy(), Yc[:, -1].copy(), dX.copy(), dY.copy()
    for k in range(1, NG + 1):
        sx *= r; sy *= r; px, py = px + sx, py + sy
        XX[NG:NG + nx, NG + ny - 1 + k] = px; YY[NG:NG + nx, NG + ny - 1 + k] = py

    for k in range(1, NG + 1):                                # csi ghosts: exact wrap
        XX[NG - k, :], YY[NG - k, :] = XX[NG + nx - k, :], YY[NG + nx - k, :]
        XX[NG + nx - 1 + k, :], YY[NG + nx - 1 + k, :] = XX[NG + k - 1, :], YY[NG + k - 1, :]
    return XX, YY, nx, ny


def min_detJ(XX, YY, nx, ny):
    """URANOS' own fd6 detJ2D over 1-GN..n+GN (mesh_module.f90:594-604)."""
    c = (3. / 4., -3. / 20., 1. / 60.)
    i = np.arange(NG - GN, NG + nx + GN)[:, None]
    j = np.arange(NG - GN, NG + ny + GN)[None, :]
    xx = xe = yx = ye = 0.0
    for l in (1, 2, 3):
        xx = xx + c[l - 1] * (XX[i + l, j] - XX[i - l, j]); yx = yx + c[l - 1] * (YY[i + l, j] - YY[i - l, j])
        xe = xe + c[l - 1] * (XX[i, j + l] - XX[i, j - l]); ye = ye + c[l - 1] * (YY[i, j + l] - YY[i, j - l])
    return float((xx * ye - xe * yx).min())


def convert(src, out, aoa, verbose=True):
    X, Y = read_mesh(src)
    ni, nj = X.shape
    X, Y, dup = drop_wrap(X, Y)
    X, Y = rotate(X, Y, aoa)
    X, Y, flipped = to_clockwise(X, Y)
    XX, YY, nx, ny = build(X, Y)
    mn = min_detJ(XX, YY, nx, ny)
    if mn <= 1e-12:
        raise SystemExit('min detJ2D = %.3e <= 1e-12: URANOS would abort at grid setup' % mn)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as f:
        f.write('%d %d %d\n' % (nx, ny, 1))
        for A in (XX, YY):
            f.write(''.join('%25.16E\n' % t for t in A.T.ravel()))
    if verbose:
        print('  source    %s  ni=%d nj=%d%s' % (os.path.basename(src), ni, nj,
              '  (duplicated wrap line removed)' if dup else '  (unique nodes, no wrap line to remove)'))
        print('  verse     %s' % ('counter-clockwise -> REVERSED' if flipped else 'already clockwise, unchanged'))
        print('  alpha     %.3f deg baked into the mesh' % aoa)
        print('  grid      nx=%d  ny=%d   (cells %d)' % (nx, ny, nx * ny))
        print('  min detJ2D over 1-GN..n+GN = %.4e   > 1e-12  ok' % mn)
        print('  written   %s' % out)
    return nx, ny, mn


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        raise SystemExit('usage: python3 mesh_to_uranos.py <in.x|in.msh> <out.txt> [aoa_deg]')
    src, out = sys.argv[1], sys.argv[2]
    aoa = float(sys.argv[3]) if len(sys.argv) > 3 else AOA
    convert(src, out, aoa)
