#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# Generate a smooth, PERIODIC, genuinely-curved 2D body-fitted mesh for URANOS.
# A uniform grid X(i),Y(j) is warped so that the off-diagonal metrics
# x_eta, y_csi are non-zero (true curvilinear coupling):
#     xx(i,j) = X(i) + amp*sin(2*pi*(Y(j)-ymin)/Ly)
#     yy(i,j) = Y(j) + amp*sin(2*pi*(X(i)-xmin)/Lx)
# The warp vanishes at the periodic images (sin over a full period), so the
# mesh is periodic-compatible.  Used for the freestream/GCL test: a uniform
# flow must stay uniform on this mesh if the metric identities hold discretely.
#
# [!] That warp is ADDITIVELY SEPARABLE: x_csi and y_eta are constant, x_eta
#     depends on j alone and y_csi on i alone.  Weighted by the cell volume the
#     Jacobian cancels and every cofactor left in front of a difference is
#     CONSTANT along the index it differentiates, so discrete integration by
#     parts is EXACT there -- for every scheme.  Any CONSERVATION test on this
#     mesh (and on a Cartesian one) therefore cannot fail, and says nothing.
#     For those, use mode 'nonsep':
#         xx(i,j) = X(i) + amp*sin(kx*X(i))*sin(ky*Y(j))
#         yy(i,j) = Y(j) + amp*sin(kx*X(i))*sin(2*ky*Y(j))
#     whose cofactors depend on BOTH indices.  The separable warp stays the
#     right choice for ORDER and decomposition-invariance tests.
#
# usage: gen_curved_mesh.py nx ny nz xmin xmax ymin ymax amp out [sep|nonsep]
# ----------------------------------------------------------------------------
import sys, math

def main():
    a = sys.argv[1:]
    nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    xmin, xmax, ymin, ymax = (float(v) for v in a[3:7])
    amp = float(a[7])
    out = a[8]
    mode = a[9] if len(a) > 9 else 'sep'
    if mode not in ('sep', 'nonsep'):
        sys.exit("mode must be 'sep' or 'nonsep'")
    ng = 7                      # GN(4) + mesh_gh(3)

    Lx, Ly = xmax - xmin, ymax - ymin
    X = [xmin + Lx*((i-0.5)/nx) for i in range(1-ng, nx+ng+1)]
    Y = [ymin + Ly*((j-0.5)/ny) for j in range(1-ng, ny+ng+1)]

    if mode == 'sep':
        # arithmetic kept EXACTLY as before (same association) so the separable
        # meshes already on disk are reproduced byte for byte
        def wx(xv, yv): return amp*math.sin(2.0*math.pi*(yv-ymin)/Ly)   # x-shift, function of Y alone
        def wy(xv, yv): return amp*math.sin(2.0*math.pi*(xv-xmin)/Lx)   # y-shift, function of X alone
    else:
        def wx(xv, yv): return amp*math.sin(2.0*math.pi*(xv-xmin)/Lx)*math.sin(2.0*math.pi*(yv-ymin)/Ly)
        def wy(xv, yv): return amp*math.sin(2.0*math.pi*(xv-xmin)/Lx)*math.sin(4.0*math.pi*(yv-ymin)/Ly)

    with open(out, "w") as f:
        f.write(f"{nx} {ny} {nz}\n")
        for j in range(len(Y)):                 # x block: i fastest, j slowest
            for i in range(len(X)):
                f.write(f"{X[i] + wx(X[i], Y[j])!r}\n")
        for j in range(len(Y)):                 # y block
            for i in range(len(X)):
                f.write(f"{Y[j] + wy(X[i], Y[j])!r}\n")

    print(f"wrote {out}: curved {nx}x{ny} (ng={ng}), amp={amp}, mode={mode}, "
          f"max|x-shift|={amp:.3g}, off-diagonal metrics active")

if __name__ == "__main__":
    main()
