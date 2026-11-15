#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# Generate a 2D "fake-Cartesian" body-fitted mesh for URANOS: a tensor-product
# grid whose nodes coincide exactly with the 1D uniform mapping
# x(i) = xmin + Lx*(i-0.5)/nx ,  y(j) = ymin + Ly*(j-0.5)/ny .
#
# Feeding this mesh through the curvilinear path must reproduce the native 1D
# Cartesian result bit-for-bit, so it is the regression reference for the
# curvilinear port.
#
# File layout expected by mesh_module::read_curvilinear_mesh:
#   line 1 : nx ny nz
#   then   : all x values, i fastest then j,  over [1-ng, n+ng] per side
#   then   : all y values, same loop order
# with ng = GN + mesh_gh = 4 + 3 = 7 ghost nodes per side.
# ----------------------------------------------------------------------------
import sys

def main():
    # usage: gen_fake_cartesian_mesh.py [nx ny nz xmin xmax ymin ymax out]
    # no args -> shock_tube_x.dat geometry
    a = sys.argv[1:]
    if len(a) == 8:
        nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
        xmin, xmax, ymin, ymax = (float(v) for v in a[3:7])
        out = a[7]
    else:
        nx, ny, nz = 101, 5, 5
        xmin, xmax = -5.0, 5.0
        ymin, ymax = -1.0, 1.0
        out = "grids/fake_cartesian_shock_tube_x.txt"
    ng = 7                      # GN(4) + mesh_gh(3)

    Lx = xmax - xmin
    Ly = ymax - ymin

    # match uniformGrid_f0 operation order exactly (xi=(i-0.5)/nx ; x=xmin+Lx*xi)
    # so the mesh nodes coincide bit-for-bit with the 1D tensor extrusion
    xv = [xmin + Lx*((i-0.5)/nx) for i in range(1-ng, nx+ng+1)]
    yv = [ymin + Ly*((j-0.5)/ny) for j in range(1-ng, ny+ng+1)]

    with open(out, "w") as f:
        f.write(f"{nx} {ny} {nz}\n")
        # x block: i fastest, j slowest
        for j in range(1-ng, ny+ng+1):
            for i in range(1-ng, nx+ng+1):
                f.write(f"{xv[i-(1-ng)]!r}\n")
        # y block
        for j in range(1-ng, ny+ng+1):
            for i in range(1-ng, nx+ng+1):
                f.write(f"{yv[j-(1-ng)]!r}\n")

    nrow = nx + 2*ng
    ncol = ny + 2*ng
    print(f"wrote {out}: {nrow} x {ncol} nodes per block (ng={ng}), "
          f"x in [{xv[0]:.6g},{xv[-1]:.6g}], y in [{yv[0]:.6g},{yv[-1]:.6g}]")

if __name__ == "__main__":
    main()
