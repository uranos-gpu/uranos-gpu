#!/usr/bin/env python3
# ============================================================================
#  make_ramp_mesh.py  --  ramp mesh driver
#
#  Compiles and runs the Fortran generator ramp_mesh_gen.f90 (Adams two-boundary
#  blend, rounded corner, constant-thickness ribbon) and writes into
#  mesh/<retau-dir>/ :
#      MeshXY_ramp.txt        <- geometry read by URANOS (mesh_file)
#      ramp_grid.namelist.txt <- grid block of the deck (nx, ny, nz, extents, ReTau)
#
#  NOTE ON THE TWO ReTau's (they legitimately differ):
#    * MESH ReTau (in the generator &input_list) -> sets the wall-unit CLUSTERING
#      (y+, Dx+, Dz+).  Cluster the mesh for the resolution you want.
#    * INFLOW ReTau (--inflow-retau, written into ramp_grid.namelist.txt's ReTau)
#      -> the deck 'Reynolds' field = digital-filter INFLOW Re_tau (GetRetau ->
#      free-stream Re).  This is the CALIBRATION knob (e.g. mesh@340, inflow@440).
#    If --inflow-retau is omitted it defaults to the generator's ReTau.
#
#  USAGE (from anywhere; paths are resolved relative to this script):
#    python3 make_ramp_mesh.py --input <gen_namelist.dat> --retau-dir retau440 \
#            [--inflow-retau 440] [--rebuild] [--view]
#
#  The DNS mesh of the paper: python3 make_ramp_mesh.py --retau-dir dns --inflow-retau 380
# ============================================================================
import os, sys, argparse, subprocess, shutil, re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))   # repo root (DATA/<case>/mesh)
GEN_F90 = os.path.join(HERE, "ramp_mesh_gen.f90")
GEN_EXE = os.path.join(HERE, "ramp_mesh_gen.exe")
VIEWER  = os.path.join(HERE, "view_mesh.py")


def parse_namelist(path):
    """Flat &input_list parser: key -> raw string value (comments/commas stripped)."""
    vals = {}
    for line in open(path):
        line = line.split("!")[0].strip().rstrip(",").strip()
        if "=" not in line or line.startswith("&") or line.startswith("/"):
            continue
        k, v = line.split("=", 1)
        vals[k.strip()] = v.strip().strip(",").strip()
    return vals


def build_exe(rebuild=False):
    if rebuild or not os.path.exists(GEN_EXE):
        print(f"[make_ramp_mesh] compiling {os.path.basename(GEN_F90)} -> {os.path.basename(GEN_EXE)}")
        subprocess.run(["gfortran", "-O2", GEN_F90, "-o", GEN_EXE], check=True, cwd=HERE)
    return GEN_EXE


def num(v):
    """'-29.0' -> -29.0 ; '896' -> 896 ; keep as-is on failure."""
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def main():
    ap = argparse.ArgumentParser(description="Generate a ramp mesh and its deck grid block in one folder.")
    ap.add_argument("--retau-dir", required=True,
                    help="output subfolder under mesh/ (e.g. retau440, dns, iles). Created if absent.")
    ap.add_argument("--input", default=None,
                    help="generator &input_list .dat (domain knobs). "
                         "Default: <retau-dir>/mesh_gen_input_ramp_mesh_gen.dat")
    ap.add_argument("--inflow-retau", type=float, default=None,
                    help="INFLOW Re_tau written to the deck (default = generator's ReTau)")
    ap.add_argument("--rebuild", action="store_true", help="force recompile of ramp_mesh_gen.exe")
    ap.add_argument("--view", action="store_true", help="also render a sanity PNG via view_mesh.py")
    a = ap.parse_args()

    if a.input is None:
        a.input = os.path.join(HERE, a.retau_dir, "mesh_gen_input_ramp_mesh_gen.dat")
    if not os.path.exists(a.input):
        sys.exit(f"ERROR: generator input not found: {a.input}")

    exe = build_exe(a.rebuild)
    gen = parse_namelist(a.input)
    for req in ("nx", "ny", "nz", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax", "Mach", "ReTau"):
        if req not in gen:
            sys.exit(f"ERROR: '{req}' missing from {a.input}")

    target = os.path.join(HERE, a.retau_dir)
    os.makedirs(target, exist_ok=True)

    # 1) working generator input in the target dir, out_file forced to MeshXY_ramp.txt
    work_in = os.path.join(target, "mesh_gen_input_ramp_mesh_gen.dat")
    with open(a.input) as f:
        txt = f.read()
    if re.search(r"out_file", txt):
        txt = re.sub(r"out_file\s*=\s*'[^']*'", "out_file = 'MeshXY_ramp.txt'", txt)
    else:
        txt = txt.replace("&end", "    out_file = 'MeshXY_ramp.txt'\n&end")
    open(work_in, "w").write(txt)

    # 2) run the Fortran generator IN the target dir -> writes MeshXY_ramp.txt there
    print(f"[make_ramp_mesh] running generator -> {a.retau_dir}/MeshXY_ramp.txt")
    r = subprocess.run([exe, "mesh_gen_input_ramp_mesh_gen.dat"], cwd=target,
                       capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0 or not os.path.exists(os.path.join(target, "MeshXY_ramp.txt")):
        sys.exit("ERROR: generator failed:\n" + r.stderr)

    # 3) INFLOW Re_tau for the deck (calibration knob), may differ from mesh ReTau
    inflow_retau = a.inflow_retau if a.inflow_retau is not None else num(gen["ReTau"])
    mesh_rel = os.path.relpath(os.path.join(target, "MeshXY_ramp.txt"), REPO).replace(os.sep, "/")

    # 4) ramp_grid.namelist.txt -- grid block of the deck
    nml = os.path.join(target, "ramp_grid.namelist.txt")
    with open(nml, "w") as f:
        f.write("! ---- grid block for the URANOS .dat (curvilinear ramp) ----\n")
        f.write("! generated by make_ramp_mesh.py wrapping ramp_mesh_gen.f90 (Adams two-boundary blend)\n")
        f.write(f"! mesh clustered for Re_tau={num(gen['ReTau'])}; INFLOW Re_tau (deck)={inflow_retau}\n")
        f.write(f"nx = {gen['nx']}\n")
        f.write(f"ny = {gen['ny']}\n")
        f.write(f"nz = {gen['nz']}\n")
        f.write(f"xmin = {float(gen['xmin']):.6f}\n")
        f.write(f"xmax = {float(gen['xmax']):.6f}\n")
        f.write(f"ymin = {float(gen['ymin']):.6f}\n")
        f.write(f"ymax = {float(gen['ymax']):.6f}\n")
        f.write(f"zmin = {float(gen['zmin']):.6f}\n")
        f.write(f"zmax = {float(gen['zmax']):.6f}\n")
        f.write("gridpoint_x = 'uniform'\n")
        f.write("gridpoint_y = 'uniform'\n")
        f.write("gridpoint_z = 'uniform'\n")
        f.write(f"mesh_file   = '{mesh_rel}'\n")
        f.write("! ReTau = INFLOW value (GetRetau -> free-stream Re); the calibration knob\n")
        f.write(f"ReTau = {float(inflow_retau):.1f}\n")
        f.write(f"Mach  = {float(gen['Mach']):.3f}\n")
    print(f"[make_ramp_mesh] wrote {a.retau_dir}/ramp_grid.namelist.txt (nx={gen['nx']} ny={gen['ny']} "
          f"nz={gen['nz']}, inflow ReTau={inflow_retau}, mesh_file={mesh_rel})")

    # 5) optional visual sanity check
    if a.view:
        png = os.path.join(target, "_mesh_view.png")
        subprocess.run(["python3", VIEWER, "MeshXY_ramp.txt", "--H", str(gen["ymax"]),
                        "--retau", str(inflow_retau), "--save", png], cwd=target)
        print(f"[make_ramp_mesh] view -> {png}")

    print(f"\n[make_ramp_mesh] DONE. The decks read mesh/{a.retau_dir}/{{MeshXY_ramp.txt, ramp_grid.namelist.txt}}")


if __name__ == "__main__":
    main()
