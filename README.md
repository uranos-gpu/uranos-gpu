# URANOS — A GPU-Accelerated Navier–Stokes Solver for Compressible Wall-Bounded Flows

[![License: BSD](https://img.shields.io/badge/License-BSD-blue.svg)](#license)
[![Fortran](https://img.shields.io/badge/language-Fortran90-informational.svg)]()
[![GPU](https://img.shields.io/badge/GPU-OpenACC-success.svg)]()
[![MPI](https://img.shields.io/badge/Parallel-MPI-informational.svg)]()

URANOS is a massively parallel solver for DNS/LES/WMLES of **compressible** wall-bounded flows, designed for modern pre-exascale systems and workstations. It supports both **CPU** and **GPU** execution with **OpenACC** + MPI, and ships with post-processing utilities to produce VTK and line data for analysis and visualization. :contentReference[oaicite:0]{index=0}

---

## Table of Contents

- [Features](#features)
- [Getting Started](#getting-started)
- [Build](#build)
  - [Dependencies](#dependencies)
  - [Build Modes](#build-modes)
  - [Examples](#examples)
- [Run](#run)
  - [Local (CPU/GPU)](#local-cpugpu)
  - [SLURM Examples](#slurm-examples)
- [Test Cases](#test-cases)
- [Input File (`file.dat`) Cheat-Sheet](#input-file-filedat-cheat-sheet)
- [Outputs & Post-Processing](#outputs--post-processing)
- [Performance & Portability](#performance--portability)
- [Citing URANOS](#citing-uranos)
- [Publications Using URANOS](#publications-using-uranos)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [FAQ / Troubleshooting](#faq--troubleshooting)

---

## Features

- **High-fidelity compressible flow**: DNS, LES, and WMLES with energy-preserving central schemes and hybrid WENO/TENO families.
- **Wall models & SGS**: WALE, Smagorinsky, Sigma; wall-temperature control via `Trat` (adiabatic/cold/hot walls).
- **Shock capturing**: density, density-jump, and Ducros shock sensors.
- **GPU portability** via OpenACC (NVIDIA & AMD through vendor toolchains), and MPI domain decomposition. :contentReference[oaicite:1]{index=1}
- **Post-processing** tool (`PostUranos.exe`) to export **VTK** (3D fields) and **VTK2D** (planes/lines) and GNUPLOT-ready data.

---

## Getting Started

Clone:

```bash
git clone https://github.com/uranos-gpu/uranos-gpu.git
cd uranos-gpu
````

Recommended first run: the **Shock Tube (1D)** or **Channel DNS** tests (see [Test Cases](#test-cases)).

---

## Build

### Dependencies

* **Fortran compiler**: `gfortran`, `nvfortran` (NVIDIA HPC SDK), or Cray Fortran
* **MPI library**: OpenMPI, MPICH, or Cray MPI
* **(GPU optional)**: vendor toolchain (e.g., NVIDIA HPC SDK or Cray environment)

### Build Modes

Use the provided `Makefile`:

```
make -j <nproc> comp="<compiler>" mode="<build_mode>"
```

* `comp`: `gnu` | `nvhpc` | `cray`
* `mode`: `cpu` | `cpu_debug` | `gpu` | `gpu_profiling`

### Examples

**GNU, optimized CPU**

```bash
make -j comp=gnu mode=cpu
```

**GNU, debug CPU (bounds & runtime checks)**

```bash
make -j comp=gnu mode=cpu_debug
```

**NVIDIA GPU (V100/A100/H100), optimized**

```bash
make -j comp=nvhpc mode=gpu
```

**NVIDIA GPU, profiling with NVTX**

```bash
make -j comp=nvhpc mode=gpu_profiling
```

**Cray/AMD GPU (e.g., LUMI), optimized**

```bash
make -j comp=cray mode=gpu
```

---

## Run

URANOS uses an input file `file.dat` and a standard MPI launcher.

### Local (CPU/GPU)

**Basic (no restart)**

```bash
mpirun -np <nprocs> ./Uranos.exe path/to/file.dat
```

**With restart**

```bash
mpirun -np <nprocs> ./Uranos.exe path/to/file.dat path/to/restart.bin
```

**Single-GPU workstation**

```bash
mpirun -np 1 ./Uranos.exe ./tests/flat_plate/input.dat
```

**Multi-GPU workstation (2 ranks → 2 GPUs)**

```bash
mpirun -np 2 ./Uranos.exe ./tests/flat_plate/input.dat
```

> **Tip:** On GPU nodes, set ranks = GPUs per node and load your site’s `nvhpc`/MPI modules.

### SLURM Examples

**Generic CPU (1 node, 32 ranks)**

```bash
#!/bin/bash
#SBATCH -J uranos_cpu
#SBATCH -p cpu_partition
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

module purge
module load gcc/12.1.0 openmpi/4.1.5

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat
```

**Generic GPU (1 node, 4 GPUs → 4 ranks)**

```bash
#!/bin/bash
#SBATCH -J uranos_gpu
#SBATCH -p gpu_partition
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

module purge
module load nvhpc/24.3 openmpi/4.1.5   # adjust to your site

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat
```

**CINECA Leonardo (NVIDIA, 4 GPUs/node; 2 nodes → 8 GPUs)**

```bash
#!/bin/bash
#SBATCH -J uranos_leonardo
#SBATCH -p boost_usr_prod
#SBATCH --time=24:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

module purge
module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8
module load nvhpc/23.1

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat
```

**LUMI (AMD MI250/MI300 via Cray stack)**

```bash
#!/bin/bash
#SBATCH -J uranos_lumi
#SBATCH -A <project>
#SBATCH -p <partition>
#SBATCH -N 2
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=8

module purge
module load PrgEnv-cray
module load cray-mpich

srun -n ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat
```

> **Note:** LUMI/Leonardo specifics may evolve; check your center docs for current module versions.

---

## Test Cases

All examples live in `tests/` and can be post-processed with `PostUranos.exe`.

### 1) Shock Tube (1D)

```bash
mpirun -np 4 ./Uranos.exe tests/shock_tube_x.dat
./PostUranos.exe tests/shock_tube_x.dat
```

Output directory: `DATA/SHOCK_TUBE_1D/`

### 2) Channel DNS (3D, M\_b≈1.5)

```bash
mpirun -np 4 ./Uranos.exe tests/channel_dns.dat
./PostUranos.exe tests/channel_dns.dat
```

Wall-normal stats during runtime:

* `DATA/CHANNEL_DNS/VEL_STATS`
* `DATA/CHANNEL_DNS/BUD_STATS`

### 3) Hypersonic Boundary Layer (WMLES)

```bash
mpirun -np 4 ./Uranos.exe tests/hypersonic_boundary_layer.dat
./PostUranos.exe tests/hypersonic_boundary_layer.dat
```

Stations/lines are saved per downstream location for both line plots and contours.

---

## Input File (`file.dat`) Cheat-Sheet

> **Grid & Domain**

* `xmin,xmax,ymin,ymax,zmin,zmax` — domain extents
* `uniform` — use built-in uniform grid (no external files)
* `gridpoint_x/y/z` — paths to grid files (if not uniform)
* `dims` — problem dimensionality (`2` or `3`)
* `nx,ny,nz` — grid points (must match grid files)

> **Time Integration**

* `CFL` — Courant number (RK stable up to 1; keep ≤ 0.8 for safety)
* `logical_CFL` — `.true.` for adaptive Δt; `.false.` for fixed `Dt`
* `itmax`, `tmax` — iteration/time limits

> **Physics & Models**

* `Reynolds`, `Mach`, `Prandtl`
* `Trat` — wall-to-adiabatic temperature ratio (`=1` adiabatic; `<1` cold; `>1` hot)
* `LES` — enable SGS; `sgs_model = WALE | Smagorinsky | sigma_model`
* `viscous` — enable viscous fluxes

> **Schemes & Sensors**

* `scheme` — `energy_preserving | hybrid_wenoEP | hybrid_tenoEP | hybrid_tenoaEP`
* `sensor` — `density | density_jump | ducros`
* `fd_order` — central order (2,4,6); `weno_order` — WENO/TENO (3,5,7)

> **I/O & Decomposition**

* `cart_dims(1:3)` — MPI cartesian split (X,Y,Z)
* `output_file_name`, `data_dir`
* `itout` — restart cadence; `StOut` — statistics cadence; `StFlg` — enable stats

> **BCs & Extras**

* `bc(1..6)` — BC per face (see `src/bc_module.f90`)
* `sponge(1..6)` — activate sponge zones per face
* `ic` — initial condition (see `src/ic_module.f90`)
* `inflow_profile`, `smooth_inflow`, `turb_inflow`

---

## Outputs & Post-Processing

URANOS writes binary data under:

```
DATA/<data_dir>/BINARY/*.bin
```

Post-process with:

```bash
./PostUranos.exe path/to/file.dat [path/to/restart.bin]
```

Generated products:

* `DATA/<data_dir>/VTK/`   → 3D fields (VTK)
* `DATA/<data_dir>/VTK2D/` → 2D fields/planes (VTK)
* GNUPLOT-ready line data for quick profiling

> If you don’t see VTK outputs after running `PostUranos.exe`, verify the `data_dir` in your `file.dat` matches the simulation output path and that `BINARY/*.bin` exists for the requested times. (Common issue: mismatched `output_file_name`/`data_dir` between run and post.)

---

## Performance & Portability

URANOS emphasizes portability across NVIDIA and AMD GPU architectures via OpenACC (tested on EuroHPC Leonardo and LUMI in recent releases). See the CPC articles for detailed performance analyses and design choices.

---

## Citing URANOS

If you use URANOS in academic work, please cite:

* **URANOS-2.0** — *Computer Physics Communications*, 2024. [https://doi.org/10.1016/j.cpc.2024.109285](https://doi.org/10.1016/j.cpc.2024.109285)
* **URANOS (original)** — *Computer Physics Communications*, 2023. [https://doi.org/10.1016/j.cpc.2023.108717](https://doi.org/10.1016/j.cpc.2023.108717)

BibTeX snippets are available at the DOI pages.

---

## Publications Using URANOS

A selection (recent first):

* De Vanna & Benini (2025). *WMLES of a Transonic Gas Turbine Vane – Part I.* **J. Turbomach.** [https://doi.org/10.1115/1.4069131](https://doi.org/10.1115/1.4069131)
* De Vanna (2025). *Entropy losses in transonic shock–boundary-layer interaction.* **Phys. Fluids.** [https://doi.org/10.1063/5.0278759](https://doi.org/10.1063/5.0278759)
* De Vanna & Benini (2025). *Impact of Wall Cooling on Transonic Gas Turbine Stators Aerothermodynamics.* **Appl. Therm. Eng.** [https://doi.org/10.1016/j.applthermaleng.2025.126396](https://doi.org/10.1016/j.applthermaleng.2025.126396)
* De Vanna, Baldan, Picano, Benini (2023). *WMLES + IBM for compressible flows.* **Computers & Fluids.** [https://doi.org/10.1016/j.compfluid.2023.106058](https://doi.org/10.1016/j.compfluid.2023.106058)
* De Vanna et al. (2022). *Effect of convective schemes in WR/WMLES.* **Computers & Fluids.** [https://doi.org/10.1016/j.compfluid.2022.105710](https://doi.org/10.1016/j.compfluid.2022.105710)
* De Vanna, Bernardini, Picano, Benini (2022). *WMLES of SWBLI.* **Int. J. Heat Fluid Flow.** [https://doi.org/10.1016/j.ijheatfluidflow.2022.109071](https://doi.org/10.1016/j.ijheatfluidflow.2022.109071)
* De Vanna et al. (2021). *Unified WR/WMLES method.* **Phys. Rev. Fluids.** [https://doi.org/10.1103/PhysRevFluids.6.034614](https://doi.org/10.1103/PhysRevFluids.6.034614)
* De Vanna, Picano, Benini, Quinn (2021). *LES of hypersonic intake (M5).* **AIAA J.** [https://doi.org/10.2514/1.J060160](https://doi.org/10.2514/1.J060160)
* De Vanna, Benato, Picano, Benini (2021). *High-order viscous terms (variable viscosity).* **Acta Mechanica.** [https://doi.org/10.1007/s00707-021-02937-2](https://doi.org/10.1007/s00707-021-02937-2)

---

## Contributing

We welcome issues and PRs! Please open an Issue or Discussion to propose changes before submitting a PR. See `CONTRIBUTING.md` (if absent, include testing notes in your PR).

---

## License

This project is released under a BSD-style license. See [`LICENSE`](./LICENSE) for details.

---

## Acknowledgments

Developed at the **Department of Industrial Engineering, University of Padova**, with runs on EuroHPC systems (Leonardo, LUMI). We thank collaborators and HPC centers for support.

---

## FAQ / Troubleshooting

**Q1. `PostUranos.exe` runs but VTK folders are empty.**
A: Ensure `DATA/<data_dir>/BINARY/*.bin` exists and matches the `file.dat` you pass to `PostUranos.exe`. If multiple runs share a `data_dir`, point to the intended `restart.bin` explicitly when post-processing.

**Q2. How many MPI ranks should I use on GPUs?**
A: Typically 1 rank per GPU on a node. For multi-node runs, set `--ntasks-per-node` to the number of GPUs per node and use the site-provided CUDA/ROCm/MPI modules.

**Q3. Which scheme should I use for shocked flows?**
A: Prefer hybrid WENO/TENO variants with an appropriate shock sensor; avoid purely energy-preserving central schemes for strong shocks.

**Q4. Typical CFL values?**
A: The RK scheme is stable up to CFL ≈ 1; for safety, use ≤ 0.8. Start conservative and increase once stable.

**Q5. Where do statistics go in channel/BL cases?**
A: See the test-case sections—velocity and budget stats are written during runtime in dedicated subfolders.

