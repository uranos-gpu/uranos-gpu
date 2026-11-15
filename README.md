![logo](images/logo_2026_bis.png)



# URANOS — Exascaling the Future of Aerospace

[![License: BSD](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](#license)
[![Fortran](https://img.shields.io/badge/language-Fortran90-informational.svg)]()
[![GPU](https://img.shields.io/badge/GPU-OpenACC-success.svg)]()
[![MPI](https://img.shields.io/badge/Parallel-MPI-informational.svg)]()


**URANOS** is a **GPU-accelerated Navier–Stokes solver** for **compressible flows** with direct relevance to **aerospace propulsion systems** and **aerothermodynamics**.

**URANOS** solves the compressible Navier–Stokes equations in generalized curvilinear coordinates on body-fitted, extruded structured grids. It runs on CPUs and on NVIDIA GPUs with **OpenACC** + **MPI**, and scales from a workstation to hundreds of GPU nodes.

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
- [Outputs](#outputs)
- [Performance](#performance)
- [Citing URANOS](#citing-uranos)
- [Publications Using URANOS](#publications-using-uranos)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [FAQ / Troubleshooting](#faq--troubleshooting)

---

## Features

- **Curvilinear framework**: generalized curvilinear coordinates on extruded body-fitted grids, either built in (uniform Cartesian) or read from a `mesh_file`; metric terms satisfy the discrete identities to round-off.
- **Structure-preserving convective fluxes**: kinetic-energy and entropy preserving central fluxes `JP`, `KEEP`, `AEC` and `PEP`, of order 2, 4 or 6.
- **Shock capturing**: hybrid schemes that blend the central flux with WENO-Z (3/5/7), TENO (5/7) or adaptive TENO-A (5) reconstructions in characteristic variables, switched by a density, density-jump or Ducros sensor.
- **Thermodynamics**: calorically perfect gas, or thermally perfect gas with NASA polynomials for air, N2 and O2; power-law, Sutherland or constant viscosity.
- **Viscous operators**: Laplacian, staggered or conservative formulation (the conservative one makes total energy telescope exactly).
- **LES**: Smagorinsky, WALE, sigma and mixed-time-scale subgrid models, optional Yoshizawa isotropic term; wall-resolved LES and implicit LES.
- **Time integration**: three-stage SSP Runge–Kutta with adaptive time step from CFL and Fourier limits, or a fixed `Dt`.
- **Boundary conditions**: periodic, zero-gradient, extrapolation, adiabatic and isothermal no-slip walls, non-reflecting characteristic inflow and outflow, digital-filter turbulent inflow.
- **Parallelism**: MPI 3D domain decomposition with halo-minimizing block choice and node-compact rank placement, CUDA-aware halo exchange, results bitwise independent of the decomposition.

---

## Getting Started

Clone:

```bash
git clone https://github.com/uranos-gpu/uranos-gpu.git
cd uranos-gpu
```

Recommended first run: the **Sod shock tube** or the **isentropic vortex** (see [Test Cases](#test-cases)); both run on a laptop CPU in seconds.

---

## Build

### Dependencies

* **Fortran compiler**: `gfortran` or `nvfortran` (NVIDIA HPC SDK)
* **MPI library**: OpenMPI or MPICH (CUDA-aware for multi-GPU runs)
* **(GPU optional)**: NVIDIA HPC SDK (OpenACC, curand)

### Build Modes

Use the provided `Makefile`:

```
make comp="<compiler>" mode="<build_mode>"
```

* `comp=gnu`: `mode` = `cpu` | `cpu_debug` | `cpu_profiling`
* `comp=nvhpc`: `mode` = `gpu` | `gpu_debug` | `gpu_profiling`
* `gpu_cc=<ccXX>` (optional, `nvhpc` only): target compute capability; the default `ccnative` detects the GPU of the build node, so set it (e.g. `gpu_cc=cc80`) when building on a login node whose GPU differs from the compute nodes.

The Makefile detects the number of cores and compiles in parallel.

### Examples

**GNU, optimized CPU**

```bash
make comp=gnu mode=cpu
```

**GNU, debug CPU (bounds & runtime checks)**

```bash
make comp=gnu mode=cpu_debug
```

**NVIDIA GPU (V100/A100/H100), optimized**

```bash
make comp=nvhpc mode=gpu
```

**NVIDIA GPU, profiling with NVTX**

```bash
make comp=nvhpc mode=gpu_profiling
```

`make clean` removes the build artifacts (`.objs/`, `.mods/`).

---

## Run

URANOS uses an input file `file.dat` (Fortran namelist `&input_list`) and a standard MPI launcher. Run from the repository root: paths inside a deck (`mesh_file`, `data_dir`) are relative to it.

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
mpirun -np 1 ./Uranos.exe DATA/03_sod_shock_tube/teno5.dat
```

**Multi-GPU workstation (2 ranks → 2 GPUs)**

```bash
mpirun -np 2 ./Uranos.exe path/to/file.dat
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
module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8   # adjust to the current Leonardo stack
module load nvhpc/23.1

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat
```

> **Note:** module names and versions change over time; check your center documentation.

---

## Test Cases

The cases of the URANOS-3.0 paper live in `DATA/`, one self-contained folder each with the input decks, the meshes, the scripts that generate them and the scripts that post-process the runs into the paper figures. See [`DATA/README.md`](./DATA/README.md).

| folder | content |
|---|---|
| `01_isentropic_vortex` | isentropic vortex, four KEP fluxes, Cartesian and wavy grid |
| `02_order_of_accuracy` | order of accuracy of the convective and viscous operators |
| `03_sod_shock_tube` | Sod shock tube, six shock-capturing reconstructions, CP and TP gas |
| `04_kolmogorov_flow` | steady Kolmogorov flow, three viscous schemes |
| `05_taylor_green_vortex` | Taylor–Green vortex at Re = 1600 |
| `06_time_integration` | temporal order and stability of the Runge–Kutta scheme |
| `07_supersonic_compression_ramp` | 24° compression ramp at Ma = 2.91, body-fitted DNS |
| `08_multi_gpu_scaling` | strong and weak scaling, triperiodic flow at rest |

Example:

```bash
make comp=gnu mode=cpu
python3 DATA/03_sod_shock_tube/run_sod.py     # six Sod runs, one per reconstruction
python3 DATA/03_sod_shock_tube/plot_sod.py    # density, velocity and pressure profiles
```

Each folder has its own `README.md` with the exact commands.

The Python scripts need `numpy`, `matplotlib` and, for some figures, `scipy` and a LaTeX installation.

---

## Input File (`file.dat`) Cheat-Sheet

> **Grid & Domain**

* `xmin,xmax,ymin,ymax,zmin,zmax` — domain extents
* `nx,ny,nz` — grid points; `dims` — `2` or `3`
* `gridpoint_x/y/z = 'uniform'` — built-in uniform grid
* `mesh_file` — curvilinear in-plane (x-y) grid, extruded uniformly in z; `utils/mesh_to_uranos.py` converts a structured PLOT3D or gmsh O-grid to this format

> **Time Integration**

* `logical_CFL` — `.true.` for adaptive Δt from `CFL` (and `Fourier` for viscous runs); `.false.` for fixed `Dt`
* `n_step = 3` — Runge–Kutta stages
* `itmax`, `tmax` — iteration/time limits

> **Physics & Models**

* `Reynolds`, `Mach`, `Prandtl`, `Tref` (reference temperature in K)
* `thermally_perfect_gas` — `0` calorically perfect, `1` thermally perfect; `species_name = 'air' | 'N2' | 'O2'`
* `vis_flag` — `0` power law, `1` Sutherland, `2` constant viscosity
* `viscous` — enable viscous fluxes
* `LES` — enable SGS; `sgs_model = 'Smagorinsky' | 'WALE' | 'sigma_model' | 'MixedTimeScale'`; `yoshizawa`

> **Schemes & Sensors**

* `scheme` — `energy_preserving | hybrid_wenoEP | hybrid_tenoEP | hybrid_tenoaEP`
* `ep_type` — `JP | KEEP | AEC | PEP` (`PEP` requires the thermally perfect gas)
* `fd_order` — central order (2, 4, 6); `weno_order` — WENO-Z (3, 5, 7), TENO (5, 7), TENO-A (5)
* `diffusion_scheme` — `staggered | laplacian | conservative`
* `sensor` — `density | density_jump | ducros`; `ducrosToll` — Ducros threshold
* `enable_limiter` — optional positivity limiter

> **Boundary & Initial Conditions**

* `bc(1..6)` — one per face (W, E, S, N, B, F): `periodic`, `neumann`, `extrapolation`, `adiabatic_wall`, `isothermal_wall`, `diabatic_wall` (wall at `Trat` times the recovery temperature), `nscbc_inflow_generalized`, `nscbc_outflow`
* `inflow_profile = 'TurbulentBoundaryLayerInflow'`, `turb_inflow` — digital-filter turbulent inflow
* `ic` — `ambient`, `TGV`, `KolmogorovFlow`, `periodic_euler_x/y/z`, `shock_tube_x/y/z`, `isentropic_vortex_x/y`, `turbulent_channel`, `turbulent_BL`

> **I/O & Decomposition**

* `cart_dims(1:3)` — MPI Cartesian split (X, Y, Z); entries left at `0` are chosen to minimize the halo surface
* `output_file_name`, `data_dir`
* `itout` — field/restart cadence; `rsOut` — residual-monitor cadence; `StFlg`, `StOut` — statistics on/off and cadence

---

## Outputs

URANOS writes everything under `DATA/<data_dir>/`:

* `BINARY/<output_file_name><iteration>.bin` — conservative fields (MPI-IO), also used as restart files
* `RESIDUALS_MONITOR.txt` — L2 norm of the right-hand side and field extrema
* `BINARY_STAT/` — running statistics, when `StFlg = .true.`
* case diagnostics such as `TGV_MONITOR.txt` (Taylor–Green vortex) or `TCH_MONITOR.txt`, `VEL_STATS`, `BUD_STATS` (turbulent channel)
* `BINARY_DF/` — digital-filter state, so that a restart continues the same inflow turbulence

The post-processing scripts shipped with each case in `DATA/` show how to read these files.

---

## Performance

On one NVIDIA A100 the reference hypersonic boundary layer (50 million nodes) advances at about 0.76 s per iteration. On the Leonardo system the weak-scaling efficiency stays at 88–98 % up to 256 nodes (1024 A100). See the URANOS-3.0 paper and the CPC articles for the detailed performance analyses.

---

## Citing URANOS

If you use URANOS in academic work, please cite:

* **URANOS-3.0** — F. De Vanna, S. Regazzo, A. Aiello, S. Baaziz, G. Baldan, G. Coppola, E. Benini, *URANOS-3.0: GPU-accelerated structure-preserving solver for hypersonic flows on curvilinear/body-fitted grids* (DOI to be added).
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

We welcome issues and PRs! Please open an Issue or Discussion to propose changes before submitting a PR, and include testing notes in your PR.

---

## License

This project is released under the BSD 3-Clause license. See [`LICENSE`](./LICENSE) for details.

---

## Acknowledgments

Developed at the **Department of Industrial Engineering, University of Padova**, with runs on EuroHPC systems (Leonardo). We thank collaborators and HPC centers for support.

---

## FAQ / Troubleshooting

**Q1. The run stops with "The number of processors specified is not compatible".**
A: All three `cart_dims` entries are set and their product differs from the number of MPI ranks. Launch with a matching `-np`, or set some entries to `0` and let URANOS choose them.

**Q2. How many MPI ranks should I use on GPUs?**
A: One rank per GPU. For multi-node runs, set `--ntasks-per-node` to the number of GPUs per node and use a CUDA-aware MPI.

**Q3. Which scheme should I use for shocked flows?**
A: A hybrid scheme (`hybrid_wenoEP`, `hybrid_tenoEP`, `hybrid_tenoaEP`) with a shock sensor; pure `energy_preserving` fluxes are meant for smooth flows.

**Q4. Typical CFL values?**
A: The paper cases use `CFL` between 0.3 and 0.8, with 0.5 as a safe start; the solver warns when `CFL` exceeds the stability limit of the Runge–Kutta scheme for the chosen `fd_order`.

**Q5. Where do the statistics go?**
A: In `DATA/<data_dir>/BINARY_STAT/` (running means) and, for the channel and the Taylor–Green vortex, in the monitor files listed in [Outputs](#outputs).
