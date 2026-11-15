![logo](images/logo_2026_bis.png)



# URANOS — Exascaling the Future of Aerospace

[![License: BSD](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](#license)
[![Fortran](https://img.shields.io/badge/language-Fortran90-informational.svg)]()
[![GPU](https://img.shields.io/badge/GPU-OpenACC-success.svg)]()
[![MPI](https://img.shields.io/badge/Parallel-MPI-informational.svg)]()


**URANOS** (Unsteady Robust All-around Navier–StOkes Solver) is a **GPU-accelerated solver for compressible flows**, built for direct and large-eddy simulation of turbulent, high-speed and hypersonic flows of aerospace interest.

It solves the compressible Navier–Stokes equations in **generalized curvilinear coordinates** on body-fitted structured grids, with **kinetic-energy and entropy preserving** central fluxes, **high-order shock capturing**, calorically or **thermally perfect** gas models, and runs on CPUs and **NVIDIA GPUs** with **OpenACC + MPI**, from a laptop to hundreds of GPU nodes.

---

## Table of Contents

- [Highlights](#highlights)
- [Quick Start](#quick-start)
- [Build](#build)
- [Run](#run)
- [Writing an Input Deck](#writing-an-input-deck)
- [Input Reference](#input-reference)
- [Physical Model](#physical-model)
- [Numerical Methods](#numerical-methods)
- [Grids](#grids)
- [Boundary and Initial Conditions](#boundary-and-initial-conditions)
- [Parallel Execution](#parallel-execution)
- [Outputs](#outputs)
- [Verification and Validation Cases](#verification-and-validation-cases)
- [Code Layout](#code-layout)
- [Performance](#performance)
- [Citing URANOS](#citing-uranos)
- [Publications Using URANOS](#publications-using-uranos)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [FAQ / Troubleshooting](#faq--troubleshooting)

---

## Highlights

- **Body-fitted grids**: curvilinear in-plane grids extruded in the third direction, with metric terms that satisfy the discrete geometric identities to round-off.
- **Structure-preserving convection**: four kinetic-energy preserving flux families (`JP`, `KEEP`, `AEC`, `PEP`) with entropy-conserving and pressure-equilibrium-preserving properties, at 2nd, 4th or 6th order.
- **Shock capturing**: the central flux is blended with WENO-Z (3/5/7), TENO (5/7) or adaptive TENO-A (5) reconstructions in characteristic variables, switched by a shock sensor.
- **Real-gas thermodynamics**: calorically perfect gas, or thermally perfect gas with NASA polynomials for air, N₂ and O₂.
- **Three viscous operators**: Laplacian, staggered and a conservative formulation in which total energy telescopes exactly.
- **LES**: Smagorinsky, WALE, sigma and mixed-time-scale subgrid models.
- **Inflow turbulence**: digital-filter synthetic turbulence with non-reflecting characteristic boundaries.
- **Scalable**: MPI 3D domain decomposition, CUDA-aware halo exchange, 87–98 % weak-scaling efficiency up to 1024 NVIDIA A100; results are bitwise independent of the decomposition.

---

## Quick Start

```bash
# 1. get the code
git clone https://github.com/uranos-gpu/uranos-gpu.git
cd uranos-gpu

# 2. build (CPU)
make comp=gnu mode=cpu

# 3. run a paper case: the Sod shock tube with six shock-capturing reconstructions
python3 DATA/03_sod_shock_tube/run_sod.py
python3 DATA/03_sod_shock_tube/plot_sod.py      # -> sod_profiles.pdf
```

On a GPU machine, build with `make comp=nvhpc mode=gpu` and launch one MPI rank per GPU. To set up your own simulation, start from the [minimal deck](#writing-an-input-deck) below.

---

## Build

### Requirements

| | CPU build | GPU build |
|---|---|---|
| Compiler | `gfortran` | `nvfortran` (NVIDIA HPC SDK, OpenACC, curand) |
| MPI | OpenMPI or MPICH (`mpif90` wrapper) | CUDA-aware MPI (e.g. OpenMPI, HPC-X) |
| Make | GNU make | GNU make |

The Python scripts of the test cases need `numpy` and `matplotlib` (some figures also `scipy` and LaTeX).

### Build modes

```bash
make comp=<compiler> mode=<mode> [gpu_cc=<ccXX>]
```

| `comp` | `mode` | purpose |
|---|---|---|
| `gnu` | `cpu` | optimized CPU build (`-O3`) |
| `gnu` | `cpu_debug` | bounds and runtime checks, warnings, core dump |
| `gnu` | `cpu_profiling` | `gprof` instrumentation |
| `nvhpc` | `gpu` | optimized OpenACC build for NVIDIA GPUs |
| `nvhpc` | `gpu_debug` | GPU build with line info and floating-point traps |
| `nvhpc` | `gpu_profiling` | GPU build with NVTX ranges for Nsight Systems |

- `gpu_cc` sets the target compute capability (default `ccnative`, i.e. the GPU of the build node). Set it explicitly, e.g. `gpu_cc=cc80` for A100 or `gpu_cc=cc90` for H100, when you compile on a login node without the target GPU.
- The Makefile compiles in parallel on all available cores; objects go to `.objs/`, modules to `.mods/`, the executable is `Uranos.exe`.
- `make clean` removes all build artifacts.
- The working precision is double (`rp = dp` in `src/parameters_module.f90`).

---

## Run

```bash
mpirun -np <nprocs> ./Uranos.exe <deck.dat> [<restart>]
```

- Run from the **repository root**: the solver writes to `DATA/<data_dir>/`, and paths in a deck (`mesh_file`) are relative to the current directory.
- `<restart>` is either a field file (`DATA/<data_dir>/BINARY/<name><it>.bin`) or the keyword `restart`, which picks the latest field file in `DATA/<data_dir>/BINARY/`. The run continues from the time and iteration stored in the file up to `itmax`/`tmax`.
- On GPU nodes, use one MPI rank per GPU; ranks are bound to the GPUs of their node automatically.
- `source .environment.sh` sets `OMP_NUM_THREADS=1` and the OpenMPI handle/parameter checks used during development.

### SLURM examples

**Generic GPU node (4 GPUs → 4 ranks)**

```bash
#!/bin/bash
#SBATCH -J uranos
#SBATCH -p <gpu_partition>
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

module purge
module load nvhpc openmpi        # site-specific names and versions

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/deck.dat
```

**CINECA Leonardo (4 A100 per node, 2 nodes → 8 GPUs)**

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

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/deck.dat
```

**Generic CPU node (32 ranks)**

```bash
#!/bin/bash
#SBATCH -J uranos_cpu
#SBATCH -p <cpu_partition>
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

module purge
module load gcc openmpi

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/deck.dat
```

---

## Writing an Input Deck

A deck is a Fortran namelist `&input_list ... /`. Comments start with `!`, strings are quoted, arrays can be given element by element (`bc(1) = 'periodic'`) or as a list.

This complete deck advects an isentropic vortex on a periodic 2D box; it runs in under a second on two CPU cores:

```fortran
&input_list
  ! domain and grid
  xmin = -5.0,  xmax = 5.0
  ymin = -5.0,  ymax = 5.0
  zmin =  0.0,  zmax = 1.0
  nx = 64, ny = 64, nz = 1
  dims = 2
  gridpoint_x = 'uniform', gridpoint_y = 'uniform', gridpoint_z = 'uniform'
  mesh_file   = ''
  cart_dims   = 0, 0, 1

  ! time integration
  logical_CFL = .true.
  CFL   = 0.5
  Dt    = 1.0e-3
  n_step = 3
  tmax  = 10.0
  itmax = 200

  ! physics
  Mach     = 0.5
  Reynolds = 1000.0
  Prandtl  = 0.71
  thermally_perfect_gas = 0
  viscous  = .false.
  LES      = .false.
  sgs_model = 'WALE'

  ! numerics
  scheme     = 'energy_preserving'
  ep_type    = 'KEEP'
  fd_order   = 6
  weno_order = 5
  sensor     = 'ducros'

  ! boundary and initial conditions
  bc = 'periodic', 'periodic', 'periodic', 'periodic', 'periodic', 'periodic'
  inflow_profile = 'none'
  ic = 'isentropic_vortex_x'

  ! output
  output_file_name = 'vortex'
  data_dir = 'my_first_run'
  itout = 100
  rsOut = 10
  StFlg = .false.
  StOut = 0
/
```

```bash
mpirun -np 2 ./Uranos.exe vortex.dat
```

The run prints one line every `itout` iterations (iteration, percentage of `tmax`, time step, maximum Mach number, extrema of density, temperature and pressure, wall time per iteration) and writes `DATA/my_first_run/`. The decks in `DATA/` are complete, tested examples for every feature listed below.

---

## Input Reference

Keys marked *required* have no default and must be set in every deck.

### Domain, grid and decomposition

| key | default | description |
|---|---|---|
| `xmin, xmax, ymin, ymax, zmin, zmax` | *required* | domain extents |
| `nx, ny, nz` | *required* | number of nodes; in 2D `nz` is set to 1 |
| `dims` | *required* | `2` or `3` |
| `gridpoint_x/y/z` | *required* | `'uniform'` |
| `mesh_file` | `''` | curvilinear in-plane grid file (see [Grids](#grids)); overrides `gridpoint_x/y` |
| `cart_dims(1:3)` | *required* | MPI ranks per direction; `0` entries are chosen by the solver |

### Time integration

| key | default | description |
|---|---|---|
| `logical_CFL` | `.false.` | `.true.`: adaptive Δt from `CFL` and `Fourier`; `.false.`: fixed `Dt` |
| `CFL` | *required* | convective stability number (0.3–0.8 in the paper cases) |
| `Fourier` | `0.1` | diffusive stability number, used when `viscous = .true.` |
| `Dt` | *required* | time step when `logical_CFL = .false.` (initial step otherwise) |
| `n_step` | *required* | Runge–Kutta stages, `3` |
| `itmax`, `tmax` | *required* | iteration and time limits (the first reached stops the run) |

### Physics

| key | default | description |
|---|---|---|
| `Mach`, `Reynolds`, `Prandtl` | *required* | reference Mach, Reynolds and Prandtl numbers (meaning of `Reynolds` per initial condition [below](#initial-conditions-ic)) |
| `thermally_perfect_gas` | `0` | `0` calorically perfect (γ = 1.4), `1` thermally perfect |
| `species_name` | `'air'` | `'air'`, `'N2'`, `'O2'` (thermally perfect gas) |
| `Tref` | `273.15` | reference (free-stream) temperature in K; sets the thermally perfect tables and Sutherland's law |
| `vis_flag` | `1` | `0` power law, `1` Sutherland, `2` constant viscosity |
| `viscous` | `.false.` | enable the viscous fluxes |
| `diffusion_scheme` | `'staggered'` | `'staggered'`, `'laplacian'`, `'conservative'` |
| `LES` | `.false.` | enable a subgrid model |
| `sgs_model` | *required* | `'Smagorinsky'`, `'WALE'`, `'sigma_model'`, `'MixedTimeScale'` |
| `yoshizawa` | `.true.` | add the Yoshizawa isotropic subgrid term (LES only) |
| `Trat` | `1.0` | wall-to-recovery temperature ratio for `diabatic_wall` |
| `kolmogorov_forcing` | `0.1` | body-force amplitude for `ic = 'KolmogorovFlow'` |

### Convective scheme and shock capturing

| key | default | description |
|---|---|---|
| `scheme` | *required* | `'energy_preserving'` (central only), `'hybrid_wenoEP'`, `'hybrid_tenoEP'`, `'hybrid_tenoaEP'`, `'none'` (no convection) |
| `ep_type` | `'JP'` | `'JP'`, `'KEEP'`, `'AEC'`, `'PEP'` (`PEP` needs `thermally_perfect_gas = 1`) |
| `fd_order` | *required* | order of the central flux: `2`, `4`, `6` |
| `weno_order` | `5` | WENO-Z `3/5/7`, TENO `5/7`, TENO-A `5` |
| `sensor` | *required* | `'ducros'`, `'density'`, `'density_jump'`, `'weno'` (shock capturing everywhere) |
| `ducrosToll` | `0.1` | threshold of the Ducros sensor |
| `enable_limiter` | `.false.` | positivity limiter on density and internal energy |
| `limiter_rho_floor`, `limiter_e_floor` | `1e-3`, `1e-4` | floors of the limiter |
| `limiter_abort_frac` | `1e-4` | stop the run when a larger fraction of cells is limited (`≤ 0` disables) |

### Boundary and initial conditions

| key | default | description |
|---|---|---|
| `bc(1:6)` | *required* | one per face: W (x min), E (x max), S (y min), N (y max), B (z min), F (z max) |
| `ic` | *required* | initial condition |
| `inflow_profile` | — | `'TurbulentBoundaryLayerInflow'`; read only when `bc(1) = 'nscbc_inflow_generalized'` |
| `turb_inflow` | `.false.` | add digital-filter turbulence to the inflow profile |

### Output

| key | default | description |
|---|---|---|
| `data_dir` | *required* | output folder, `DATA/<data_dir>/` |
| `output_file_name` | *required* | prefix of the field files |
| `itout` | *required* | field/restart output and screen-print interval |
| `rsOut` | `10` | interval of `RESIDUALS_MONITOR.txt` |
| `StFlg`, `StOut` | `.true.`, `0` | statistics on/off and sampling interval (`StOut > 0` needed) |
| `restartStat` | `.false.` | continue the running statistics of a restarted run |
| `printres` | `.true.` | `.false.` disables all output (timing runs) |
| `itPrb` | `0` | probe sampling interval (`0` = off), see [Outputs](#outputs) |
| `itWallMap` | `0` | wall-pressure map interval (`0` = off) |

### Parallel communication

| key | default | description |
|---|---|---|
| `mpi_opt_level` | `2` | halo exchange: `1` blocking send-receive, `2` non-blocking with buffers (recommended), `3` MPI derived datatypes |
| `cuda_aware` | `.true.` | pass device buffers directly to MPI (needs a CUDA-aware MPI on GPUs) |

---

## Physical Model

URANOS solves the compressible Navier–Stokes equations for the conservative variables ρ, ρu, ρv, ρw, ρE, closed by the ideal-gas law p = ρRT.

- **Calorically perfect gas**: constant cp, γ = 1.4.
- **Thermally perfect gas**: cp(T) from NASA polynomials of the selected species, with the temperature recovered from the internal energy by Newton iterations. Below the validity range of the polynomials a smooth quartic tail keeps cp physical.
- **Transport**: viscosity from a power law, Sutherland's law (default) or a constant; thermal conductivity from a constant Prandtl number.

**Non-dimensional form.** Density, temperature and pressure are scaled by their free-stream values, so the free stream has ρ = T = p = 1, and velocities by √(R T_ref), so the free-stream velocity is u∞ = √γ_ref · Mach. All quantities in the outputs are in these units, and `time` in the field files is the solver time.

---

## Numerical Methods

**Convective fluxes.** Central, locally conservative fluxes in split form, of order `fd_order`:

| `ep_type` | properties | gas |
|---|---|---|
| `JP` | kinetic-energy preserving (Jameson–Pirozzoli) | CP and TP |
| `KEEP` | kinetic-energy and entropy preserving (truncated asymptotic expansion) | CP and TP |
| `AEC` | asymptotically entropy conservative and pressure-equilibrium preserving; the solver picks the CP or TP form | CP and TP |
| `PEP` | pressure-equilibrium preserving, derived for thermally perfect gases | TP only |

**Shock capturing.** With a `hybrid_*` scheme, the sensor flags the interfaces near discontinuities; there the central flux is replaced by a Lax–Friedrichs split flux reconstructed in characteristic variables with WENO-Z, TENO or TENO-A. Elsewhere the central flux is used unchanged.

**Viscous fluxes.** `staggered` (default) and `laplacian` use compact second-derivative stencils; `conservative` builds the full stress tensor and heat flux on cell faces so that the kinetic-energy and internal-energy exchanges cancel and total energy is conserved exactly; it keeps its formal order with variable viscosity on curvilinear grids.

**Subgrid models.** Algebraic eddy-viscosity models (Smagorinsky, WALE, sigma, mixed time scale) with an optional Yoshizawa isotropic term; with `LES = .false.` the scheme runs as DNS or implicit LES.

**Time integration.** Three-stage strong-stability-preserving Runge–Kutta (third order). With `logical_CFL = .true.` the step is the minimum of the convective (`CFL`) and diffusive (`Fourier`) limits over the domain, evaluated every step; the solver warns when `CFL` exceeds the stability limit of the chosen `fd_order`.

---

## Grids

- **Built-in uniform grid**: `gridpoint_x/y/z = 'uniform'` and `mesh_file = ''`.
- **Curvilinear grid**: `mesh_file` gives the in-plane (x, y) node coordinates; z stays uniform between `zmin` and `zmax`. The file is plain text: a header `nx ny nz`, then all x coordinates, then all y coordinates, one value per line, with `i` running fastest over `1-7 … nx+7` and `j` over `1-7 … ny+7` (seven ghost layers per side).
- `utils/mesh_to_uranos.py` converts a structured 2D O-grid (PLOT3D `*.x` or gmsh `*.msh`) to this format, including orientation and ghost-layer construction. The case folders contain generators for wavy, warped and compression-ramp grids.
- **2D runs**: `dims = 2`, one node in z, `cart_dims(3) = 1`.
- The metric terms are computed with the same central operator as the fluxes, so a uniform flow is preserved to round-off on any grid.

---

## Boundary and Initial Conditions

### Boundary conditions (`bc`)

| value | type |
|---|---|
| `periodic` | periodic (handled by the halo exchange; set it on both opposite faces) |
| `neumann` | zero-gradient |
| `extrapolation` | extrapolation (e.g. top boundary of a boundary layer) |
| `adiabatic_wall` | no-slip adiabatic wall |
| `isothermal_wall` | no-slip wall at the free-stream temperature |
| `diabatic_wall` | no-slip wall at `Trat` × recovery temperature |
| `nscbc_inflow_generalized` | non-reflecting characteristic inflow (W face), with the `inflow_profile` and optional digital-filter turbulence |
| `nscbc_outflow` | non-reflecting characteristic outflow (W/E/S/N faces) |

Walls are applied by mirroring into the ghost layers; the characteristic conditions act on the right-hand side at the boundary nodes.

### Initial conditions (`ic`)

| value | flow | notes |
|---|---|---|
| `ambient` | fluid at rest, ρ = p = 1 | |
| `isentropic_vortex_x`, `isentropic_vortex_y` | isentropic vortex advected along x or y | |
| `periodic_euler_x/y/z` | harmonic density wave advected along one direction | order-of-accuracy tests |
| `shock_tube_x/y/z` | Sod shock tube along one direction | |
| `TGV` | Taylor–Green vortex | writes `TGV_MONITOR.txt` with statistics on |
| `KolmogorovFlow` | fluid at rest driven by a sinusoidal body force | amplitude `kolmogorov_forcing` |
| `turbulent_channel` | perturbed turbulent channel, walls on S and N | `Reynolds` = target Re_τ, `Mach` = bulk Mach number; constant mass flow forcing |
| `turbulent_BL` | turbulent boundary layer with a wall on S | `Reynolds` = Re_τ at the inflow; use with `nscbc_inflow_generalized`, `TurbulentBoundaryLayerInflow` and `cart_dims(2) = 1` |

---

## Parallel Execution

- **Decomposition**: the domain is split into a 3D Cartesian grid of blocks. Entries of `cart_dims` left at `0` are chosen to minimize the halo surface; ranks sharing a node are given a compact tile of blocks.
- **Halo exchange**: four ghost layers, exchanged one direction after the other so that edge and corner ghosts are always current. On GPUs the buffers stay on the device and go directly to a CUDA-aware MPI.
- **Reproducibility**: the result is bitwise independent of the number of ranks and of the decomposition.
- **GPU offload**: the whole time step runs on the device with OpenACC; data move to the host only for output.

---

## Outputs

Everything is written to `DATA/<data_dir>/`:

| file / folder | content |
|---|---|
| `BINARY/<output_file_name><it>.bin` | conservative fields every `itout` iterations; also the restart files |
| `GRID/` | grid coordinates (`Xgrid.txt`, `Ygrid.txt`, `Zgrid.txt`) |
| `RESIDUALS_MONITOR.txt` | every `rsOut` iterations: time step, L2 norm of the right-hand side of each equation, field extrema |
| `BINARY_STAT/` | running statistics (`StFlg = .true.`, `StOut > 0`), averaged over the homogeneous directions detected from the grid and the periodic boundaries |
| `TGV_MONITOR.txt` | Taylor–Green vortex: kinetic energy, enstrophy, dissipation rates, entropy and total energy |
| `TCH_MONITOR.txt`, `VEL_STATS/`, `BUD_STATS/` | turbulent channel: wall quantities, velocity statistics, budgets |
| `BINARY_DF/` | digital-filter state, so that a restarted run continues the same inflow turbulence |
| `PROBES.txt` | probe time series (`itPrb > 0`) |
| `WALL_PRESSURE.bin` | wall-pressure map on the S wall (`itWallMap > 0`) |

**Field file format.** Little-endian binary: four 32-bit integers `nx, ny, nz, it`, two 64-bit reals `time, dt`, then the five conservative fields ρ, ρu, ρv, ρw, ρE as 64-bit reals, Fortran order (x fastest), interior nodes only. In Python:

```python
import struct, numpy as np

def read_field(path):
    with open(path, 'rb') as f:
        nx, ny, nz, it = struct.unpack('<4i', f.read(16))
        time, dt       = struct.unpack('<2d', f.read(16))
        q = np.frombuffer(f.read(5*nx*ny*nz*8), dtype='<f8').reshape(5, nz, ny, nx)
    return it, time, q          # q[0] = rho, q[1:4] = momentum, q[4] = rho E

it, t, q = read_field('DATA/my_first_run/BINARY/vortex0000200.bin')
rho = q[0]; u = q[1]/rho
```

**Probes.** Field probes (`prb_nfield`, `prb_x`, `prb_y`) record ρ, u, v, w, p, T at the node nearest to each point, averaged in z. Wall probes (`prb_nwall`, `prb_xc`, `prb_side`) record the pressure and skin-friction coefficients at a fraction `prb_xc` of the chord of a closed body wall on the S face (O-grid with periodic W/E faces). Up to 16 probes of each kind.

---

## Verification and Validation Cases

The `DATA/` folder holds the cases of the URANOS-3.0 paper. Each folder is self-contained: input decks, meshes, the scripts that generate them, and the scripts that run the campaigns and turn the results into the paper figures. See [`DATA/README.md`](./DATA/README.md) and the `README.md` of each folder.

| folder | case | what it verifies |
|---|---|---|
| `01_isentropic_vortex` | isentropic vortex, Cartesian and wavy grid | convective fluxes, free-stream preservation on curvilinear grids |
| `02_order_of_accuracy` | density wave and Kolmogorov flow | formal order of every convective, shock-capturing and viscous operator |
| `03_sod_shock_tube` | Sod shock tube, CP and TP gas | shock capturing, thermally perfect gas |
| `04_kolmogorov_flow` | steady Kolmogorov flow | viscous operators |
| `05_taylor_green_vortex` | Taylor–Green vortex at Re = 1600 | transition and decay, all fluxes and viscous schemes |
| `06_time_integration` | temporal order and stability | Runge–Kutta scheme, CFL and Fourier limits |
| `07_supersonic_compression_ramp` | 24° ramp at Ma = 2.91, body-fitted DNS | shock/boundary-layer interaction, digital-filter inflow, NSCBC |
| `08_multi_gpu_scaling` | triperiodic flow at rest | strong and weak scaling on many GPUs |

Folders 01, 03, 04 and 06 run on a CPU; the others are meant for GPUs.

---

## Code Layout

```
src/                      solver
  main.f90                time loop
  parameters_module.f90   namelist, constants, precision
  input_module.f90        input reading and checks
  mpi_module.f90          decomposition and rank placement
  mpi_comm_module.f90     halo exchange backends
  mesh_module.f90         grids and metric terms
  thermodynamics_module   CP/TP gas, NASA polynomials
  advection_module.f90    convective-flux dispatcher
  advection_schemes/      JP, KEEP, AEC, PEP central and hybrid kernels, WENO/TENO reconstruction
  viscous_module.f90      viscous-flux dispatcher
  viscous_schemes/        laplacian, staggered, conservative operators
  sgs_module.f90          subgrid models (WALE in sgs_schemes/)
  shock_detection_module  shock sensors
  bc_module.f90           wall, zero-gradient and extrapolation boundaries
  nscbc_boundary_conditions.f90   characteristic inflow/outflow
  inflow_module.f90, df_module.f90   inflow profile and digital filter
  ic_module.f90           initial conditions
  rhs_module.f90          right-hand side and body forces
  time_module.f90         time step and Runge-Kutta stages
  update_variable_module  primitive variables, positivity limiter
  output/                 restart, residuals, statistics, probes
libs/                     shared utilities (I/O, math, interpolation, statistics, random numbers)
utils/mesh_to_uranos.py   structured O-grid converter
DATA/                     paper cases
```

---

## Performance

- **One GPU**: a spatially developing hypersonic boundary layer (Ma = 5.86, 50 million nodes) advances at about 0.76 s per iteration on one NVIDIA A100 with the reference setup (JP flux, WENO-Z 5, staggered viscous scheme, WALE).
- **Many GPUs**: on CINECA Leonardo (4 A100 per node) the weak-scaling efficiency is 87–97 % at 256 nodes (1024 GPUs) and the strong-scaling efficiency 79–86 % at 64 nodes for the shock-capturing configurations.

The full analysis, including the cost of each numerical choice, is in the URANOS-3.0 paper.

---

## Citing URANOS

If you use URANOS in academic work, please cite the software papers:

* F. De Vanna, S. Regazzo, A. Aiello, S. Baaziz, G. Baldan, G. Coppola, E. Benini. *URANOS-3.0: GPU-accelerated structure-preserving solver for hypersonic flows on curvilinear/body-fitted grids* (DOI to be added).
* F. De Vanna, G. Baldan (2024). *URANOS-2.0: Improved performance, enhanced portability, and model extension towards exascale computing of high-speed engineering flows.* **Computer Physics Communications** 303, 109285. [doi:10.1016/j.cpc.2024.109285](https://doi.org/10.1016/j.cpc.2024.109285)
* F. De Vanna, F. Avanzi, M. Cogo, S. Sandrin, M. Bettencourt, F. Picano, E. Benini (2023). *URANOS: A GPU accelerated Navier-Stokes solver for compressible wall-bounded flows.* **Computer Physics Communications** 287, 108717. [doi:10.1016/j.cpc.2023.108717](https://doi.org/10.1016/j.cpc.2023.108717)

`CITATION.cff` is read by GitHub's "Cite this repository" button; BibTeX entries are available at the DOI pages.

---

## Publications Using URANOS

Journal articles whose results were obtained with URANOS, most recent first.

**2026**

* F. De Vanna, E. Benini. *A design-oriented IBM–WMLES framework for rapid parametric mapping of turbine-cascade wake aerothermodynamics.* **Aerospace Science and Technology** 176, 112625. [doi:10.1016/j.ast.2026.112625](https://doi.org/10.1016/j.ast.2026.112625)
* F. De Vanna. *Transonic film cooling for future aviation gas turbines: A high-fidelity large eddy simulations reference.* **Aerospace Science and Technology** 173, 111801. [doi:10.1016/j.ast.2026.111801](https://doi.org/10.1016/j.ast.2026.111801)

**2025**

* F. De Vanna, E. Benini. *Wall-Modeled Large-Eddy Simulation of a Transonic Gas Turbine Vane—Part I: Model Setup and Assessment of Turbulent Length Scales.* **Journal of Turbomachinery** 147(11), 111015. [doi:10.1115/1.4069131](https://doi.org/10.1115/1.4069131)
* F. De Vanna, E. Benini. *Wall-Modeled Large-Eddy Simulation of a Transonic Gas Turbine Vane—Part II: Mach Number Effect and Losses Prediction.* **Journal of Turbomachinery** 147(11), 111016. [doi:10.1115/1.4069384](https://doi.org/10.1115/1.4069384)
* F. De Vanna. *Entropy losses in transonic shock–boundary-layer interaction.* **Physics of Fluids** 37, 076159. [doi:10.1063/5.0278759](https://doi.org/10.1063/5.0278759)
* F. De Vanna, E. Benini. *Impact of wall cooling on transonic gas turbine stators aerothermodynamics: Insights from Wall-Modeled LES.* **Applied Thermal Engineering** 272, 126396. [doi:10.1016/j.applthermaleng.2025.126396](https://doi.org/10.1016/j.applthermaleng.2025.126396)

**2024**

* F. De Vanna, G. Baldan. *URANOS-2.0: Improved performance, enhanced portability, and model extension towards exascale computing of high-speed engineering flows.* **Computer Physics Communications** 303, 109285. [doi:10.1016/j.cpc.2024.109285](https://doi.org/10.1016/j.cpc.2024.109285)

**2023**

* F. De Vanna, G. Baldan, F. Picano, E. Benini. *On the coupling between wall-modeled LES and immersed boundary method towards applicative compressible flow simulations.* **Computers & Fluids** 266, 106058. [doi:10.1016/j.compfluid.2023.106058](https://doi.org/10.1016/j.compfluid.2023.106058)
* F. De Vanna, G. Baldan, F. Picano, E. Benini. *Effect of convective schemes in wall-resolved and wall-modeled LES of compressible wall turbulence.* **Computers & Fluids** 250, 105710. [doi:10.1016/j.compfluid.2022.105710](https://doi.org/10.1016/j.compfluid.2022.105710)
* F. De Vanna, F. Avanzi, M. Cogo, S. Sandrin, M. Bettencourt, F. Picano, E. Benini. *URANOS: A GPU accelerated Navier-Stokes solver for compressible wall-bounded flows.* **Computer Physics Communications** 287, 108717. [doi:10.1016/j.cpc.2023.108717](https://doi.org/10.1016/j.cpc.2023.108717)

**2022**

* F. De Vanna, M. Bernardini, F. Picano, E. Benini. *Wall-modeled LES of shock-wave/boundary layer interaction.* **International Journal of Heat and Fluid Flow** 98, 109071. [doi:10.1016/j.ijheatfluidflow.2022.109071](https://doi.org/10.1016/j.ijheatfluidflow.2022.109071)

**2021**

* F. De Vanna, F. Picano, E. Benini, M. K. Quinn. *Large-Eddy Simulations of the Unsteady Behavior of a Hypersonic Intake at Mach 5.* **AIAA Journal** 59(10), 3859–3872. [doi:10.2514/1.J060160](https://doi.org/10.2514/1.J060160)
* F. De Vanna, M. Cogo, M. Bernardini, F. Picano, E. Benini. *Unified wall-resolved and wall-modeled method for large-eddy simulations of compressible wall-bounded flows.* **Physical Review Fluids** 6, 034614. [doi:10.1103/PhysRevFluids.6.034614](https://doi.org/10.1103/PhysRevFluids.6.034614)
* F. De Vanna, A. Benato, F. Picano, E. Benini. *High-order conservative formulation of viscous terms for variable viscosity flows.* **Acta Mechanica** 232, 2115–2133. [doi:10.1007/s00707-021-02937-2](https://doi.org/10.1007/s00707-021-02937-2)

**2020**

* F. De Vanna, F. Picano, E. Benini. *A sharp-interface immersed boundary method for moving objects in compressible viscous flows.* **Computers & Fluids** 201, 104415. [doi:10.1016/j.compfluid.2019.104415](https://doi.org/10.1016/j.compfluid.2019.104415)

---

## Contributing

Issues and pull requests are welcome. Please open an Issue to discuss a change before submitting a PR, and describe in the PR how you tested it (for example, the `DATA/` cases you ran and how the results compare).

---

## License

URANOS is released under the BSD 3-Clause license. See [`LICENSE`](./LICENSE).

---

## Acknowledgments

Developed at the **Department of Industrial Engineering, University of Padova**, with simulations on the EuroHPC system Leonardo at CINECA. We thank our collaborators and the HPC centers for their support.

---

## FAQ / Troubleshooting

**The run stops with "The number of processors specified is not compatible".**
All three `cart_dims` entries are set and their product differs from the number of MPI ranks. Launch with a matching `-np`, or set some entries to `0` and let URANOS choose them.

**How many MPI ranks per GPU?**
One. For multi-node runs set `--ntasks-per-node` to the number of GPUs per node and use a CUDA-aware MPI; with a non-CUDA-aware MPI set `cuda_aware = .false.`.

**Which scheme for flows with shocks?**
A hybrid scheme (`hybrid_wenoEP`, `hybrid_tenoEP`, `hybrid_tenoaEP`) with the Ducros sensor. The pure `energy_preserving` fluxes are meant for smooth flows; for strong shocks also consider `enable_limiter = .true.`.

**Which time step?**
`logical_CFL = .true.` with `CFL = 0.5` is a safe start; the paper cases use 0.3–0.8 and `Fourier = 0.1`. The solver warns when `CFL` exceeds the stability limit of the chosen `fd_order`.

**The run stops with "restarting time is greater than total simulation time".**
The restart file was written after `tmax`. Increase `tmax` (and `itmax`) in the deck.

**Where are the statistics?**
In `DATA/<data_dir>/BINARY_STAT/` (running means over the homogeneous directions) and, for the channel and the Taylor–Green vortex, in the monitor files listed in [Outputs](#outputs).
