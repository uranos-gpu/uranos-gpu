<img src="images/logo6.svg" alt="image" width="100%" height="auto">

# References

When using **URANOS** for your computations, please cite the following papers:

- Francesco De Vanna, Giacomo Baldan (2024). *URANOS-2.0: Improved performance, enhanced portability, and model extension towards exascale computing of high-speed engineering flows*. Computer Physics Communications, 2024. https://doi.org/10.1016/j.cpc.2024.109285

- Francesco De Vanna, Filippo Avanzi, Michele Cogo, Simone Sandrin, Matt Bettencourt, Francesco Picano, Ernesto Benini (2023). *URANOS: A GPU accelerated Navier–Stokes solver for compressible wall-bounded flows*. Computer Physics Communications, 2023. https://doi.org/10.1016/j.cpc.2023.108717

# Scientific Papers Obtained with URANOS

- Francesco De Vanna, Ernesto Benini (2025). *Wall-Modeled LES of a Transonic Gas Turbine Vane – Part II: Mach Number Effect and Losses Prediction*. Journal of Turbomachinery, 2025.

- Francesco De Vanna, Ernesto Benini (2025). *Wall-Modeled LES of a Transonic Gas Turbine Vane – Part I: Model Setup and Assessment of Turbulent Length Scales*. Journal of Turbomachinery, 2025. https://doi.org/10.1115/1.4069131

- Francesco De Vanna (2025). *Entropy losses in transonic shock–boundary-layer interaction*. Physics of Fluids, 2025. https://doi.org/10.1063/5.0278759

- Francesco De Vanna, Ernesto Benini (2025). *Impact of Wall Cooling on Transonic Gas Turbine Stators Aerothermodynamics: Insights from Wall-Modeled LES*. Applied Thermal Engineering, 2025. https://doi.org/10.1016/j.applthermaleng.2025.126396

- Francesco De Vanna, Giacomo Baldan, Francesco Picano, Ernesto Benini (2023). *On the coupling between wall-modeled LES and immersed boundary method towards applicative compressible flow simulations*. Computers & Fluids, 2023. https://doi.org/10.1016/j.compfluid.2023.106058

- Francesco De Vanna, Giacomo Baldan, Francesco Picano, Ernesto Benini (2022). *Effect of convective schemes in wall-resolved and wall-modeled LES of compressible wall turbulence*. Computers & Fluids, 2022. https://doi.org/10.1016/j.compfluid.2022.105710

- Francesco De Vanna, Matteo Bernardini, Francesco Picano, Ernesto Benini (2022). *Wall-modeled LES of shock-wave/boundary layer interaction*. International Journal of Heat and Fluid Flow, 2022. https://doi.org/10.1016/j.ijheatfluidflow.2022.109071

- Francesco De Vanna, Francesco Picano, Ernesto Benini, Mark Kenneth Quinn (2021). *Large-Eddy Simulations of the unsteady behaviour of a hypersonic intake at Mach 5*. AIAA Journal, 2021. https://doi.org/10.2514/1.J060160

- Francesco De Vanna, Michele Cogo, Matteo Bernardini, Francesco Picano, Ernesto Benini (2021). *Unified wall-resolved and wall-modeled method for large-eddy simulations of compressible wall-bounded flows*. Physical Review Fluids, 2021. https://doi.org/10.1103/PhysRevFluids.6.034614

- Francesco De Vanna, Alberto Benato, Francesco Picano, Ernesto Benini (2021). *High-order conservative formulation of viscous terms for variable viscosity flows*. Acta Mechanica, 2021. https://doi.org/10.1007/s00707-021-02937-2

- Francesco De Vanna, Francesco Picano, Ernesto Benini (2020). *A sharp-interface immersed boundary method for moving objects in compressible viscous flows*. Computers & Fluids, 2020. https://doi.org/10.1016/j.compfluid.2019.104415


# Compiling

URANOS can be compiled for both CPU and GPU architectures.  
The compilation process requires the following dependencies:

1. **Fortran compiler** – such as `gfortran`, `nvfortran`, or the Cray Fortran compiler.
2. **MPI library** – such as OpenMPI, MPICH, or Cray MPI.
3. **(Optional, for GPU)** NVIDIA HPC SDK (`nvfortran`) or other supported GPU compilers.

A `Makefile` with predefined settings is provided to facilitate the compilation process.  
Different compiling modes can be selected by specifying `make` options.

---

## Syntax

```bash
make -j <nproc> comp="<compiler>" mode="<build_mode>"


# Running

URANOS runs with a standard MPI launcher (`mpirun`) and requires an input file (`file.dat`) defining the physical and numerical setup, plus an optional restart file (`restart.bin`) if resuming from previous results; examples of `file.dat` are available in the `tests` folder. On GPU systems, MPI ranks should match the number of GPUs per node, and the appropriate compiler/MPI/CUDA modules must be loaded according to the cluster’s documentation.

## Examples

### Basic (no restart)
mpirun -np <nprocs> ./Uranos.exe path/to/file.dat

### With restart
mpirun -np <nprocs> ./Uranos.exe path/to/file.dat path/to/restart.bin

### Local workstation (single GPU)
mpirun -np 1 ./Uranos.exe ./tests/flat_plate/input.dat

### Local workstation (multi-GPU, 2 ranks → 2 GPUs)
mpirun -np 2 ./Uranos.exe ./tests/flat_plate/input.dat

## SLURM Scripts

### CPU cluster (1 node, 32 ranks)
#!/bin/bash
#SBATCH -J uranos_cpu
#SBATCH -p cpu_partition
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

module purge
module load gcc/12.1.0 openmpi/4.1.5

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat

### GPU cluster (generic, 1 node, 4 GPUs → 4 ranks)
#!/bin/bash
#SBATCH -J uranos_gpu
#SBATCH -p gpu_partition
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

module purge
module load nvhpc/24.3 openmpi/4.1.5    # adjust to your site

mpirun -np ${SLURM_NTASKS} ./Uranos.exe path/to/file.dat

### CINECA Leonardo (4 GPUs/node, use 2 nodes → 8 GPUs)
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

### Restart on GPU cluster
#!/bin/bash
#SBATCH -J uranos_restart
#SBATCH -p gpu_partition
#SBATCH --time=06:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

module purge
module load nvhpc/24.3 openmpi/4.1.5

mpirun -np ${SLURM_NTASKS} ./Uranos.exe ./cases/mycase.dat ./results/restart.bin






# Basics tests

SHOCK TUBE

To become familiar with URANOS it is recommended to first launch a one-dimensional test consisting of a shock tube. The commands are as follows:
```
mpirun -np 4 ./Uranos.exe tests/shock_tube_x.dat
```
The run produces the temporal evolution of a shock tube which can be post-processed downstream using PostUranos.exe. The command is as follows:
```
./PostUranos.exe tests/shock_tube_x.dat
```

The process produces as autoput a `DATA/SHOCK_TUBE_1D` directory within which the results of the simulation and the post-treatments. In particular results can be visualized retrospectively in terms of line graphs (GNUPLOT) or two-dimensional fields (VTK).

CHANNEL DNS

As the use of the solver becomes more complex, the user is encouraged to try the `CHANNEL_DNS` test. The text consists of a turbulent channel flow with a bulk Mach number of 1.5. The following is the operating command using 4 computing units:

```
mpirun -np 4 ./Uranos.exe tests/channel_dns.dat
```

Statistics are produced while the solver is running and collected in `DATA/CHANNEL_DNS/VEL_STATS` and `DATA/CHANNEL_DNS/BUD_STATS` for the wall normal and the mmomentum budget statistics, respectively. Using PostUranos.exe allows to derive contours plots. 

BOUNDARY LAYER

Similarly to the channel flow test case, the third test case comprises of a hypersonic turbulent boundary layer with the lower-wall modeled according to a Wall-Modeled LES framework. The test can be run according to the following command:

```
mpirun -np 4 ./Uranos.exe tests/hypersonic_boundary_layer.dat
```

Again, wall normal statistics are collected while the code is running and saved in specific directories over some discrete stations of the boundary layer.


# Interpreting the file.dat file

`xmin` defines the x-left boundary of the domain

`xmax` defines the x-right boundary of the domain

`ymin` defines the y-left boundary of the domain

`ymax` defines the y-right boundary of the domain

`zmin` defines the z-left boundary of the domain

`zmax` defines the z-right boundary of the domain


`gridpoint_x` specifies the path/to/the/gridX/file

`gridpoint_y` specifies the path/to/the/gridY/file

`gridpoint_z` specifies the path/to/the/gridZ/file

`uniform` grid is readely available without requiring an external grid file

`dims` specifies the problem dimensions and could 2 or 3

`tmax` is the total simulation time (in code units). 

`itmax` is the maximum number of iterations

`nx`, `ny`, `nz` are the number of point along the X, Y, Z coordinates. 
Those must be consistent with input grids.

`cart_dims(1)`, `cart_dims(2)`, `cart_dims(3)` specify the number of procs along the
X, Y, and Z directions. Some problems (channel flow and boundary layer) split the domain only along X and Z, but not along Y. 

`CFL` is the Courant-Friedrichs-Lewy parameter. The implemented Runge-Kutta method is stable up to CFL = 1. For safaty margin keep less than 0.8. 

`Dt` is a fixed time step is you dont want adaptive time-stepping (logical\_CFL must be .false.)

`bc(1)` = boundary condition at the x-left side

`bc(2)` = boundary condition at the x-right side

`bc(3)` = boundary condition at the y-left side

`bc(4)` = boundary condition at the y-right side

`bc(5)` = boundary condition at the z-left side

`bc(6)` = boundary condition at the z-right side

bc\_module.f90 provide a list of several boundary conditions

`sponge(1)` activate a sponge zone for the x-left side

`sponge(2)` activate a sponge zone for the x-right side

`sponge(3)` activate a sponge zone for the y-left side

`sponge(4)` activate a sponge zone for the x-right side

`sponge(5)` activate a sponge zone for the z-left side

`sponge(6)` activate a sponge zone for the x-right side

`Trat` fix the wall-to-adiabatic temperature ratio (i.e., `Trat = 1` the wall is adibatic, `Trat < 1` the wall is cold, `Trat > 1` the wall is hot)

`inflow_profile` select the inflow profile (look at `inflow_module.f90`)

`smooth_inflow` provides a time gradually increasing inflow condition

`turb_inflow` adds turbulence at the inflow

`Reynolds` is the Reynolds number of the flow

`Mach` is the Mach number of the flow

`Prandtl` is the Mach number of the flow

`logical_CFL` is `.true.` for adaptive time-stepping, `.false.` for static time-step

`scheme` select the convective scheme among: 
- `energy_preserving`, central energy-stable scheme, not to be used with shocked flows
- `hybrid_wenoEP` WENO scheme hybridized with energy preserving
- `hybrid_tenoEP` TENO scheme hybridized with energy preserving
- `hybrid_tenoaEP` TENO-A scheme hybrided with energy preserving

`sensor` defines the shock-sensor among:
- `density` shock sensor based on density gradient
- `density_jump` shock sensor based on density jump in cell
- `ducros` shock sensor based on a modified version of Ducros sensor

`fd_order` defines the order of accurary of central schemes (2,4,6). Default 6
`weno_order` defines the order of accurary of weno schemes (3,5,7). Default 5

`LES` activates the SGS model
`sgs_model` select the  SGS model among: 
- `WALE`
- `Smagorinsky`
- `sigma_model`

`viscous` activates viscous flux computation

`itout` defines the output iteration of restart.bin files

`StOut` defines the output iteration of statistics.bin files

`StFlg` activates the output of statistics

`ic` defines the initial condition (see `src/ic_module.f90`)

`output_file_name` defines the prefix  of all output files

`data_dir` specifies the output directory




# Interpreting the outputs
As default, URANOS outputs a `DATA/data_dir/BINARY` directory where `.bin` files are stored. Each file represents a possible restart of a simulation.
Data in the `DATA/data_dir/BINARY` can be post-processed with `PostUranos.exe` via the following command

```
./PostUranos.exe path/to/file.dat (path/to/restart.bin optional)
```

PostUranos is a serial code which provide the contours plots and some video related to a simulation. You can find the visual results in `DATA/data_dir/VTK` (for 3D fields) and `DATA/data_dir/VTK2D` (for 2D fields)

# Contributing
We appreciate any contributions and feedback that can improve URANOS. If you wish to contribute to the tool, please get in touch with the maintainers or open an Issue in the repository / a thread in Discussions. Pull Requests are welcome, but please propose/discuss the changes in an linked Issue first.

# Licencing
Please refer to the licence file. 



