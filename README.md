<img src="images/logo6.svg" alt="image" width="100%" height="auto">

# Reference

Please, while using URANOS for your computations cite the following papers:

- Francesco De Vanna, Giacomo Baldan (2024). URANOS-2.0: Improved performance, enhanced portability, and model extension towards exascale computing of high-speed engineering flows. Under review in Computer Physics Communications.

- De Vanna, F., Avanzi, F., Cogo, M., Sandrin, S., Bettencourt, M., Picano, F., & Benini, E. (2023). URANOS: A GPU accelerated Navier-Stokes solver for compressible wall-bounded flows. Computer Physics Communications, 108717, doi: https://doi.org/10.1016/j.cpc.2023.108717

# Scientific papers obtained with URANOS

- Francesco De Vanna, Giacomo Baldan, Francesco Picano, Ernesto Benini, On the coupling between wall-modeled LES and immersed boundary method towards applicative compressible flow simulations, Computer & Fluids, 2023, https://doi.org/10.1016/j.compfluid.2023.106058

- Francesco De Vanna, Giacomo Baldan, Francesco Picano, Ernesto Benini, Effect of convective schemes in wall-resolved and wall-modeled LES of compressible wall turbulence, Computer & Fluids, 2022, https://doi.org/10.1016/j.compfluid.2022.105710

- Francesco De Vanna, Matteo Bernardini, Francesco Picano, Ernesto Benini, Wall-modeled LES of shock-wave/boundary layer interaction, International Journal of Heat and Fluid Flow, 2022, https://doi.org/10.1016/j.ijheatfluidflow.2022.109071

- Francesco De Vanna, Francesco Picano, Ernesto Benini, Mark Kenneth Quinn, Large- Eddy-Simulations of the unsteady behaviour of a hypersonic intake at Mach 5, AIAA Journal, 2021, https://doi.org/10.2514/1.J060160

- Francesco De Vanna, Michele Cogo, Matteo Bernardini, Francesco Picano, Ernesto Benini, Unified wall-resolved and wall-modelled method for large-eddy simulations of compressible wall-bounded flows, Physical Review Fluids, 2021, https://doi.org/10.1103/PhysRevFluids.6.034614

- Francesco De Vanna, Alberto Benato, Francesco Picano, Ernesto Benini, High order conservative formulation of viscous terms for variable viscosity flows, Acta Mechanica, 2021, https://doi.org/10.1007/s00707-021-02937-2

- Francesco De Vanna, Francesco Picano, Ernesto Benini, A sharp-interface immersed boundary method for moving objects in compressible viscous flows, Computers & Fluids, 2020, https://doi.org/10.1016/j.compfluid.2019.104415


# Compiling

URANOS requires (1) a Fortran compiler and (2) an MPI library. For the GPU version, the NVIDIA compiler. A Makefile with
predefined settings is available to facilitate compiling.
Different compiling modes can be selected by changing "make" options. The following is the syntax

```
make -j comp="option1" mode="option2"
```

`comp` currently supports these choices:
- `gnu`: GNU compiler
- `nvhpc`: NVHPC compiler (NVIDIA GPUs)
- `cray`: Cray compiler (AMD GPUs)

`"option2"` specifies the compiler options. Each compiler has dedicated options

Example:

```
make -j comp=gnu mode=cpu_debug
```

compiles the code with the gnu compiler using the debugging flags.

```
make -j comp=nvhpc mode=gpu
```

!CAUTION! Each cluster requires specific module to be loaded before the compiling process. 
You are address to cluster-speficic documentation concerning such issue. 

# Running

URANOS can be executed with a standard MPI launcher, e.g. `mpirun`.
In addition two files have to specified:
* `file.dat`: file defining the physical and numerical setup of the simulation, to be customized
according to the desired needs. Examples of input.dat files are available in the `tests` folder.
* `restart.bin`: only required for restart a simulation from previous results (optional)

Thus, to run a simulation without restart, type, e.g.:
```
mpirun -np "number_of_procs" ./Uranos.exe path/to/file.dat
```
if you want to restart a the simulation from previous data just type:
```
mpirun -np "number_of_procs" ./Uranos.exe path/to/file.dat path/to/restart.bin
```

For GPU runs you must distribute MPI processes according to the number of GPUs
available for each node. E.g., for CINECA Leonardo cluster -- 4 GPUS per node --  a possible submission script using
8 GPUs is:

```
#!/bin/bash
#SBATCH -p boost_usr_prod
#SBATCH --time 24:00:00
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8
module load nvhpc/23.1

mpirun -n $SLURM_NTASKS ./Uranos.exe path/to/file.dat path/to/restart.bin
```

For people not disposing of a HPC architecture, the solver can be compiled and run within the nvhpc framework by installing the suite from `https://developer.nvidia.com/hpc-sdk-downloads`.

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



