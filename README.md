# References

De Vanna, F., Avanzi, F., Cogo, M., Sandrin, S., Bettencourt, M., Picano, F., Benini, E. (2022) URANOS: a GPU accelerated Navier-Stokes solver for compressible wall-bounded flows, under review Computer Physics Communications.



# Compiling

URANOS requires (1) a Fortran compiler and (2) an MPI library. For the GPU version, the NVIDIA compiler. A Makefile with
predefined settings is available to facilitate compiling.
Different compiling modes can be selected by changing "make" options. The following is the syntax

```
make -j comp="option1" mode="option2"
```

`comp` currently supports these choices:
- `gnu`: GNU compiler
- `gnuch`: GNU compiler in CH version
- `intel`: INTEL compiler
- `pgi`: NVIDIA-PGI compiler

`"option2"` specifies the compiler options. Each compiler has dedicated options

Example:

```
make -j comp=intel model=debug
```

compiles the code with the intel compiler using the debugging flags.

# Running

URANOS can be executed with a standard MPI launcher, e.g. `mpirun`.
In addition two files have to specified:
* `file.dat`: file defining the physical and numerical setup of the simulation, to be customized
according to the desired needs. Examples of input.dat files are available in the `test` folder.
* `restart.bin`: only required for restart a simulation from previous results (optional)

Thus, to run a simulation without restart, type, e.g.:
```
mpirun -np "number_of_procs" ./Uranos.exe path/to/file.dat
```
if you want to restart a the simulation from previous data just type:
```
mpirun -np "number_of_procs" ./Uranos.exe path/to/file.dat path/to/restart.bin
```


