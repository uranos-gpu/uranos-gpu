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


