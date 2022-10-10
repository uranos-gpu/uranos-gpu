#!/bin/bash

# showing leakege in MPI
OMPI_MCA_mpi_show_handle_leaks=1	
export OMPI_MCA_mpi_show_handle_leaks

# pedentatic check of MPI functions
OMPI_MCA_mpi_param_check=1
export OMPI_MCA_mpi_param_check

#OpenMP num threads
export OMP_NUM_THREADS=1

