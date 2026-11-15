#Programe Name
MAINPROG = Uranos.exe

# Directories
DIRECTORIES = .objs .mods DATA
OBJDIR = .objs/
MODDIR = .mods/
SRCDIR = src/
SCHDIR = src/advection_schemes/
VISDIR = src/viscous_schemes/
SGSDIR = src/sgs_schemes/
OUTDIR = src/output/
LIBDIR = libs/

# Compiler
FC = mpif90

# Libraries
LIBS =

# DEFAULTS!!!
# -----------------------------
OPT          = #-O3 -mtune=native -mcmodel=medium -fdefault-real-8
DEBUG        = #-g -DDEBUG -pedantic -fbounds-check -fcheck=all -Wall -Waliasing -Wextra -fmax-errors=5 #-ffpe-trap=zero,overflow,underflow #-Werror
OPENMP       = #-fopenmp -DOPENMP -lgomp -lpthread
TIME 	     = 
MODULES      = #-J$(MODDIR)
PREPROC      = -cpp
INCMODULES   =

# Parallel build: detect core count (Linux: nproc, macOS: sysctl)
NPROC := $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
MAKEFLAGS += -j$(NPROC) --output-sync=target

# SMART COMPILING OPTIONS
# -----------------------------

ifeq ($(comp),gnu)
     FC = mpif90
     MODULES = -J$(MODDIR)
     ifeq ($(mode),cpu)
     	OPT   = -O3 -pipe
     else ifeq ($(mode),cpu_debug)
     	OPT   = -O3 -pipe
     	DEBUG = -g -DDEBUG -pedantic -fbounds-check -fcheck=all -Wall -Waliasing -Wextra -fmax-errors=5 -fdump-core
	#-finit-real=inf -ffpe-trap=invalid,zero,overflow -Wuninitialized
     	TIME    =
     else ifeq ($(mode),cpu_profiling)
     	OPT   = -O3 -pipe
     	DEBUG =
     	TIME  = -DPROFCPU
     endif
endif
ifeq ($(comp),gnuch)
     FC = mpif90.mpich
     MODULES = -J$(MODDIR)
     OPT     = -O3 -pipe
     TIME    =
     OPENMP  =
     ifeq ($(mode),debug)
     	DEBUG  = -g -DDEBUG -pedantic -fbounds-check -fcheck=all -Wall -Waliasing -Wextra -fmax-errors=5 -fdump-core
     else ifeq ($(mode),time)
     	TIME   = -DTIME
     else ifeq ($(mode),openmp)
     	OPENMP = -fopenmp -DOPENMP -lgomp -lpthread
     endif
endif

ifeq ($(comp),intel)
	FC = mpif90
	MODULES = -module $(MODDIR)
	ifeq ($(mode),debug)
		OPT    = -O3
		DEBUG  = -check bounds,uninit -g -fpe0 -traceback
		TIME   =
		OPENMP =
	else ifeq ($(mode),time)
		OPT    = -O3
		DEBUG  =
		TIME   = -DTIME
		OPENMP =
	else ifeq ($(mode),report)
		OPT    = -O3
		DEBUG  = -qopt-report -qopt-report-phase=vec
		TIME   =
		OPENMP =
	else ifeq ($(mode),knl)
		OPT    = -O3 -xMIC-AVX512
		DEBUG  =
		TIME   =
		OPENMP =
	else 
		OPT    = -O3
		DEBUG  =
		TIME   =
		OPENMP =
	endif
endif

ifeq ($(comp),nvhpc)
	FC = mpif90
	MODULES = -module $(MODDIR)
	# GPU compute capability: ccnative autodetects the local GPU (workstation).
	# On machines whose build node has no GPU (e.g. Leonardo login nodes),
	# override explicitly: make comp=nvhpc mode=gpu gpu_cc=cc80
	gpu_cc ?= ccnative
	ifeq ($(mode),gpu)
		OPT    = -O3 -acc=gpu -gpu=$(gpu_cc),deepcopy,safecache -Minline=maxsize:1000  -DNVIDIA
		DEBUG  = -cudalib=curand
		TIME   = 
		OPENMP =
	else ifeq ($(mode),gpu_debug)
                OPT    = -O3 -acc=gpu -gpu=$(gpu_cc),deepcopy,lineinfo,safecache,ptxinfo -Minline=maxsize:1000 -cudalib=curand -DNVIDIA
                DEBUG  = -g -Minfo=accel -Minit-real=snan -Minit-msg -Ktrap=inv
                TIME   = 
                OPENMP =
	else ifeq ($(mode),gpu_profiling)
                OPT    = -O3 -acc=gpu -gpu=$(gpu_cc),deepcopy,lineinfo,safecache,ptxinfo -Minline=maxsize:1000 -cudalib=curand,nvtx3 -DNVIDIA
                DEBUG  = 
                TIME   = -DNVTX
                OPENMP =
	endif
endif

ifeq ($(fast_pgi),true)
	FC      = mpif90
	OPT     = -O3 -r8 -mcmodel=medium
	OPENMP  = 
	DEBUG   =
	TIME    = 
	MODULES = -module $(MODDIR)
endif

ifeq ($(mpidex),buffer)
	MPIDEX = -DMPIBFR
endif

FFLAGS := $(PREPROC) $(OPENMP) $(TIME) $(MPIDEX) $(DEBUG) $(OPT)

# Flag-tracking stamp: rebuild all objects when FFLAGS changes, not on every Makefile edit
FLAGFILE := $(OBJDIR).flags
$(shell mkdir -p $(OBJDIR) $(MODDIR))
$(shell printf '%s' '$(FFLAGS)' | cmp -s - $(FLAGFILE) 2>/dev/null || printf '%s' '$(FFLAGS)' > $(FLAGFILE))

# Per-object extra flags (empty for every object except where overridden below).
EXTRA_FFLAGS :=
# WALE occupancy cap: compile ONLY wale_module.o with -gpu=maxregcount:128.
# On A100 this repacks compute_wale 162->128 registers, 0 spill, ~12->16 warps/SM
# (~-13.5% kernel, byte-identical).  nvfortran -gpu sub-flag -> only the OpenACC
# GPU builds; a no-op for CPU builds.
WALE_REGCAP :=
ifeq ($(comp),nvhpc)
    ifneq (,$(filter gpu gpu_debug gpu_profiling,$(mode)))
        WALE_REGCAP := -gpu=maxregcount:128
    endif
endif

# weno_reconstruction occupancy cap: compile ONLY weno_reconstruction.o with
# -gpu=maxregcount:$(wenorec_cap).  weno_reconstruction lives in its own file, so
# -Minline (intra-file only on nvfortran) does NOT inline it into the blended
# convective sweeps: it is a real ABI device call whose register count therefore
# CAPS the launch occupancy of every hybrid_wenoEP_* kernel.  Uncapped it needs
# 220 registers (= 4 gangs/SM) while the X and Z entry functions need only 176
# and 138, so the shared callee -- not the kernel -- is what limits them.
# Capping it repacks the callee and lets each sweep run at its own limit.
# Register allocation does not change arithmetic -> byte-identical.
wenorec_cap ?= 168
WENOREC_REGCAP :=
ifeq ($(comp),nvhpc)
    ifneq (,$(filter gpu gpu_debug gpu_profiling,$(mode)))
        WENOREC_REGCAP := -gpu=maxregcount:$(wenorec_cap)
    endif
endif

# Blended-KEEP occupancy cap: compile ONLY hybrid_wenoEP_KEEP.o with
# -gpu=maxregcount:$(keep_cap).  Uncapped, ptxas picks 176/248/138 for the csi/eta/zeta
# sweeps; 168 puts all three at the 6-gang band (65536/(168*64) = 6.09) at the cost of
# 84/188/32 B of spill confined to the rarely-taken WENO shock branch.  Measured on A100
# 128^3/10 iter over five decks (Sod x/y/z, TGV Cartesian, TGV wavy): convective
# -3.7% .. -9.2%, negative on every deck, and better than the weno-only cap everywhere.
# 176 is NOT a safe value (it lifts zeta 138 -> 176 and loses its gang: +1% .. +6%).
# Register allocation does not change arithmetic -> byte-identical (cmp, 5 vars, 5 decks).
# Set keep_cap=0 to disable if a build must carry zero spill in the sweeps themselves.
keep_cap ?= 168
KEEP_REGCAP :=
ifeq ($(comp),nvhpc)
    ifneq (,$(filter gpu gpu_debug gpu_profiling,$(mode)))
        ifneq ($(keep_cap),0)
            KEEP_REGCAP := -gpu=maxregcount:$(keep_cap)
        endif
    endif
endif

# Objects shared by all executables (solver + tools + libs)
COMMON_OBJS = $(addprefix $(OBJDIR), \
	  parameters_module.o \
	  file_module.o \
	  allocate_module.o \
	  mpi_module.o \
	  mpi_comm_module.o \
	  real_to_integer_module.o \
	  math_tools_module.o \
	  matrix_inversion_module.o \
	  norm_module.o \
	  statistics_module.o \
	  random_module.o \
	  interpolation_module.o \
	  reynolds_averaged_module.o \
	  nvtx.o \
	  profiling_module.o \
	  npy.o \
	  vtk_utils_module.o \
	  thermodynamics_module.o \
	  input_module.o \
	  mesh_module.o \
	  storage_module.o \
	  fluid_functions_module.o \
	  onlineStats.o \
	  df_module.o \
	  inflow_module.o \
	  shock_detection_module.o \
	  euler_matrices_module.o \
	  nscbc_boundary_conditions.o \
	  GetRetau_module.o \
	  energy_preserving_JP.o \
	  energy_preserving_KEEP.o \
	  energy_preserving_AEC2CP.o \
	  energy_preserving_AECTP.o \
	  energy_preserving_PEPTP.o \
	  weno_reconstruction.o \
	  hybrid_wenoEP_JP.o \
	  hybrid_wenoEP_KEEP.o \
	  hybrid_wenoEP_AEC2CP.o \
	  hybrid_wenoEP_AECTP.o \
	  hybrid_wenoEP_PEPTP.o \
	  advection_module.o \
	  viscous_shared.o \
	  viscous_staggered.o \
	  viscous_laplacian.o \
	  viscous_conservative.o \
	  viscous_module.o \
	  bc_module.o \
	  sgs_bc_module.o \
	  wale_module.o \
	  sgs_module.o \
	  output_probes.o \
	  rhs_module.o \
	  update_variable_module.o \
	  output_residuals.o \
	  ic_module.o \
	  init_bc_module.o \
	  time_module.o \
	  output_restart.o \
	  output_screen.o \
	  output_stats_channel.o \
	  output_stats_tgv.o \
	  output_module.o)

vpath %.f90 $(SRCDIR) $(SCHDIR) $(VISDIR) $(SGSDIR) $(OUTDIR) $(LIBDIR)

# make all
.PHONY: all
all : $(DIRECTORIES) $(MAINPROG)

# Create the directory $@ if this does not exist
$(DIRECTORIES) :
	mkdir -p $@
	
.PHONY: solver
solver: $(DIRECTORIES) $(MAINPROG)

# -------------------------------------------------------
# Compile Uranos main program
# -------------------------------------------------------
$(MAINPROG): $(COMMON_OBJS) $(SRCDIR)main.f90
	@echo "-------------------------------------------"
	@echo "Creating the executable: $(MAINPROG)"
	@echo "-------------------------------------------"
	$(FC) $(FFLAGS) $(SRCDIR)main.f90 $(MODULES) $(INCMODULES) $(COMMON_OBJS) -o $(MAINPROG) $(LIBS)

# -------------------------------------------------------
# Compile any .f90 found via vpath into .objs/
# -------------------------------------------------------
$(OBJDIR)%.o : %.f90 $(FLAGFILE)
	$(FC) $(FFLAGS) $(EXTRA_FFLAGS) -c $< -o $@ $(MODULES) $(LIBS)

# WALE kernel: occupancy-cap this one object (see WALE_REGCAP above).
$(OBJDIR)wale_module.o : EXTRA_FFLAGS := $(WALE_REGCAP)

# Shared WENO/TENO reconstruction: occupancy-cap this one object (see WENOREC_REGCAP above).
$(OBJDIR)weno_reconstruction.o : EXTRA_FFLAGS := $(WENOREC_REGCAP)

# Blended KEEP sweeps: occupancy-cap this one object (see KEEP_REGCAP above).
$(OBJDIR)hybrid_wenoEP_KEEP.o : EXTRA_FFLAGS := $(KEEP_REGCAP)

	
# Clean everything up
.PHONY: clean
clean: 
	@echo "---------------------------------------"
	@echo "Cleaning everything up"
	@echo "---------------------------------------"
	rm -rf *.exe *.o *.mod *.i *dSYM
	rm -rf *.optrpt
	rm -rf $(OBJDIR)
	rm -rf $(MODDIR)

#---------------------------------------------------------------------------
# Dependency chains (source files omitted: handled by vpath pattern rule)
#---------------------------------------------------------------------------
$(OBJDIR)allocate_module.o                   : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)file_module.o                       : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)mpi_module.o                        : $(addprefix $(OBJDIR), parameters_module.o file_module.o allocate_module.o)
$(OBJDIR)thermodynamics_module.o             : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o)
$(OBJDIR)mpi_comm_module.o                   : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o profiling_module.o)
$(OBJDIR)real_to_integer_module.o            : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)math_tools_module.o                 : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)matrix_inversion_module.o           : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)norm_module.o                       : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)statistics_module.o                 : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)random_module.o                     : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)interpolation_module.o              : $(addprefix $(OBJDIR), parameters_module.o real_to_integer_module.o matrix_inversion_module.o mpi_module.o)
$(OBJDIR)input_module.o                      : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o real_to_integer_module.o thermodynamics_module.o)
$(OBJDIR)mesh_module.o                       : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o allocate_module.o matrix_inversion_module.o)
$(OBJDIR)storage_module.o                    : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o mesh_module.o allocate_module.o)
$(OBJDIR)vtk_utils_module.o                  : $(addprefix $(OBJDIR), parameters_module.o file_module.o)
$(OBJDIR)profiling_module.o                  : $(addprefix $(OBJDIR), nvtx.o parameters_module.o file_module.o)
$(OBJDIR)fluid_functions_module.o            : $(addprefix $(OBJDIR), parameters_module.o storage_module.o interpolation_module.o thermodynamics_module.o)
$(OBJDIR)onlineStats.o                       : $(addprefix $(OBJDIR), parameters_module.o profiling_module.o statistics_module.o fluid_functions_module.o)
$(OBJDIR)df_module.o                         : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o fluid_functions_module.o random_module.o)
$(OBJDIR)inflow_module.o                     : $(addprefix $(OBJDIR), storage_module.o fluid_functions_module.o math_tools_module.o df_module.o)
$(OBJDIR)shock_detection_module.o            : $(addprefix $(OBJDIR), parameters_module.o storage_module.o fluid_functions_module.o mpi_comm_module.o)
$(OBJDIR)euler_matrices_module.o             : $(addprefix $(OBJDIR), parameters_module.o thermodynamics_module.o)
$(OBJDIR)nscbc_boundary_conditions.o         : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o mesh_module.o inflow_module.o euler_matrices_module.o)
$(OBJDIR)GetRetau_module.o                   : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o fluid_functions_module.o mesh_module.o thermodynamics_module.o)
# per-scheme energy-preserving modules (src/advection_schemes/, via vpath)
# Wired set (ep_id 1..5): JP, KEEP(2), AEC2-CP, AEC-TP, PEP (ep_type JP/KEEP/AEC/PEP).
# AEC-TP and PEP carry the
# thermodynamics_module dependency (their thermally-perfect caloric closure).
EP_DEPS = $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o)
$(OBJDIR)energy_preserving_JP.o              : $(EP_DEPS)
$(OBJDIR)energy_preserving_KEEP.o           : $(EP_DEPS)
$(OBJDIR)energy_preserving_AEC2CP.o          : $(EP_DEPS)
$(OBJDIR)energy_preserving_AECTP.o           : $(EP_DEPS) $(OBJDIR)thermodynamics_module.o
$(OBJDIR)energy_preserving_PEPTP.o           : $(EP_DEPS) $(OBJDIR)thermodynamics_module.o
EP_OBJS = $(addprefix $(OBJDIR), energy_preserving_JP.o energy_preserving_KEEP.o \
	  energy_preserving_AEC2CP.o energy_preserving_AECTP.o energy_preserving_PEPTP.o)
# blended (hybrid shock-capturing + energy-preserving) convective schemes.
# weno_reconstruction = shared WENO/TENO cell-boundary reconstruction, common to
# every blended sibling; hybrid_wenoEP_JP / _KEEP = the JP- / KEEP(2)-split variants.
$(OBJDIR)weno_reconstruction.o               : $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)hybrid_wenoEP_JP.o                  : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o euler_matrices_module.o weno_reconstruction.o)
$(OBJDIR)hybrid_wenoEP_KEEP.o               : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o euler_matrices_module.o weno_reconstruction.o)
$(OBJDIR)hybrid_wenoEP_AEC2CP.o              : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o euler_matrices_module.o weno_reconstruction.o)
$(OBJDIR)hybrid_wenoEP_AECTP.o               : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o euler_matrices_module.o weno_reconstruction.o thermodynamics_module.o)
$(OBJDIR)hybrid_wenoEP_PEPTP.o               : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o euler_matrices_module.o weno_reconstruction.o thermodynamics_module.o)
$(OBJDIR)advection_module.o                  : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o inflow_module.o nvtx.o nscbc_boundary_conditions.o hybrid_wenoEP_JP.o hybrid_wenoEP_KEEP.o hybrid_wenoEP_AEC2CP.o hybrid_wenoEP_AECTP.o hybrid_wenoEP_PEPTP.o) $(EP_OBJS)
# viscous-flux schemes (src/viscous_schemes/, via vpath) + orchestrator (src/).
# viscous_module = thin dispatcher over diffusion_scheme; the schemes it drives:
#   viscous_shared    = DivergencyGradient3D + SGS-isotropic (shared completion),
#   viscous_laplacian = chain-rule Laplacian form,
#   viscous_staggered = conservative Lele face-flux form,
#   viscous_conservative = compatible (energy-conserving) 2nd-order face form,
$(OBJDIR)viscous_shared.o                    : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mpi_comm_module.o mesh_module.o nvtx.o fluid_functions_module.o profiling_module.o)
$(OBJDIR)viscous_staggered.o                 : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o nvtx.o fluid_functions_module.o profiling_module.o)
$(OBJDIR)viscous_laplacian.o                 : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o nvtx.o fluid_functions_module.o profiling_module.o)
$(OBJDIR)viscous_conservative.o              : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mpi_comm_module.o mesh_module.o nvtx.o fluid_functions_module.o profiling_module.o)
$(OBJDIR)viscous_module.o                    : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o profiling_module.o nvtx.o viscous_shared.o viscous_staggered.o viscous_laplacian.o viscous_conservative.o)
$(OBJDIR)bc_module.o                         : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o fluid_functions_module.o inflow_module.o real_to_integer_module.o mpi_comm_module.o thermodynamics_module.o)
$(OBJDIR)sgs_bc_module.o                     : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o bc_module.o mesh_module.o fluid_functions_module.o real_to_integer_module.o matrix_inversion_module.o math_tools_module.o file_module.o profiling_module.o)
# WALE SGS kernel split into its own compilation unit (src/sgs_schemes/) so it
# can be compiled with -gpu=maxregcount:128 (see WALE_REGCAP below) without
# capping the sigma/MixedTimeScale kernels in sgs_module (194/234 regs, would spill).
$(OBJDIR)wale_module.o                       : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o mesh_module.o)
$(OBJDIR)sgs_module.o                        : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o bc_module.o sgs_bc_module.o wale_module.o)
$(OBJDIR)output_probes.o                     : $(addprefix $(OBJDIR), parameters_module.o profiling_module.o storage_module.o mpi_module.o mesh_module.o fluid_functions_module.o file_module.o)
$(OBJDIR)rhs_module.o                        : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o advection_module.o viscous_module.o sgs_module.o thermodynamics_module.o)
$(OBJDIR)update_variable_module.o            : $(addprefix $(OBJDIR), parameters_module.o storage_module.o profiling_module.o thermodynamics_module.o)
$(OBJDIR)output_residuals.o                  : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o profiling_module.o output_probes.o)
$(OBJDIR)ic_module.o                         : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o fluid_functions_module.o math_tools_module.o random_module.o file_module.o df_module.o thermodynamics_module.o)
$(OBJDIR)init_bc_module.o                    : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o ic_module.o bc_module.o rhs_module.o update_variable_module.o thermodynamics_module.o)
$(OBJDIR)time_module.o                       : $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o shock_detection_module.o output_probes.o)
# output modules (src/output/, via vpath); output_module is the driver
$(OBJDIR)output_restart.o                    : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o profiling_module.o file_module.o)
$(OBJDIR)output_screen.o                     : $(addprefix $(OBJDIR), parameters_module.o profiling_module.o mpi_module.o storage_module.o)
$(OBJDIR)output_stats_channel.o              : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o mesh_module.o onlineStats.o file_module.o fluid_functions_module.o thermodynamics_module.o)
$(OBJDIR)output_stats_tgv.o                  : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o onlineStats.o file_module.o)
$(OBJDIR)output_module.o                     : $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o profiling_module.o inflow_module.o onlineStats.o output_restart.o output_screen.o output_residuals.o output_probes.o output_stats_channel.o output_stats_tgv.o)
