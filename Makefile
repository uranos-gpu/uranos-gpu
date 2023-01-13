#Programe Name
MAINPROG = Uranos.exe
MAINPOST = PostUranos.exe
MAINSTAT = PostUranosStat.exe

# Directories
DIRECTORIES = .objs .mods DATA
OBJDIR = .objs/
MODDIR = .mods/
SRCDIR = src/
PSTDIR = pst_src/
STADIR = pst_src_stat/
LIBDIR = libs/

# Compiler
FC = mpif90

#LIBS = -L/usr/local/opt/lapack/lib -llapack -lblas

# Libraries
ifeq ($(HOME), /marconi/home/userexternal/fdevanna)
	LIBS = #-L/cineca/prod/opt/compilers/intel/pe-xe-2018/binary/lib/intel64/ -L${LAPACK_LIB} -L${BLAS_LIB} -lifcoremt -lblas  -llapack
else
	LIBS = #-Llibs/lapack/lib -llapack -lblas
endif

$cudaver=$(nvcc --version | sed -n -e 's/^.*release //p' | awk '{ print substr($1, 1, length($1)-1) }')

# DEFAULTS!!!
# -----------------------------
OPT          = #-O3 -mtune=native -mcmodel=medium -fdefault-real-8
DEBUG        = #-g -DDEBUG -pedantic -fbounds-check -fcheck=all -Wall -Waliasing -Wextra -fmax-errors=5 #-ffpe-trap=zero,overflow,underflow #-Werror
OPENMP       = #-fopenmp -DOPENMP -lgomp -lpthread
TIME 	     = 
MODULES      = #-J$(MODDIR)
# SMART COMPILING OPTIONS
# -----------------------------


ifeq ($(comp),gnu)
     FC = mpif90
     MODULES = -J$(MODDIR)
     OPT     = -O3
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
ifeq ($(comp),gnuch)
     FC = mpif90.mpich
     MODULES = -J$(MODDIR)
     OPT     = -O3
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



ifeq ($(comp),pgi)
	FC = mpif90
	MODULES = -module $(MODDIR)
	ifeq ($(mode),debug)
		OPT    = -O3
		DEBUG  = 
		TIME   = -DNVTX -lnvToolsExt
		OPENMP =
	else ifeq ($(mode),time)
		OPT    = -O3
		DEBUG  =
		TIME   = -DTIME
		OPENMP =
	else ifeq ($(mode),gpu)
		OPT    = -O3 -g -acc -ta=tesla:deepcopy -cudalib=curand#,lineinfo
		DEBUG  = 
		TIME   =
		OPENMP =
	else ifeq ($(mode),debug_gpu)
		OPT    = -O3 -acc -ta=tesla -cudalib=curand
		DEBUG  = -Minfo=accel -g
		TIME   = -DNVTX -lnvToolsExt
		OPENMP =
	else ifeq ($(mode),debug_gpu_marconi100)
		FC     = mpipgifort
		OPT    = -O3 -acc -ta=tesla:deepcopy
		DEBUG  = -Minfo=accel -g
		TIME   = -DNVTX -cudalib=curand,nvtx3
		OPENMP =
	else ifeq ($(mode),marconi100)
		FC     = mpipgifort
		OPT    = -O3 -acc -ta=tesla:deepcopy -cudalib=curand
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




ifeq ($(fast_pgi),true)
	FC      = mpif90
	OPT     = -O3 -r8 -mcmodel=medium
	OPENMP  = 
	DEBUG   =
	TIME    = 
	MODULES = -module $(MODDIR)
endif

ifeq ($(pgi-cuda),true)
	FC      = mpif90
	OPT     = -fast -r8 #-acc
	TIME    = 
	OPENMP	= 
	DEBUG   = -Minfo=accel
	MODULES = -module $(MODDIR)
endif

ifeq ($(mpidex),buffer)
	MPIDEX = -DMPIBFR
endif

LDFLAGS := -cpp $(OPENMP) $(TIME) $(MPIDEX) $(DEBUG) $(OPT)

# Fortran Compiler flags
FFLAGS    = $(LDFLAGS)

#Object files
OBJECTS = $(addprefix $(OBJDIR), mpi_module.o \
	  parameters_module.o \
	  mesh_module.o \
	  input_module.o \
	  output_module.o \
	  output_tch.o \
	  storage_module.o \
	  allocate_module.o \
	  time_module.o \
	  ic_module.o \
	  init_bc_module.o \
	  bc_module.o \
	  rhs_module.o \
	  advection_module.o \
	  shock_detection_module.o \
	  eigen_matrix_module.o \
	  fluid_functions_module.o \
	  viscous_module.o \
	  flux_module.o \
	  mpi_comm_module.o \
	  inflow_module.o \
	  sgs_module.o \
	  df_module.o \
	  GetRetau_module.o \
	  statistics_module.o \
	  post_storage_module.o \
	  post_input_module.o \
	  post_shell_tools_module.o \
	  post_computation_module.o \
	  post_output_gnu_module.o \
	  post_output_vtk_module.o \
	  post_output_module.o \
	  post_bc_module.o \
	  post_solutions_module.o \
	  post_statistic_module.o \
	  post_statistic_spatial_module.o \
	  post_statistic_time_module.o \
	  post_stat_input_module.o \
	  post_stat_storage_module.o \
	  post_stat_output_module.o \
	  post_stat_computation_module.o \
	  post_stat_output_bl_module.o \
	  norm_module.o \
	  integration_module.o \
	  real_to_integer_module.o \
	  interpolation_module.o \
	  dirac_delta_module.o \
	  math_tools_module.o \
	  matrix_inversion_module.o \
	  performance_module.o \
	  random_module.o \
	  reynolds_averaged_module.o \
	  file_module.o \
	  wmles_module.o \
	  onlineStats.o \
	  vtk_utils_module.o\
	  nvtx.o \
	  npy.o \
	  profiling_module.o)



# make all 
all : $(DIRECTORIES) $(MAINPROG) $(MAINPOST) $(MAINSTAT) $(DEBUGGER) $(GRSCALER)

# Create the directory $@ if this does not exist
$(DIRECTORIES) :
	mkdir -p $@
	
# -------------------------------------------------------
# Compile Uranos main program 
# -------------------------------------------------------
$(MAINPROG): $(OBJECTS) $(SRCDIR)main.f90 Makefile
	@echo "-------------------------------------------"
	@echo "Creating the executable: $(MAINPROG)"
	@echo "-------------------------------------------"
	$(FC) $(LDFLAGS) $(SRCDIR)main.f90 $(MODULES) $(OBJECTS) -o $(MAINPROG) $(LIBS)

# -------------------------------------------------------
# Compile post_Uranos main program 
# -------------------------------------------------------
$(MAINPOST): $(OBJECTS) $(PSTDIR)post_main.f90 Makefile
	@echo "-------------------------------------------"
	@echo "Creating the executable: $(MAINPOST)"
	@echo "-------------------------------------------"
	$(FC) $(LDFLAGS) $(PSTDIR)post_main.f90 $(MODULES) $(OBJECTS) -o $(MAINPOST) $(LIBS)

# -------------------------------------------------------
# Compile post_Uranos_stat main program 
# -------------------------------------------------------
$(MAINSTAT): $(OBJECTS) $(STADIR)post_stat_main.f90 Makefile
	@echo "-------------------------------------------"
	@echo "Creating the executable: $(MAINSTAT)"
	@echo "-------------------------------------------"
	$(FC) $(LDFLAGS) $(STADIR)post_stat_main.f90 $(MODULES) $(OBJECTS) -o $(MAINSTAT) $(LIBS)


# -------------------------------------------------------
# Compiling modules in SRC
# -------------------------------------------------------
$(OBJDIR)%.o : $(SRCDIR)%.f90 Makefile 
	$(FC) $(FFLAGS) -c $< -o $@ $(MODULES) $(LIBS)

# -------------------------------------------------------
# Compiling modules in PST_SRC
# -------------------------------------------------------
$(OBJDIR)%.o : $(PSTDIR)%.f90 Makefile
	$(FC) $(FFLAGS) -c $< -o $@ $(MODULES) $(LIBS)

# -------------------------------------------------------
# Compiling modules in PST_SRC_STAT
# -------------------------------------------------------
$(OBJDIR)%.o : $(STADIR)%.f90 Makefile
	$(FC) $(FFLAGS) -c $< -o $@ $(MODULES) $(LIBS)

# -------------------------------------------------------
# Compiling libraries
# -------------------------------------------------------
$(OBJDIR)%.o : $(LIBDIR)%.f90 Makefile
	$(FC) $(FFLAGS) -c $< -o $@ $(MODULES) $(LIBS)

	
# Clean everything up
.PHONY: clean
clean: 
	@echo "---------------------------------------"
	@echo "Cleaning everything up"
	@echo "---------------------------------------"
	rm -rf *.exe *.o *.mod rm -rf *dSYM
	rm -rf *.optrpt
	rm -rf $(OBJDIR)
	rm -rf $(MODDIR)

#---------------------------------------------------------------------------
# Dependency chains
#---------------------------------------------------------------------------
$(OBJDIR)parameters_module.o			: $(SRCDIR)parameters_module.f90
$(OBJDIR)allocate_module.o			: $(SRCDIR)allocate_module.f90 $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)mpi_comm_module.o			: $(SRCDIR)mpi_comm_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o)
$(OBJDIR)mpi_module.o				: $(SRCDIR)mpi_module.f90 $(addprefix $(OBJDIR), parameters_module.o performance_module.o file_module.o)
$(OBJDIR)input_module.o				: $(SRCDIR)input_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o real_to_integer_module.o)
$(OBJDIR)mesh_module.o				: $(SRCDIR)mesh_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o math_tools_module.o input_module.o matrix_inversion_module.o)
$(OBJDIR)output_module.o			: $(SRCDIR)output_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o norm_module.o statistics_module.o fluid_functions_module.o file_module.o integration_module.o wmles_module.o fluid_functions_module.o onlineStats.o output_tch.o)
$(OBJDIR)output_tch.o				: $(SRCDIR)output_tch.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o storage_module.o onlineStats.o file_module.o)
$(OBJDIR)storage_module.o			: $(SRCDIR)storage_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o mesh_module.o allocate_module.o)
$(OBJDIR)time_module.o				: $(SRCDIR)time_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o shock_detection_module.o)
$(OBJDIR)ic_module.o				: $(SRCDIR)ic_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o fluid_functions_module.o math_tools_module.o random_module.o file_module.o df_module.o)
$(OBJDIR)init_bc_module.o		 	: $(SRCDIR)init_bc_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o ic_module.o bc_module.o rhs_module.o)
$(OBJDIR)bc_module.o			 	: $(SRCDIR)bc_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o fluid_functions_module.o inflow_module.o real_to_integer_module.o mpi_comm_module.o profiling_module.o)
$(OBJDIR)rhs_module.o			 	: $(SRCDIR)rhs_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o advection_module.o viscous_module.o integration_module.o sgs_module.o profiling_module.o)
$(OBJDIR)advection_module.o			: $(SRCDIR)advection_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o eigen_matrix_module.o flux_module.o inflow_module.o profiling_module.o)
$(OBJDIR)eigen_matrix_module.o			: $(SRCDIR)eigen_matrix_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o)
$(OBJDIR)shock_detection_module.o		: $(SRCDIR)shock_detection_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o fluid_functions_module.o mpi_comm_module.o)
$(OBJDIR)viscous_module.o			: $(SRCDIR)viscous_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o nvtx.o)
$(OBJDIR)flux_module.o				: $(SRCDIR)flux_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mpi_module.o)
$(OBJDIR)inflow_module.o			: $(SRCDIR)inflow_module.f90 $(addprefix $(OBJDIR), storage_module.o fluid_functions_module.o math_tools_module.o df_module.o)
$(OBJDIR)sgs_module.o				: $(SRCDIR)sgs_module.f90 $(addprefix $(OBJDIR), storage_module.o mpi_module.o parameters_module.o bc_module.o wmles_module.o)
$(OBJDIR)df_module.o				: $(SRCDIR)df_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o fluid_functions_module.o)
$(OBJDIR)performance_module.o			: $(SRCDIR)performance_module.f90 $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)GetRetau_module.o			: $(SRCDIR)GetRetau_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o fluid_functions_module.o mesh_module.o)
$(OBJDIR)wmles_module.o				: $(SRCDIR)wmles_module.f90 $(addprefix $(OBJDIR), parameters_module.o fluid_functions_module.o)
$(OBJDIR)profiling_module.o	                : $(SRCDIR)profiling_module.f90 $(addprefix $(OBJDIR), nvtx.o parameters_module.o)
$(OBJDIR)nvtx.o					: $(SRCDIR)nvtx.f90


$(OBJDIR)post_storage_module.o			: $(PSTDIR)post_storage_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o mesh_module.o)
$(OBJDIR)post_input_module.o			: $(PSTDIR)post_input_module.f90 $(addprefix $(OBJDIR), post_storage_module.o post_shell_tools_module.o)
$(OBJDIR)post_shell_tools_module.o		: $(PSTDIR)post_shell_tools_module.f90 $(addprefix $(OBJDIR), post_storage_module.o)
$(OBJDIR)post_computation_module.o		: $(PSTDIR)post_computation_module.f90 $(addprefix $(OBJDIR), post_storage_module.o storage_module.o fluid_functions_module.o shock_detection_module.o post_solutions_module.o)
$(OBJDIR)post_output_gnu_module.o		: $(PSTDIR)post_output_gnu_module.f90 $(addprefix $(OBJDIR), post_storage_module.o storage_module.o norm_module.o post_solutions_module.o post_statistic_module.o post_computation_module.o file_module.o fluid_functions_module.o statistics_module.o)  
$(OBJDIR)post_output_vtk_module.o		: $(PSTDIR)post_output_vtk_module.f90 $(addprefix $(OBJDIR), post_storage_module.o storage_module.o post_computation_module.o vtk_utils_module.o)
$(OBJDIR)post_output_module.o			: $(PSTDIR)post_output_module.f90 $(addprefix $(OBJDIR), post_storage_module.o post_output_gnu_module.o storage_module.o post_output_vtk_module.o npy.o file_module.o)
$(OBJDIR)post_bc_module.o			: $(PSTDIR)post_bc_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o post_storage_module.o bc_module.o)
$(OBJDIR)post_solutions_module.o		: $(PSTDIR)post_solutions_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o post_storage_module.o)
$(OBJDIR)post_statistic_module.o		: $(PSTDIR)post_statistic_module.f90 $(addprefix $(OBJDIR), post_storage_module.o storage_module.o post_statistic_time_module.o post_statistic_spatial_module.o)
$(OBJDIR)post_statistic_spatial_module.o	: $(PSTDIR)post_statistic_spatial_module.f90 $(addprefix $(OBJDIR), post_storage_module.o )
$(OBJDIR)post_statistic_time_module.o		: $(PSTDIR)post_statistic_time_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o)

# POST STAT MODULES
$(OBJDIR)post_stat_storage_module.o		: $(STADIR)post_stat_storage_module.f90 $(addprefix $(OBJDIR), parameters_module.o file_module.o allocate_module.o mpi_module.o storage_module.o)
$(OBJDIR)post_stat_input_module.o		: $(STADIR)post_stat_input_module.f90 $(addprefix $(OBJDIR), post_shell_tools_module.o post_stat_storage_module.o)
$(OBJDIR)post_stat_output_module.o		: $(STADIR)post_stat_output_module.f90 $(addprefix $(OBJDIR), parameters_module.o post_stat_storage_module.o mpi_module.o post_stat_output_bl_module.o)
$(OBJDIR)post_stat_computation_module.o		: $(STADIR)post_stat_computation_module.f90 $(addprefix $(OBJDIR), parameters_module.o post_stat_storage_module.o)
$(OBJDIR)post_stat_output_bl_module.o		: $(STADIR)post_stat_output_bl_module.f90 $(addprefix $(OBJDIR), parameters_module.o post_stat_storage_module.o wmles_module.o)



#---------------------------------------------------------------------------
# Libraries
#---------------------------------------------------------------------------
$(OBJDIR)norm_module.o				: $(LIBDIR)norm_module.f90 $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)fluid_functions_module.o		: $(LIBDIR)fluid_functions_module.f90 $(addprefix $(OBJDIR), parameters_module.o storage_module.o interpolation_module.o)
$(OBJDIR)integration_module.o			: $(LIBDIR)integration_module.f90 $(addprefix $(OBJDIR), parameters_module.o mpi_module.o)
$(OBJDIR)real_to_integer_module.o		: $(LIBDIR)real_to_integer_module.f90 $(addprefix $(OBJDIR), parameters_module.o performance_module.o)
$(OBJDIR)dirac_delta_module.o			: $(LIBDIR)dirac_delta_module.f90 $(addprefix $(OBJDIR), parameters_module.o real_to_integer_module.o)
$(OBJDIR)math_tools_module.o			: $(LIBDIR)math_tools_module.f90 $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)matrix_inversion_module.o		: $(LIBDIR)matrix_inversion_module.f90 $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)interpolation_module.o			: $(LIBDIR)interpolation_module.f90 $(addprefix $(OBJDIR), parameters_module.o real_to_integer_module.o matrix_inversion_module.o mpi_module.o)
$(OBJDIR)random_module.o			: $(LIBDIR)random_module.f90 $(addprefix $(OBJDIR), parameters_module.o) 
$(OBJDIR)reynolds_averaged_module.o		: $(LIBDIR)reynolds_averaged_module.f90 
$(OBJDIR)file_module.o				: $(LIBDIR)file_module.f90 $(addprefix $(OBJDIR), parameters_module.o)  
$(OBJDIR)vtk_utils_module.o			: $(LIBDIR)vtk_utils_module.f90 $(addprefix $(OBJDIR), parameters_module.o file_module.o)
$(OBJDIR)statistics_module.o			: $(LIBDIR)statistics_module.f90 $(addprefix $(OBJDIR), parameters_module.o)
$(OBJDIR)onlineStats.o	   			: $(LIBDIR)onlineStats.f90 $(addprefix $(OBJDIR), parameters_module.o statistics_module.o fluid_functions_module.o)
$(OBJDIR)npy.o	                     		: $(LIBDIR)npy.f90 








