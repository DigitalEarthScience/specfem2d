
# REVISION
#   HNG, Jul 12,2011; HNG, Apr 09,2010
# binary directory
BIN = ../bin
include ${PETSC_DIR}/lib/petsc/conf/variables

# serial fortran compiler
#FC = gfortran
#FFLAGS = -g -std=f2003 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wline-truncation -fno-trapping-math -fno-backslash

# parallel fortran compiler
FC = ifort #~/openmpi-1.4.3/bin/mpif90 # OR path to appropriate MPI forrtan compiler
#FFLAGS = -g -std=f2008 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wline-truncation -fno-trapping-math -fno-backslash
FFLAGS = -std08 -xHost -O3 -mcmodel medium -shared-intel -check all -check bounds -fp-stack-check -C -g -traceback -ftz -funroll-loops -unroll5 -align sequence -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage 

#FFLAGS = -std08 -xHost -O3 -g -mcmodel medium -shared-intel -C -ftz -funroll-loops -unroll5 -align sequence -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage 
FC_COMPILE = $(FC) $(FFLAGS)

# C compiler
CC = icc
CFLAGS = -g -O3 -Werror -Wall -pedantic-errors
CC_COMPILE = $(CC) $(CFLAGS)

LIB = 
LIBMPI = -lmpi #-L/~/openmpi-1.4.3/lib/lmpi # OR path to appropriate MPI library
#LIBSCOTCH = -L/home/hgharti/lsoft/scotch_5.1.12/lib -lscotch -lscotcherr
LIBSCOTCH = -L$(HOME)/lsoft/scotch_6.0.4-intel16/lib -lscotch -lscotcherr        
INCSCOTCH = -I$(HOME)/lsoft/scotch_6.0.4-intel16/include

# object tree
OBJ = $(BIN)/global.o           \
      $(BIN)/string_library.o   \
      $(BIN)/math_library.o     \
      $(BIN)/shape_library.o    \
      $(BIN)/read_input.o       \
      $(BIN)/read_mesh.o       \
      $(BIN)/gll_library.o      \
      $(BIN)/mesh_spec.o        \
      $(BIN)/apply_bc.o         \
      $(BIN)/blur_kernel.o      \
      $(BIN)/gindex.o		        \
      $(BIN)/preprocess.o      \
      $(BIN)/visual.o           \
      $(BIN)/postprocess.o      \
	    $(EMPTY_MACRO)

# object tree for SEM serial
OBJ_SEM = $(BIN)/serial_library.o \
          $(BIN)/math_library_serial.o   \
	        $(BIN)/prepare_sparse_serial.o	\
          $(BIN)/solver.o         \
          $(BIN)/solver_petsc.o   \
	        $(EMPTY_MACRO)

# object tree for SEMMPI parallel
OBJ_SEMMPI = $(BIN)/set_precision_mpi.o \
             $(BIN)/mpi_library.o       \
             $(BIN)/math_library_mpi.o  \
             $(BIN)/ghost_library_mpi.o \
	           $(BIN)/prepare_sparse.o	\
             $(BIN)/parsolver.o	\
	    			 $(BIN)/parsolver_petsc.o \
	           $(EMPTY_MACRO)

# object tree for PARTMESH
OBJ_PARTMESH = $(BIN)/global.o           \
               $(BIN)/string_library.o   \
							 $(BIN)/serial_library.o   \
							 $(BIN)/math_library.o     \
               $(BIN)/partmesh_library.o \
							 $(BIN)/read_input.o       \
	             $(BIN)/partmesh_scotch.o  \
               $(BIN)/partmesh.o \
	             $(EMPTY_MACRO)
				
OBJ_semgeotech = $(OBJ) $(OBJ_SEM) $(BIN)/semgeotech.o $(BIN)/maxwell.o $(BIN)/semimage2d.o
OBJ_psemgeotech = $(OBJ) $(OBJ_SEMMPI) $(BIN)/psemgeotech.o $(BIN)/maxwell.o $(BIN)/semimage2d.o

OBJ_semimage = $(OBJ) $(OBJ_SEM) $(BIN)/semimage2d.o $(BIN)/semgeotech.o
OBJ_psemimage = $(OBJ) $(OBJ_SEMMPI) $(BIN)/psemimage2d.o $(BIN)/psemgeotech.o

# main executable for spectral-element method
semgeotech = $(BIN)/semgeotech
psemgeotech = $(BIN)/psemgeotech
semimage = $(BIN)/semimage
psemimage = $(BIN)/psemimage
PARTMESH = $(BIN)/partmesh

default: semimage psemimage
# default: semgeotech psemgeotech

clean:
	rm -f *.mod $(BIN)/*.o core

#all: clean partmesh semimage semgeotech psemgeotech
all: clean partmesh semimage psemimage

semgeotech: $(OBJ_semgeotech)
	$(FC_COMPILE) -o $(semgeotech) $(OBJ_semgeotech) $(LIB)

psemgeotech: $(OBJ_psemgeotech)  
	$(FC_COMPILE) -o $(psemgeotech) $(OBJ_psemgeotech) $(LIBMPI)

semimage: $(OBJ_semimage)
	$(FC_COMPILE) -o $(semimage) $(OBJ_semimage) $(LIB)

psemimage: $(OBJ_psemimage)  
	$(FC_COMPILE) -o $(psemimage) $(OBJ_psemimage) $(LIBMPI) ${PETSC_KSP_LIB}
	
partmesh: $(OBJ_PARTMESH)
	$(FC_COMPILE) -o ${PARTMESH} $(OBJ_PARTMESH) $(LIBSCOTCH)

# compile file which can be compiled as both serial and parllel versions
# serial
$(BIN)/%.serial.o: %.F90
	${FC_COMPILE} -c -o $@ $<
# parallel
$(BIN)/%.parallel.o: %.F90
	${FC_COMPILE} -DUSE_MPI -c -o $@ $<

# compile instruction for each source
$(BIN)/set_precision_mpi.o: set_precision_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/set_precision_mpi.o set_precision_mpi.f90
	
$(BIN)/global.o: global.f90
	$(FC_COMPILE) -c -o $(BIN)/global.o global.f90

$(BIN)/string_library.o: string_library.f90
	$(FC_COMPILE) -c -o $(BIN)/string_library.o string_library.f90

$(BIN)/math_library.o: math_library.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library.o math_library.f90
	
$(BIN)/mpi_library.o: mpi_library.f90
	$(FC_COMPILE) -c -o $(BIN)/mpi_library.o mpi_library.f90
	
$(BIN)/serial_library.o: serial_library.f90
	$(FC_COMPILE) -c -o $(BIN)/serial_library.o serial_library.f90

$(BIN)/math_library_serial.o: math_library_serial.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library_serial.o math_library_serial.f90

$(BIN)/math_library_mpi.o: math_library_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library_mpi.o math_library_mpi.f90

$(BIN)/read_input.o: read_input.f90
	$(FC_COMPILE) -c -o $(BIN)/read_input.o read_input.f90

$(BIN)/read_mesh.o: read_mesh.f90
	$(FC_COMPILE) -c -o $(BIN)/read_mesh.o read_mesh.f90

$(BIN)/gll_library.o: gll_library.f90
	$(FC_COMPILE) -c -o $(BIN)/gll_library.o gll_library.f90

$(BIN)/mesh_spec.o: mesh_spec.f90
	$(FC_COMPILE) -c -o $(BIN)/mesh_spec.o mesh_spec.f90

$(BIN)/shape_library.o: shape_library.f90
	$(FC_COMPILE) -c -o $(BIN)/shape_library.o shape_library.f90

$(BIN)/gauss_library.o: gauss_library.f90
	$(FC_COMPILE) -c -o $(BIN)/gauss_library.o gauss_library.f90
	
$(BIN)/apply_bc.o: apply_bc.f90
	$(FC_COMPILE) -c -o $(BIN)/apply_bc.o apply_bc.f90
	
$(BIN)/blur_kernel.o: blur_kernel.f90
	$(FC_COMPILE) -c -o $(BIN)/blur_kernel.o blur_kernel.f90

$(BIN)/map_location.o: map_location.f90
	$(FC_COMPILE) -c -o $(BIN)/map_location.o map_location.f90

$(BIN)/compute_pressure.o: compute_pressure.f90
	$(FC_COMPILE) -c -o $(BIN)/compute_pressure.o compute_pressure.f90

$(BIN)/prepare_sparse_serial.o: prepare_sparse_serial.f90
	$(FC_COMPILE) -c -o $(BIN)/prepare_sparse_serial.o prepare_sparse_serial.f90

$(BIN)/prepare_sparse.o: prepare_sparse.f90
	$(FC_COMPILE) -c -o $(BIN)/prepare_sparse.o prepare_sparse.f90

$(BIN)/gindex.o: gindex.f90
	$(FC_COMPILE) -c -o $(BIN)/gindex.o gindex.f90

$(BIN)/ghost_library_mpi.o: ghost_library_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/ghost_library_mpi.o ghost_library_mpi.f90

$(BIN)/preprocess.o: preprocess.f90
	$(FC_COMPILE) -c -o $(BIN)/preprocess.o preprocess.f90

$(BIN)/visual.o: visual.f90
	$(FC_COMPILE) -c -o $(BIN)/visual.o visual.f90
	
$(BIN)/postprocess.o: postprocess.f90
	$(FC_COMPILE) -c -o $(BIN)/postprocess.o postprocess.f90
	
$(BIN)/solver.o: solver.f90
	$(FC_COMPILE) -c -o $(BIN)/solver.o solver.f90

$(BIN)/solver_petsc.o: solver_petsc.f90
	$(FC_COMPILE) -c -o $(BIN)/solver_petsc.o solver_petsc.f90
	
$(BIN)/parsolver.o: parsolver.f90
	$(FC_COMPILE) -c -o $(BIN)/parsolver.o parsolver.f90	

$(BIN)/parsolver_petsc.o: parsolver_petsc.f90
	$(FC_COMPILE) ${PETSC_FC_INCLUDES} -cpp -c -o $(BIN)/parsolver_petsc.o parsolver_petsc.f90 ${PETSC_KSP_LIB}

# following objective files are compiled from the same file but with 
# different compiler preprocessing directives options.
# -DUSE_MPI provides the flag to use the MPI library for the compilation.
# the flag -DUSE_MPI is equivalent to #define USE_MPI 1.
$(BIN)/semimage2d.o: semimage2d.F90
	$(FC_COMPILE) -c -o $(BIN)/semimage2d.o semimage2d.F90
	
$(BIN)/psemimage2d.o: semimage2d.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/psemimage2d.o semimage2d.F90
	
# following objective files are compiled from the same file but with 
# different compiler preprocessing directives options.
$(BIN)/semgeotech.o: semgeotech.F90
	$(FC_COMPILE) -c -o $(BIN)/semgeotech.o semgeotech.F90
	
$(BIN)/psemgeotech.o: semgeotech.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/psemgeotech.o semgeotech.F90
	
# mesh partition
$(BIN)/partmesh_library.o: partmesh_library.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh_library.o partmesh_library.f90

$(BIN)/partmesh_scotch.o: partmesh_scotch.f90 partmesh_library.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh_scotch.o partmesh_scotch.f90 $(INCSCOTCH)

$(BIN)/partmesh.o: partmesh.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh.o partmesh.f90

