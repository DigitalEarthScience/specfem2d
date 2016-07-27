# Makefile for SPECFEM3D_MAXWELL
# REVISION
#   HNG, Jan 17,2014

default: createdir semimage

all: createdir default partmesh semimage psemimage

createdir:
	(mkdir -p bin; mkdir -p input; mkdir -p output; mkdir -p partition; mkdir -p tmp)

clean:
	(cd src; make clean)

semimage: 
	(cd src; make $@)

psemimage: 
	(cd src; make $@)

partmesh: 
	(cd src; make $@)

