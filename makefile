
include $(PETSC_DIR)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/MPICH

OBJS = simpleFEM main.o generateesm.o readmesh.o
LD_LIBRARY_PATH= $(PETSC_DIR)/$(PETSC_ARCH)/lib
FFLAGS = -x f95-cpp-input -I$(PETSC_DIR)/include/petsc/finclude/ -I$(PETSC_DIR)/include/ -I$(PETSC_DIR)/$(PETSC_ARCH)/include/ -Dmypetscpath=$(PETSC_DIR)


all: main.o generateesm.o readmesh.o
	$(FC) -L$(LD_LIBRARY_PATH) -o simpleFEM main.o generateesm.o readmesh.o -lpetsc

main.o: main.f95
	$(FC) $(FFLAGS) -c main.f95

generateesm.o: generateesm.f95
	$(FC) $(FFLAGS) -c generateesm.f95

readmesh.o: readmesh.f95
	$(FC) $(FFLAGS) -c readmesh.f95

clean:
	rm -rf $(OBJS)


