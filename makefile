
include $(PETSC_DIR)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/MPICH

OBJS = simpleFEM main.o generateesm.o readmesh.o readBC.o recoverstress.o
LD_LIBRARY_PATH= $(PETSC_DIR)/$(PETSC_ARCH)/lib
FFLAGS = -cpp -I$(PETSC_DIR)/include/petsc/finclude/ -I$(PETSC_DIR)/include/ -I$(PETSC_DIR)/$(PETSC_ARCH)/include/ -Dmypetscpath=$(PETSC_DIR) 
FC= mpifort
all: main.o generateesm.o readmesh.o readBC.o recoverstress.o
	$(FC) -Wl,-rpath=$(LD_LIBRARY_PATH) -L$(LD_LIBRARY_PATH) -o simpleFEM main.o generateesm.o readmesh.o readBC.o recoverstress.o -lpetsc

main.o: main.f08
	$(FC) $(FFLAGS) -c main.f08

generateesm.o: generateesm.f08
	$(FC) $(FFLAGS) -c generateesm.f08

readmesh.o: readmesh.f08
	$(FC) $(FFLAGS) -c readmesh.f08

readBC.o: readBC.f08
	$(FC) $(FFLAGS) -c readBC.f08

recoverstress.o: recoverstress.f08
	$(FC) $(FFLAGS) -c recoverstress.f08

clean:
	rm -rf $(OBJS)


