

FC = mpiifort
# GIT tag usage
# >>> git tag -a v1.4 -m "my version 1.4"
HDF5_FFLAG = -I${HDF5INCLUDE}
HDF5_LDFLAG = -L${HDF5DIR} -lhdf5_fortran

# add hdf5 support
LDFLAGS = $(HDF5_LDFLAG) #`pkg-config hdf5 --libs` -lhdf5_fortran
FFLAGS = $(HDF5_FFLAG) #`pkg-config hdf5 --cflags`

ifeq ($(FC),gfortran)

  #FFLAGS = -O2
  FFLAGS += -O0 -Wall -fcheck=all
  #FFLAGS = -O0 -Wall -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,underflow
  #FFLAGS = -O0 -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
  #FFLAGS = -O3 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal

  FFLAGS += -cpp

endif
ifeq ($(FC),mpiifort)

  FFLAGS += -O0
  #FFLAGS = -Ofast

  FFLAGS += -fpp

endif

FFLAGS += -g# -Wall


.PHONY: clean docs

all: hdf5_test_1.x hdf5_test_2.x

hdf5_test_1.x: parallel.o hdf5_utils_mpi.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

hdf5_test_2.x: parallel_change_core.o hdf5_utils_mpi.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

parallel_change_core.o: hdf5_utils_mpi.o

parallel.o: hdf5_utils_mpi.o

%.o : %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# hdf5_utils_mpi.f90 : FORCE
# 	cd core && python gen_code.py

# FORCE:

docs:
	cd docs; doxygen HDF5_utils.doxy

clean:
	rm -f *.o *.mod *.x
