

FC = mpif90.openmpi

# GIT tag usage
# >>> git tag -a v1.4 -m "my version 1.4"
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

# add hdf5 support
LDFLAGS = `pkg-config hdf5 --libs` -lhdf5_fortran
FFLAGS = `pkg-config hdf5 --cflags`

ifeq ($(FC),gfortran)

  #FFLAGS = -O2
  FFLAGS += -O0 -Wall -fcheck=all
  #FFLAGS = -O0 -Wall -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,underflow
  #FFLAGS = -O0 -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
  #FFLAGS = -O3 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal

  FFLAGS += -cpp -DVERSION=\"$(GIT_VERSION)\"

endif
ifeq ($(FC),ifort)

  FFLAGS += -O0
  #FFLAGS = -Ofast

  FFLAGS += -fpp -DVERSION=\"$(GIT_VERSION)\"

endif

FFLAGS += -g

$(info $(FFLAGS))



.PHONY: clean docs


hdf5_test.x: parallel.o hdf5_utils_mpi.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

parallel.o: hdf5_utils_mpi.o

%.o : %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

docs:
	cd docs; doxygen HDF5_utils.doxy

clean:
	rm -f *.o *.mod *.x
