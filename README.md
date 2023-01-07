# HDF5_utils_mpi

This library is to provide a high level interface into [HDF5](https://portal.hdfgroup.org/display/support).
This is, of course, based off HDF5's own High-Level module ([HDLT](https://support.hdfgroup.org/HDF5/doc/HL/RM_H5LT.html)) but customized to my own needs.

The aim of this library is to:
 - **MPI Parallel data reading/writing with HDF5** (new feature in this fork)
 - Abstract most of the HDF5 tyes.
 - Read/write datasets of multiple types and ranks.
   - Assumed reading/writing full arrays.
   - Remove need to pass dimensions (allocation and consistency of dimensions is left to the user).
 - Write/read attributes to file/group/datasets.
 - Provide the ability to create groups, and access datasets by either absolute or relative path.
 - Some auxiliary functions:
   - Check if dataset/object exists.
   - Get rank and size of dataset, either to check dimensions or to allocate array before reading.

## Module summary

The doxygen generated document can be accessed at [https://luohancfd.github.io/HDF5_utils_mpi/namespacehdf5__utils__mpi.html](https://luohancfd.github.io/HDF5_utils_mpi/namespacehdf5__utils__mpi.html)

## Examples

Check the `example` folder

## Copyright

Copyright (c) 2017-2019 Justin Erwin

Copyright (c) 2020      Han Luo
