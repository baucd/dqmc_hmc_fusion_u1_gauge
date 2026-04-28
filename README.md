# Build Instructions

This project uses CMake for configuration and compilation. Below is a minimal example and a description of key options.

## Example

```bash
cmake -DCMAKE_Fortran_COMPILER=mpiifx \
      -DENABLE_PIFLUX=ON \
      -DENABLE_DYNERR=ON \
      ../
```

## Options

### CMAKE_Fortran_COMPILER

* Description: Specifies the Fortran compiler used for building the project.
* Example: mpiifx
* Notes:
    * Typically set to an MPI-enabled compiler for parallel execution.
    * Alternatives may include gfortran, ifort, or other MPI wrapper compilers such as mpif90.

### ENABLE_PIFLUX

* Type: Boolean (ON or OFF)
* Description: Enables initialization of the gauge field with a π-flux (pi-flux) configuration.
* Default: OFF (unless otherwise specified in CMakeLists)
* Usage:
    * ON: Use π-flux initial configuration.
    * OFF: Use default or random initial configuration.
* Notes:
    * Useful for simulations where specific topological or symmetry properties are required.
    * May affect convergence behaviour and physical observables.

### ENABLE_DYNERR

* Type: Boolean (ON or OFF)
* Description: Enables dynamic correlation error statistics during runtime.
* Default: OFF
* Usage:
    * ON: Track and update error estimates dynamically.
    * OFF: Disable dynamic error tracking.
* Notes:
    * Can provide more accurate uncertainty estimates.
    * May introduce additional computational overhead.

## General Notes

* All options follow standard CMake conventions: -D<OPTION>=<VALUE>.
* Boolean options accept ON, OFF, TRUE, FALSE.
* It is recommended to build in a separate directory:
```bash
mkdir build && cd build
cmake [options] ..
make -j
```
Additional project-specific options may be available; check CMakeLists.txt for a complete list.
