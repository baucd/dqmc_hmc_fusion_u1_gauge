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

### ENABLE_HMC
* Description: Use HMC code.

### ENABLE_CUDA_CG
* Description: Use CUDA and nvidia GPU computation for HMC CG. Not CUDA, because DQMC code also use CUDA flag but not very usable.

### ENABLE_CUDA_MEAS
* Description: Use CUDA and nvidia GPU computation for HMC measurements.
* Usage: Must with **ENABLE_CUDA_CG** before **ENABLE_CUDA_MEAS**.

### ENABLE_CGOPT
* Description: Use optimized CG instead of naive one.

### ENABLE_PCG
* Description: Use preconditioned CG.

### ENABLE_DIRINV
* Description: Use inversion of matrix to get Green's function.

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
