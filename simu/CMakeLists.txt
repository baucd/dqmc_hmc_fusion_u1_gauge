cmake_minimum_required(VERSION 3.10)
enable_language(Fortran)

# hmc
if(ENABLE_HMC)
    add_compile_options(-DHMC)
    add_compile_options(-DCGOPT)
endif()

if(ENABLE_PIFLUX)
    add_compile_options(-DINI_PIFLUX)
endif()
if(ENABLE_DYNERR)
    add_compile_options(-DDYNERROR)
endif()

file(GLOB_RECURSE DIR_SIMU *.f90)
include_directories("${PROJECT_BINARY_DIR}/fqmc/fqmc_math")
include_directories("${PROJECT_BINARY_DIR}/fqmc/fqmc_modl")
include_directories("${PROJECT_BINARY_DIR}/fqmc/fqmc_algo")
add_library(simu STATIC ${DIR_SIMU})

# set exe
add_executable(exe fthmc_main.f90)

# set MKL or other blas/lapack libraries
find_library(mkl_found libmkl_core.a $ENV{MKLROOT}/lib/intel64)
option(DISABLE_MKL "Disable MKL" OFF)
if(mkl_found)
    message(STATUS "MKL libraries found!")
    add_library(mkl_intel_lp64 STATIC IMPORTED)
    add_library(mkl_core STATIC IMPORTED)
    add_library(mkl_sequential STATIC IMPORTED)
    set_property(TARGET mkl_intel_lp64 PROPERTY IMPORTED_LOCATION $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a)
    set_property(TARGET mkl_core PROPERTY IMPORTED_LOCATION $ENV{MKLROOT}/lib/intel64/libmkl_core.a)
    set_property(TARGET mkl_sequential PROPERTY IMPORTED_LOCATION $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a)
    if(ENABLE_OPENMP)
        add_library(mkl_intel_thread STATIC IMPORTED)
        set_property(TARGET mkl_intel_thread PROPERTY IMPORTED_LOCATION $ENV{MKLROOT}/libmkl_intel_thread.a)
        set(mkl_libs mkl_intel_thread mkl_intel_lp64 mkl_core)
    else()
        set(mkl_libs mkl_sequential mkl_intel_lp64 mkl_core)
    endif()
else()
    message(STATUS "MKL libraries not found! Please set the correct MKLROOT or use other blas and lapack library!")
    if(DISABLE_MKL)
        # use other blas lapack library, please add them to dyn lib search path
        #set(mkl_libs blas lapack)
        set(mkl_libs openblas)
    endif()
endif()

# maybe not necessary, set -ldl
if(ENABLE_OPENMP)
    set(omp_lib "dl")
endif()
# set linked library base on different platforms
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    target_link_libraries(exe simu fqmc_algo fqmc_modl fqmc_math -Wl,--start-group ${mkl_libs} -Wl,--end-group -lpthread -lm ${omp_lib})
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(exe simu fqmc_algo fqmc_modl fqmc_math ${mkl_libs} -lpthread -lm ${omp_lib})
else()
    message(FATAL_ERROR "Unknown platform!")
endif()

# install to bin/
install(TARGETS exe DESTINATION ${PROJECT_SOURCE_DIR}/bin)
