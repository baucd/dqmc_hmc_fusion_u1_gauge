cmake_minimum_required(VERSION 3.18)

# ===== 编译选项（对应 make.sys）=====
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")

include_directories("${PROJECT_BINARY_DIR}/fqmc/fqmc_math")
include_directories("${PROJECT_BINARY_DIR}/fqmc/fqmc_modl")
include_directories("${PROJECT_BINARY_DIR}/fqmc/fqmc_algo")
include_directories("${PROJECT_BINARY_DIR}/simu")

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

# ===== 分析程序列表 =====
set(PROGS
    rebin
    jackv5
    jackv5_sq
    jackv5_sq_ori
    jackv5_sq_0to2pi
    jackv5_ising
    jackv5_coratio
    jackv5_coratio_y
    jackv5_coratio_spin
    jackv5_coratio_gamma
    jackv5_gtau
    jackv5_dynsus
    jackv5_dynsus_swave
)


# ===== 批量生成可执行文件 =====
foreach(p ${PROGS})
    add_executable(${p} ${p}.f90)

    target_link_libraries(${p}
        PRIVATE
        fqmc_math
        -Wl,--start-group ${mkl_libs} -Wl,--end-group -lpthread -lm
    )

    install(TARGETS ${p} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach()

# ===== 特殊程序 =====

# binder（多源文件）
add_executable(binder
    errcalc_binder.f90
    binder.f90
)

target_link_libraries(binder
    PRIVATE
    fqmc_math
    -Wl,--start-group ${mkl_libs} -Wl,--end-group -lpthread -lm
)

# ===== 特殊程序 part2 =====
set(PROGS
    inv_fourier
    inv_fourier_dimer
)


# ===== 批量生成可执行文件 =====
foreach(p ${PROGS})
    add_executable(${p} ${p}.f90)

    target_link_libraries(${p}
        PRIVATE
        simu
        fqmc_modl
        fqmc_algo
        fqmc_math
        -Wl,--start-group ${mkl_libs} -Wl,--end-group -lpthread -lm
    )

    install(TARGETS ${p} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach()

install(TARGETS binder DESTINATION ${PROJECT_SOURCE_DIR}/bin)
