# Check blacklisted versions first and abort build
if(${BLA_VENDOR} MATCHES "ATLAS")
    message(FATAL_ERROR "The RMG threading model is incompatible with the ATLAS blas threading model. If you really wish to use it comment out the line in PlatformChecks/BlasFilter.inc.")
endif()
if(${BLA_VENDOR} MATCHES "Goto")
    message(FATAL_ERROR "Goto blas is obsolete and not recommended. If you really wish to use it comment out the line in PlatformChecks/BlasFilter.inc.")
endif()

# Check for openblas library supporting changing number of threads
if(${BLAS_LIBRARIES} MATCHES "openblas")
    include("PlatformChecks/CheckOpenBlasSetNumThreads.inc")
    if(RMG_OPENBLAS_SET_NUM_THREADS_INTERNAL)
        add_definitions(-DOPENBLAS_SET_NUM_THREADS=1)
        set(RMG_BLAS_MSG "openblas-parallel.")
    else()
        add_definitions(-DOPENBLAS_SET_NUM_THREADS=0)
        set(RMG_BLAS_MSG "openblas-serial.")
    endif()
endif()

# Check for essl library variants
if(DEFINED BLAS_essl_LIBRARY)
    if(${BLAS_essl_LIBRARY} MATCHES "essl")
        set(BLAS_LIBRARIES ${BLAS_essl_LIBRARY})
        set(BLAS_blas_LIBRARY ${BLAS_essl_LIBRARY})
        if(${BLAS_LIBRARIES} MATCHES "smp" OR ${BLAS_essl_LIBRARY} MATCHES "smp")
            set(RMG_BLAS_MSG "essl-parallel.")
        else()
            set(RMG_BLAS_MSG "essl-serial.")
        endif()
    endif()
endif()

# Mkl library supports changing number of threads
if(RMG_MKLLIBS)
    add_definitions(-DMKLBLAS_SET_NUM_THREADS=1)
    set(RMG_BLAS_MSG "Intel mkl-parallel.")
endif()

# Cray-libsci supports changing number of threads
if(NOT DEFINED RMG_BLAS_MSG)
    if(CRAY_HOST)
        set(RMG_BLAS_MSG "Cray-libsci-parallel.")
    else()
        set(RMG_BLAS_MSG "unknown.")
    endif()
endif()
add_definitions(-DRMG_BLAS_MSG="${RMG_BLAS_MSG}")

