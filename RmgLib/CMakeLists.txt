

option (SEPARATE_BUILD OFF)

if (SEPARATE_BUILD)
cmake_minimum_required (VERSION 2.6)
project(RmgLib)
set (PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include})
include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
    if (SEPARATE_BUILD STREQUAL "ON")
        find_package(Boost REQUIRED COMPONENTS thread system)
        include_directories(${Boost_INCLUDE_DIRS})
        # Search for an MPI implementation
        find_package(MPI REQUIRED)
        if(MPI_CXX_FOUND)
            include_directories(${MPI_CXX_INCLUDE_PATH})
        else (MPI_CXX_FOUND)
            include(CMakeForceCompiler)
            CMAKE_FORCE_CXX_COMPILER(CC ${CMAKE_CXX_COMPILER_ID})
        endif (MPI_CXX_FOUND)
        if(MPI_C_FOUND)
            include_directories(${MPI_C_INCLUDE_PATH})
        else(MPI_C_FOUND)
          CMAKE_FORCE_C_COMPILER(cc ${CMAKE_C_COMPILER_ID})
        endif (MPI_C_FOUND)

        # Compiler flags
        if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "-Wall -O3 -std=c++11")        ## Optimize
        endif()

    endif()
endif (SEPARATE_BUILD)

set (RmgLibSources 
src/TradeImages.cpp
src/Lattice.cpp
src/BaseThread.cpp
src/BaseGrid.cpp
src/Mgrid.cpp
src/MpiQueue.cpp
src/app_cir_driver.cpp
src/app_cil_driver.cpp
src/app_cil_driver_threaded.cpp
src/app_grad_driver.cpp
src/app_del2_driver.cpp
src/FiniteDiff.cpp
src/RmgTimer.cpp
src/RmgPrintTimings.cpp
src/RmgSumAll.cpp
src/RmgMaxAll.cpp
src/ThreadBarrier.cpp
src/rmg_error.cpp
src/get_vh.cpp
src/packfuncs.cpp
)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_library (RmgLibShared SHARED ${RmgLibSources})
endif()
add_library (RmgLib STATIC ${RmgLibSources})
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

    add_executable (poisson_pbc examples/Poisson/poisson_pbc.cpp)
    target_link_libraries (poisson_pbc ${Boost_LIBRARIES} RmgLibShared)
    if(DEFINED ENV{CRAY_LIBSCI_VERSION})
        target_link_libraries (poisson_pbc ${Boost_LIBRARIES} RmgLibShared)
    else(DEFINED ENV{CRAY_LIBSCI_VERSION})
        target_link_libraries (poisson_pbc ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} RmgLibShared)
    endif()

endif()


