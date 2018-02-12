set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} "-g -DSPGDEBUG -DSPGWARNING")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/spglib/src/version.h version_file)
string(REGEX MATCH "SPGLIB_MAJOR_VERSION ([0-9]+)" spglib_major_version ${version_file})
set(spglib_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "SPGLIB_MINOR_VERSION ([0-9]+)" spglib_minor_version ${version_file})
set(spglib_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "SPGLIB_MICRO_VERSION ([0-9]+)" spglib_micro_version ${version_file})
set(spglib_micro_version ${CMAKE_MATCH_1})
set(serial "${spglib_major_version}.${spglib_minor_version}.${spglib_micro_version}")
set(soserial "1")

# Source code
include_directories("${PROJECT_SOURCE_DIR}/spglib/src")
set(SOURCES ${PROJECT_SOURCE_DIR}/spglib/src/arithmetic.c
            ${PROJECT_SOURCE_DIR}/spglib/src/cell.c
            ${PROJECT_SOURCE_DIR}/spglib/src/debug.c
            ${PROJECT_SOURCE_DIR}/spglib/src/delaunay.c
            ${PROJECT_SOURCE_DIR}/spglib/src/determination.c
            ${PROJECT_SOURCE_DIR}/spglib/src/hall_symbol.c
            ${PROJECT_SOURCE_DIR}/spglib/src/kgrid.c
            ${PROJECT_SOURCE_DIR}/spglib/src/kpoint.c
            ${PROJECT_SOURCE_DIR}/spglib/src/mathfunc.c
            ${PROJECT_SOURCE_DIR}/spglib/src/niggli.c
            ${PROJECT_SOURCE_DIR}/spglib/src/overlap.c
            ${PROJECT_SOURCE_DIR}/spglib/src/pointgroup.c
            ${PROJECT_SOURCE_DIR}/spglib/src/primitive.c
            ${PROJECT_SOURCE_DIR}/spglib/src/refinement.c
            ${PROJECT_SOURCE_DIR}/spglib/src/site_symmetry.c
            ${PROJECT_SOURCE_DIR}/spglib/src/sitesym_database.c
            ${PROJECT_SOURCE_DIR}/spglib/src/spacegroup.c
            ${PROJECT_SOURCE_DIR}/spglib/src/spg_database.c
            ${PROJECT_SOURCE_DIR}/spglib/src/spglib.c
            ${PROJECT_SOURCE_DIR}/spglib/src/spin.c
            ${PROJECT_SOURCE_DIR}/spglib/src/symmetry.c)

# Shared library
add_library(symspg SHARED ${SOURCES})
set_property(TARGET symspg PROPERTY VERSION ${serial})
set_property(TARGET symspg PROPERTY SOVERSION ${soserial})
install(TARGETS symspg LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Static link library
add_library(symspg_static STATIC ${SOURCES})
set_property(TARGET symspg_static PROPERTY VERSION ${serial})
set_property(TARGET symspg_static PROPERTY SOVERSION ${soserial})
set_property(TARGET symspg_static PROPERTY OUTPUT_NAME symspg)
install(TARGETS symspg_static ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# make check
#enable_testing()
#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
#set(CMAKE_CTEST_COMMAND ctest -V)
#add_executable(spglibtest EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/src/test.c ${SOURCES})
#add_test(spglibtest spglibtest)
#add_dependencies(check spglibtest)
