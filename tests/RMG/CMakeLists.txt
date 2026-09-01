set(energy_tests "atomO_polarized")
list(APPEND energy_tests "Diamond2")
list(APPEND energy_tests "Diamond16")
list(APPEND energy_tests "C60")
list(APPEND energy_tests "CO_ccECP")
list(APPEND energy_tests "Pt-bulk-spinorbit")
list(APPEND energy_tests "Mg_2atom")
list(APPEND energy_tests "Re_2atom")
list(APPEND energy_tests "Fe_2atom")
list(APPEND energy_tests "AlN32")
list(APPEND energy_tests "Si-8atoms_EXX_kpoints")
list(APPEND energy_tests "NiO8")
list(APPEND energy_tests "graphite_stress")
list(APPEND energy_tests "Ni_Semilocal")
list(APPEND energy_tests "graphene")

set (bond_tests "CO")

SET (CTEST_ENVIRONMENT "OMP_NUM_THREADS=1" "RMG_NUM_THREADS=1")
SET(RMG_EXE "${CMAKE_BINARY_DIR}/${rmg}")
foreach(test ${energy_tests})
    SET(testname RMG_${test}) 
    SET(TEST_DIR ${test}) 
    SET(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}")

    COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}" )

    RMG_RUN_ENERGY_CHECK(${testname} ${TEST_DIR} ${RMG_EXE} ${num_proc})
endforeach(test ${tests})


foreach(test ${bond_tests})
    SET(testname RMG_${test}) 
    SET(TEST_DIR ${test}) 
    SET(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}")

    COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}" )

    RMG_RUN_BOND_CHECK(${testname} ${TEST_DIR} ${RMG_EXE} ${num_proc})
endforeach(test ${tests})


