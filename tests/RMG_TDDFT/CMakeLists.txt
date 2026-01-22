set(tests "Si-8atom-gamma-double-scalapack")
list(APPEND tests "Si-8atom-gamma-double-tiledMM")
list(APPEND tests "Si-8atom-2x2x2-double-tiledMM")
list(APPEND tests "Si-8atom-2x2x2-double-scalapack")

SET(RMG_EXE "${CMAKE_BINARY_DIR}/${rmg}")
foreach(test ${tests})
    SET(testname TDDFT_${test}) 
    SET(TEST_DIR ${test}) 
    SET(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}")

    COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}" )

    RMG_RUN_TDDFT_CHECK(${testname} ${TEST_DIR} ${RMG_EXE} ${num_proc})
endforeach(test ${tests})


