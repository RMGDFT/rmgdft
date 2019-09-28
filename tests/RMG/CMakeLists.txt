set(tests "atomO_polarized")
list(APPEND tests "C60")

SET(RMG_EXE "${CMAKE_BINARY_DIR}/rmg-cpu")
foreach(test ${tests})
    SET(testname RMG_${test}) 
    SET(TEST_DIR ${test}) 
    SET(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}")

    COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}" )

    RMG_RUN_CHECK(${testname} ${TEST_DIR} ${RMG_EXE} ${num_proc})
endforeach(test ${tests})

