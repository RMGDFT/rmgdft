set(tests "Au-benzene-Au_NCPP")
#list(APPEND tests "C60")

SET(RMG_ON "${CMAKE_BINARY_DIR}/${rmg-on}")
SET(RMG_NEGF "${CMAKE_BINARY_DIR}/${rmg-negf}")
foreach(test ${tests})
    SET(testname NEGF_${test}) 
    SET(TEST_DIR ${test}) 
    SET(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}")

    COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}" )

    RMG_RUN_NEGF_CHECK(${testname} ${TEST_DIR} ${RMG_ON} ${RMG_NEGF} ${num_proc})
endforeach(test ${tests})


