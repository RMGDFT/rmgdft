MESSAGE(STATUS "Creating tests for RMG")
if (NOT num_proc)
SET(num_proc 16)
endif()

SET(testname "atomO_polarized") 
SET(SRC_DIR "${CMAKE_SOURCE_DIR}/tests/${testname}")

COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${testname}" )

RMG_RUN_CHECK_CPU(${testname} ${num_proc})


