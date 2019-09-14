SET(testname "atomO_polarized") 
SET(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${testname}")

COPY_DIRECTORY( "${SRC_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${testname}" )

RMG_RUN_CHECK_CPU(${testname} ${num_proc})


