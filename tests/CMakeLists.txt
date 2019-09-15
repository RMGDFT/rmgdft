MESSAGE(STATUS "Creating tests for RMG")
if (NOT num_proc)
SET(num_proc 16)
endif()

SUBDIRS(RMG)
SUBDIRS(RMG_ON)


