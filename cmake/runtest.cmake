file(REMOVE_RECURSE "${WORKDIR}")
EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_DIR}" "${WORKDIR}" )

execute_process(COMMAND ${MPIEXEC} -n ${PROCS} ${RMG_EXE} input WORKING_DIRECTORY ${WORKDIR})
execute_process(COMMAND python3 ${SECCMD}  WORKING_DIRECTORY ${WORKDIR})
file(REMOVE_RECURSE "${WORKDIR}/Waves")
#execute_process(COMMAND mpirun -n 8 /home/luw/rmgdft/build/rmg-cpu input WORKING_DIRECTORY ${WORKDIR})
