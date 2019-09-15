file(REMOVE_RECURSE ${WORKDIR}/input.00.log ${WORKDIR}/input.00.options ${WORKDIR}/input.00.rmsdv.xmgr ${WORKDIR}/Waves)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_EXE} input WORKING_DIRECTORY ${WORKDIR})
#execute_process(COMMAND mpirun -n 8 /home/luw/rmgdft/build/rmg-cpu input WORKING_DIRECTORY ${WORKDIR})
