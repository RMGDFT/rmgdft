file(REMOVE_RECURSE ${WORKDIR}/*/input.00.log ${WORKDIR}/*/input.00.options 
    ${WORKDIR}/*/input.00.rmsdv.xmgr ${WORKDIR}/*/Waves
    ${WORKDIR}/*/input.110.00.log ${WORKDIR}/*/input.110.00.options 
    ${WORKDIR}/*/input.110.00.cond_12.dat )
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_ON} input WORKING_DIRECTORY ${WORKDIR}/lead1)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_ON} input WORKING_DIRECTORY ${WORKDIR}/center)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_NEGF} input WORKING_DIRECTORY ${WORKDIR}/3lead_lead1)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_NEGF} input.110 WORKING_DIRECTORY ${WORKDIR}/3lead_lead1)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_NEGF} input WORKING_DIRECTORY ${WORKDIR}/bias_0.0)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_NEGF} input.110 WORKING_DIRECTORY ${WORKDIR}/bias_0.0)
execute_process(COMMAND python3 ${CHECK_COND}  WORKING_DIRECTORY ${WORKDIR}/bias_0.0)
