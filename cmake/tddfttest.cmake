file(REMOVE_RECURSE ${WORKDIR}/input.00.log  ${WORKDIR}/input.00.option  
    ${WORKDIR}/input.00_spin0_current.dat
    ${WORKDIR}/Waves)
execute_process(COMMAND mpirun -n ${PROCS} ${RMG_EXE} input WORKING_DIRECTORY ${WORKDIR})
execute_process(COMMAND python3 ${CHECK_TDDFT}  WORKING_DIRECTORY ${WORKDIR})
