import os

from PyQt4 import QtGui, QtCore

from distutils.sysconfig import get_python_lib


def write_out_jobfiles(configuration, setup, grids):
    """
       write out job script files
    """

    if(setup._machine.currentText() == 'Summit'):

        jobcommon ="""
#!/bin/bash
# Begin LSF directives
#BSUB -P CHP107
#BSUB -J test
#BSUB -o test.o%J
#BSUB -W 0:30
#BSUB -alloc_flags smt1
#BSUB -nnodes 1
##BSUB -alloc_flags gpumps
# End LSF directives and begin shell commands
module load gcc
module load essl
module load boost
module load cuda
module load netlib-lapack
module load netlib-scalapack


export OMP_NUM_THREADS=7
export RMG_NUM_THREADS=7
export OMP_DYNAMIC=TRUE
export OMP_WAIT_POLICY=passive
#export GOMP_SPIN_COUNT=1000
jsrun -n1 -r1 -c1 --bind none ./rmg-cpu rmg.input
#jsrun -n48 -r6 -g1 -c7 -b none --latency_priority cpu-memory -gpu ./rmg-gpu in.cu256ft
"""
        with open('job','w') as inc_file:
            inc_file.write(jobcommon)
