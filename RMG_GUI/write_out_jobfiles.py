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
#BSUB -J cmd
#BSUB -o cmd.o%J
#BSUB -W 0:10
#BSUB -nnodes 2
#BSUB -alloc_flags "smt1"
#BSUB -alloc_flags gpumps
# End LSF directives and begin shell commands
#
#
#  RMG job script for Summit
#
#  RMG can be run in a hybrid mode where multiple CPU cores are
#  assigned to each MPI task as well as a standard mode where only one
#  CPU core is assigned to each MPI task. In the hybrid case threads
#  are used to maximize utilization of the cores assigned to an MPI task.
# 
#  The preferred mode depends on both the problem size and the number
#  of Summit nodes allocated #  for the job. For smaller problems/jobs the 
#  standard mode usually provides higher performance while for very large
#  jobs the hybrid mode is preferred.
#
#
#  Each summit node has 6 GPU's and 42 CPU cores available for applications.
#  These may be grouped into resource sets. We have been running RMG using
#  6 resource sets per node with each resource set allocated 1 GPU and 7 CPU
#  cores. In standard mode there are 6 MPI tasks per resource set with 1 thread
#  per task. In hybrid mode there is 1 MPI task per resource set with 6 threads
#  per task. Other combinations are possible and may yield better performance
#  in certain cases and this will be an area of investigation going forward.
#
#  Load modules. Use the same modules that were used to build the executable
module load gcc/6.4.0
module load boost
module load essl
module load cuda
module load fftw
module load cmake
module load netlib-scalapack/2.0.2
module unload darshan

#  No limits
ulimit -c 0

#  rmg binary in project shared directory
rmgbin=/gpfs/alpine/chp107/proj-shared/elbriggs/Tutorial/bin/rmg-gpu

#  set OMP policies
export OMP_DYNAMIC=true
export OMP_WAIT_POLICY=passive


#  Set thread counts for standard mode
export OMP_NUM_THREADS=1
export RMG_NUM_THREADS=1

#  Run using 6 MPI tasks per resource set for a total of 72 MPI tasks
jsrun -n12 -a6 -g1 -r6 -c7 -d plane:6 --latency_priority cpu-memory --smpiargs "-gpu" $rmgbin rmg.input

"""
        with open('job','w') as inc_file:
            inc_file.write(jobcommon)
