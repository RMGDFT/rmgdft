import os

from PyQt4 import QtGui, QtCore

from distutils.sysconfig import get_python_lib


def write_out_jobfiles(configuration, setup, grids):
    """
       write out job script files
    """

    elements = [ ii for ii in configuration.conf[0].elements]
    element_types = list(set(elements))

    if(setup._machine.currentText() == 'ox'):
        job1 = ''
        job1 += 'cp ~/SVN/codes/ON/rmg-on .\n'
        job1 += 'cp ~/SVN/codes/NEGF/rmg-negf .\n'
        for i in range(len(element_types)):
            dir_name = element_types[i] + '-atom'
            job1 += 'cd ' + dir_name + '\n'
            job1 += 'mpirun -n 1 ../rmg-on input&\n'
            job1 += 'cd ..\n'
        with open('job1','w') as inc_file:
            inc_file.write(job1)
            mode = 700

        job2 = ''
        job2 += 'cd lead1\n'
        job2 += 'mpirun -n %d ../rmg-on input&\n'%grids.num_proc_left
        job2 += 'cd ..\n'
        job2 += 'cd lead2\n'
        job2 += 'mpirun -n %d ../rmg-on input&\n'%grids.num_proc_right
        job2 += 'cd ..\n'
        job2 += 'cd center\n'
        job2 += 'mpirun -n %d ../rmg-on input&\n'%grids.num_proc_center
        job2 += 'cd ..\n'

        with open('job2','w') as inc_file:
            inc_file.write(job2)

        job3 = ''
        job3 += 'cd 3lead_lead1\n'
        job3 += 'fermi=$(grep FERMI ../lead1/*.log|tail -1 | cut -d"=" -f2)\n'
        job3 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat*\n'
        job3 += 'mpirun -n %d ../rmg-negf input&\n'%grids.num_proc_left3
        job3 += 'cd ..\n'
        job3 += 'cd 3lead_lead2\n'
        job3 += 'fermi=$(grep FERMI ../lead2/*.log|tail -1 | cut -d"=" -f2)\n'
        job3 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat*\n'
        job3 += 'mpirun -n %d ../rmg-negf input&\n'%grids.num_proc_right3
        job3 += 'cd ..\n'

        with open('job3','w') as inc_file:
            inc_file.write(job3)

        job4 = ''
        job4 += 'cd bias_0.0\n'
        job4 += 'fermi=$(grep FERMI ../3lead_lead1/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat1\n'
        job4 += 'fermi=$(grep FERMI ../3lead_lead2/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat2\n'
        job4 += 'fermi=$(grep FERMI ../center/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat0\n'
        job4 += 'mpirun -n %d ../rmg-negf input&\n'%grids.num_proc_negf
        job4 += 'cd ..\n'

        with open('job4','w') as inc_file:
            inc_file.write(job4)

        os.chmod('job1',0700)
        os.chmod('job2',0700)
        os.chmod('job3',0700)
        os.chmod('job4',0700)

    all_in_one ='cp ~/SVN/codes/ON/rmg-on .\n'
    all_in_one +='cp ~/SVN/codes/NEGF/rmg-negf .\n'
    all_in_one +='myid1=$(qsub job1)\n'
    all_in_one +='myid2=$(qsub -W depend=afterok:$myid1 job2)\n'
    all_in_one +='myid3=$(qsub -W depend=afterok:$myid2 job3)\n'
    all_in_one +='myid4=$(qsub -W depend=afterok:$myid3 job4)\n'
    with open('submit_all','w') as inc_file:
        inc_file.write(all_in_one)
    os.chmod('submit_all',0700)

    if(setup._machine.currentText() == 'Jaguar'):
        # cores per node
        cpn = 16
        
        tot_size = len(element_types) * cpn
        jobcommon = '#!/bin/bash\n'
        jobcommon += '#PBS -q %s\n'%setup._queue.currentText()
        jobcommon += '#PBS -j oe\n'
        jobcommon += '#PBS -A %s\n' %setup._projname.currentText()

        job1 = jobcommon
        job1 += '#PBS -N single_atom\n'
        job1 += '#PBS -l walltime=1:00:00,size=%d,gres=widow2\n\n'%tot_size
        job1 += 'cd $PBS_O_WORKDIR\n'

        for i in range(len(element_types)):
            dir_name = element_types[i] + '-atom'
            job1 += 'cd ' + dir_name + '\n'
            job1 += 'aprun -n 1 ../rmg-on input&\n'
            job1 += 'cd ..\n'

        job1 += 'wait\n'
        with open('job1','w') as inc_file:
            inc_file.write(job1)

        tot_size = (grids.num_proc_left-1)/cpn* cpn + cpn
        tot_size += (grids.num_proc_right-1)/cpn*cpn + cpn
        tot_size += (grids.num_proc_center-1)/cpn*cpn + cpn
        job2 = jobcommon
        job2 += '#PBS -N order-n\n'
        job2 += '#PBS -l walltime=4:00:00,size=%d,gres=widow2\n\n'%tot_size
        job2 += 'cd $PBS_O_WORKDIR\n'
        job2 += 'cd lead1\n'
        job2 += 'aprun -n %d ../rmg-on input&\n'%grids.num_proc_left
        job2 += 'cd ..\n'
        job2 += 'cd lead2\n'
        job2 += 'aprun -n %d ../rmg-on input&\n'%grids.num_proc_right
        job2 += 'cd ..\n'
        job2 += 'cd center\n'
        job2 += 'aprun -n %d ../rmg-on input&\n'%grids.num_proc_center
        job2 += 'cd ..\n'
        job2 += 'wait\n'

        with open('job2','w') as inc_file:
            inc_file.write(job2)

        tot_size  = (grids.num_proc_left3-1)/cpn * cpn + cpn
        tot_size += (grids.num_proc_right3-1)/cpn *cpn + cpn
        job3 = jobcommon
        job3 += '#PBS -N lead3negf\n'
        job3 += '#PBS -l walltime=4:00:00,size=%d,gres=widow2\n\n'%tot_size
        job3 += 'cd $PBS_O_WORKDIR\n'
        job3 += 'cd 3lead_lead1\n'
        job3 += 'fermi=$(grep FERMI ../lead1/*.log|tail -1 | cut -d"=" -f2)\n'
        job3 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat*\n'
        job3 += 'aprun -n %d ../rmg-negf input&\n'%grids.num_proc_left3
        job3 += 'cd ..\n'
        job3 += 'cd 3lead_lead2\n'
        job3 += 'fermi=$(grep FERMI ../lead2/*.log|tail -1 | cut -d"=" -f2)\n'
        job3 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat*\n'
        job3 += 'aprun -n %d ../rmg-negf input&\n'%grids.num_proc_right3
        job3 += 'cd ..\n'
        job3 += 'wait\n'

        with open('job3','w') as inc_file:
            inc_file.write(job3)

        tot_size = grids.num_proc_negf
        tot_size = (tot_size -1)/cpn*cpn + cpn
        job4 = jobcommon
        job4 += '#PBS -N zerobias\n'
        job4 += '#PBS -l walltime=4:00:00,size=%d,gres=widow2\n\n'%tot_size
        job4 += 'cd $PBS_O_WORKDIR\n'
        job4 += 'cd bias_0.0\n'
        job4 += 'fermi=$(grep FERMI ../3lead_lead1/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat1\n'
        job4 += 'fermi=$(grep FERMI ../3lead_lead2/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat2\n'
        job4 += 'fermi=$(grep FERMI ../center/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat0\n'
        job4 += 'aprun -n %d ../rmg-negf input&\n'%grids.num_proc_negf
        job4 += 'cd ..\n'
        job4 += 'wait\n'

        with open('job4','w') as inc_file:
            inc_file.write(job4)


    if(setup._machine.currentText() == 'chugach'):
        # cores per node
        cpn = 16
        jobcommon = '#!/bin/bash\n'
        jobcommon += '#PBS -q %s\n'%setup._queue.currentText()
        jobcommon += '#PBS -j oe\n'
        jobcommon += '#PBS -A %s\n' %setup._projname.currentText()

        job1 = jobcommon
        tot_size = len(element_types) * cpn
        job1 += '#PBS -N single_atom\n'
        job1 += '#PBS -l walltime=1:00:00,mppwidth=%d\n\n'%tot_size
        job1 += 'cd $PBS_O_WORKDIR\n'

        for i in range(len(element_types)):
            dir_name = element_types[i] + '-atom'
            job1 += 'cd ' + dir_name + '\n'
            job1 += 'aprun -n 1 ../rmg-on input&\n'
            job1 += 'cd ..\n'

        job1 += 'wait\n'
        with open('job1','w') as inc_file:
            inc_file.write(job1)

        tot_size = (grids.num_proc_left-1)/cpn* cpn + cpn
        tot_size += (grids.num_proc_right-1)/cpn*cpn + cpn
        tot_size += (grids.num_proc_center-1)/cpn*cpn + cpn
        job2 = jobcommon
        job2 += '#PBS -N order-n\n'
        job2 += '#PBS -l walltime=4:00:00,mppwidth=%d\n\n'%tot_size
        job2 += 'cd $PBS_O_WORKDIR\n'
        job2 += 'cd lead1\n'
        job2 += 'aprun -n %d ../rmg-on input&\n'%grids.num_proc_left
        job2 += 'cd ..\n'
        job2 += 'cd lead2\n'
        job2 += 'aprun -n %d ../rmg-on input&\n'%grids.num_proc_right
        job2 += 'cd ..\n'
        job2 += 'cd center\n'
        job2 += 'aprun -n %d ../rmg-on input&\n'%grids.num_proc_center
        job2 += 'cd ..\n'
        job2 += 'wait\n'

        with open('job2','w') as inc_file:
            inc_file.write(job2)

        tot_size  = (grids.num_proc_left3-1)/cpn * cpn + cpn
        tot_size += (grids.num_proc_right3-1)/cpn *cpn + cpn
        job3 = jobcommon
        job3 += '#PBS -N lead3negf\n'
        job3 += '#PBS -l walltime=4:00:00,mppwidth=%d\n\n'%tot_size
        job3 += 'cd $PBS_O_WORKDIR\n'
        job3 += 'cd 3lead_lead1\n'
        job3 += 'fermi=$(grep FERMI ../lead1/*.log|tail -1 | cut -d"=" -f2)\n'
        job3 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat*\n'
        job3 += 'aprun -n %d ../rmg-negf input&\n'%grids.num_proc_left3
        job3 += 'cd ..\n'
        job3 += 'cd 3lead_lead2\n'
        job3 += 'fermi=$(grep FERMI ../lead2/*.log|tail -1 | cut -d"=" -f2)\n'
        job3 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat*\n'
        job3 += 'aprun -n %d ../rmg-negf input&\n'%grids.num_proc_right3
        job3 += 'cd ..\n'
        job3 += 'wait\n'

        with open('job3','w') as inc_file:
            inc_file.write(job3)

        tot_size = grids.num_proc_negf
        tot_size = (tot_size -1)/cpn*cpn + cpn
        job4 = jobcommon
        job4 += '#PBS -N zerobias\n'
        job4 += '#PBS -l walltime=4:00:00,mppwidth=%d\n\n'%tot_size
        job4 += 'cd $PBS_O_WORKDIR\n'
        job4 += 'cd bias_0.0\n'
        job4 += 'fermi=$(grep FERMI ../3lead_lead1/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat1\n'
        job4 += 'fermi=$(grep FERMI ../3lead_lead2/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat2\n'
        job4 += 'fermi=$(grep FERMI ../center/*.log|tail -1 | cut -d"=" -f2)\n'
        job4 += 'sed -i "s/XXXXX/$fermi/g" LCR.dat0\n'
        job4 += 'aprun -n %d ../rmg-negf input&\n'%grids.num_proc_negf
        job4 += 'cd ..\n'
        job4 += 'wait\n'

        with open('job4','w') as inc_file:
            inc_file.write(job4)

