@ECHO OFF
doskey mpirun="C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" $1 $2 $3 $4
set RMGBINARY="C:\Program Files\rmg\bin\rmg.exe"
set RMGEXAMPLES="C:\Program Files\rmg\Examples"
cls
mode con: cols=120 lines=40
@ECHO:
@ECHO RMG v1.2 command line window. RMG requires Microsoft MPI v5 in order to run. If MS-MPI
@ECHO is installed in it's standard location you can run the sample programs using the
@ECHO following command.
@ECHO:
@ECHO        mpirun -np 4 %%RMGBINARY%% input_file
@ECHO:
@ECHO where the number of MPI processes to use (-np 4) may be adjusted depending on your
@ECHO hardware platform. The input_file should be located in your current directory or a
@ECHO subdirectory of the current directory and you must have write permissions for it.
@ECHO:
@ECHO Example input files may be obtained from 
@ECHO:
@ECHO     http://sourceforge.net/projects/rmgdft/files/Releases/
@ECHO:
@ECHO or from the Examples subdirectory of your RMG installation.
@ECHO
cmd.exe
exit /B 0
