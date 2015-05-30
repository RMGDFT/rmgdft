@ECHO OFF
cd %HOMEPATH%
doskey mpirun="C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" $1 $2 $3 $4
set RMGBINARY="C:\Program Files\rmg\bin\rmg.exe"
cls
@ECHO:
@ECHO RMG command line window. RMG requires Microsoft MPI in order to run. If MS-MPI
@ECHO is installed in it's standard location you can run the sample programs using the
@ECHO the following command.
@ECHO:
@ECHO        mpirun -np 4 %RMGBINARY% input_file
@ECHO:
@ECHO where the number of MPI processes to use (-np 4) may be adjusted depending on your
@ECHO hardware platform. The input_file should be located in your current directory or
@ECHO a subdirectory that you have write permissions for.
exit /B 0
