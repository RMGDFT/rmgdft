cd rmg
mpirun -n 16 ~/rmgdft/build/rmg-cpu input
mpirun -n 16 ~/rmgdft/build/rmg-cpu input_band
cd ../rmg-on
mpirun -n 8 ~/rmgdft/build/rmg-on-cpu input0
mpirun -n 8 ~/rmgdft/build/rmg-on-cpu input
mpirun -n 8 ~/rmgdft/build/rmg-on-cpu input_band
