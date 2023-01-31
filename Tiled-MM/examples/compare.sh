sizes=(4000 6000 8000 10000 12000 14000 16000 18000 20000)
n_rep=2

for size in ${sizes[@]}
do
    echo "Tiled MM:"
    srun -N1 ./examples/multiply -m $size -n $size -k $size -r $n_rep
    echo "cublasXt:"
    srun -N1 ./examples/cublasXt-multiply -m $size -n $size -k $size -r $n_rep
done
