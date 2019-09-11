RMG_EXEC=/home/luw/rmgdft/build/rmg-cpu
RMG_ON_EXEC=/home/luw/rmgdft/build/rmg-on-cpu
QMCPACK_PATH=/home/luw/qmcpack
CURRENT_PATH=$(pwd)
MPIEXEC=/usr/lib64/mpi/gcc/openmpi/bin/mpiexec

if [ ! -f $RMG_EXEC ]; then
    echo "$RMG_EXEC does not exist"
    exit 1
fi

if [ ! -f $RMG_ON_EXEC ]; then
    echo "$RMG_ON_EXEC does not exist"
    exit 2
fi

if [ ! -f $QMCPACK_PATH/build/bin/qmcpack ]; then
    echo "$RMG_ON_EXEC does not exist"
    exit 3
fi

rm -rf $CURRENT_PATH/atomO_pp
tar xvf atomO_pp.tar
cd $CURRENT_PATH/atomO_pp


cd $CURRENT_PATH/atomO_pp/RMG
$MPIEXEC "-n" "8" $RMG_EXEC input &> rmg.log
cd ..
ln -s $CURRENT_PATH/atomO_pp/RMG/Waves/wave.out.h5 $CURRENT_PATH/atomO_pp/RMG_noj/atomO.rmg.h5
ln -s $CURRENT_PATH/atomO_pp/RMG/Waves/wave.out.h5 $CURRENT_PATH/atomO_pp/RMG_sdj/atomO.rmg.h5


cd $CURRENT_PATH/atomO_pp/RMG_sdj
$MPIEXEC "-n" "1" "$QMCPACK_PATH/build/bin/qmcpack" "qmc_short.in.xml" &>qmc.log

echo "Test for RMG/QMCPACK interface for sdj" 
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--ke" "11.85233 0.04"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--le" "-15.870035 0.0045"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--ee" "10.368642 0.014"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--ts" "64000 0.0"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--lpp" "-40.15639 0.057"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--nlpp" "2.065389 0.020"

echo "Test for RMG/QMCPACK interface for noj"
cd $CURRENT_PATH/atomO_pp/RMG_noj
 $MPIEXEC "-n" "1" "/home/luw/qmcpack/build/bin/qmcpack" "qmc_short_noj.in.xml" &>qmc.log
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short_noj" "-e" "2" "--ke" "11.64798 0.04"

cd $CURRENT_PATH/atomO_pp/RMG_ON
$MPIEXEC "-n" "8" $RMG_ON_EXEC input &> rmg.log
cd ..
ln -s $CURRENT_PATH/atomO_pp/RMG_ON/Waves/wave.out.h5 $CURRENT_PATH/atomO_pp/RMG_ON_noj/atomO.rmg.h5
ln -s $CURRENT_PATH/atomO_pp/RMG_ON/Waves/wave.out.h5 $CURRENT_PATH/atomO_pp/RMG_ON_sdj/atomO.rmg.h5


cd $CURRENT_PATH/atomO_pp/RMG_ON_sdj
$MPIEXEC "-n" "1" "$QMCPACK_PATH/build/bin/qmcpack" "qmc_short.in.xml" &>qmc.log


echo "Test for RMG_ON/QMCPACK interface for sdj"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--ke" "11.85233 0.04"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--le" "-15.870035 0.0045"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--ee" "10.368642 0.014"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--ts" "64000 0.0"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--lpp" "-40.15639 0.057"
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short" "-e" "2" "--nlpp" "2.065389 0.020"

echo "Test for RMG_ON/QMCPACK interface for noj"
cd $CURRENT_PATH/atomO_pp/RMG_ON_noj
 $MPIEXEC "-n" "1" "/home/luw/qmcpack/build/bin/qmcpack" "qmc_short_noj.in.xml" &>qmc.log
 "$QMCPACK_PATH/tests/scripts/check_scalars.py" "--ns" "3" "--series" "0" "-p" "qmc_short_noj" "-e" "2" "--ke" "11.64798 0.04"
