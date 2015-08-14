set(CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH};${RMG_SOURCE_DIR}/Headers/my_scalapack.h")
set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES};${MPI_C_LIBRARIES};${SCALAPACK_LIBRARIES};${BLAS_LIBRARIES};${BLACS_LIBRARIES};${BLACSCINIT_LIBRARIES}")
CHECK_C_SOURCE_COMPILES("#include <mpi.h>
                           int main(int argc, char **argv) {
                           int N=1000, NB=16, myrow, izero=0, nprow;
                           numroc_ (&N, &NB, &myrow, &izero, &nprow);
                         }" RMG_SCALAPACK_INTERNAL)

