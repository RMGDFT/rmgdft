set(CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
set(CMAKE_REQUIRED_LIBRARIES "${MPI_C_LIBRARIES}")
CHECK_C_SOURCE_COMPILES("#include <mpi.h>
                           int main(int argc, char **argv) {
                           int tarray[100], n=100;
                           MPI_Comm  comm;
                           MPI_Request MPI_reqeigs;
                           MPI_Iallreduce(MPI_IN_PLACE, tarray, n, MPI_INT, MPI_SUM, comm, &MPI_reqeigs);                           
                         }" RMG_ASYNC_ALLREDUCE_INTERNAL)

