#include <unistd.h>
#include <stdlib.h>
#include <limits.h>
#include <cstring>
#include <mpi.h>
#include <omp.h>
#include <transition.h>


// Tries to determine system information that can be used to setup optimal threading.
// We assume that if the user has set OMP_NUM_THREADS manually that they know what
// they are doing so we don't try to override their choices. We will issue a warning if
// they do something that doesn't make sense though.
int GetRecommendedThreadNumber(int nthreads, int npes, int thispe, MPI_Comm comm)
{


    // If user has set nthreads manually then we don't try to adjust anything
    if(nthreads != -1) {
        std::cout << "Running with " << nthreads << " threads set manually via input file." << std::endl;
        return nthreads;
    }


    int ncpus, name_len;
    // Determine how many MPI procs there are per host
    int stride = MPI_MAX_PROCESSOR_NAME+2;
    char *hnames = new char[stride * npes]();
    MPI_Get_processor_name(&hnames[thispe * stride], &name_len);

    MPI_Allgather(&hnames[thispe * stride], stride, MPI_CHAR, hnames, stride, MPI_CHAR, comm);

    int procs_per_host = 0;
    for(int i = 0;i < npes;i++) {
        if(!std::strcmp(&hnames[thispe * stride], &hnames[i * stride])) procs_per_host++;
    }
    
    // Sysconf works on linux, hopefully C++11 works elsewhere
#if __linux__
    ncpus = sysconf( _SC_NPROCESSORS_ONLN );
#else
    ncpus = std::thread::hardware_concurrency();
#endif

    if(ncpus >= procs_per_host) {
        nthreads = ncpus / procs_per_host;
    }
    else {
        nthreads = 1;
    }
    delete [] hnames; 

    // If OMP_NUM_THREADS is not set explictly we set it
    char *tptr = getenv("OMP_NUM_THREADS");
    if(!tptr) {
        omp_set_num_threads(nthreads);
    }

    std::cout << "Running with " << procs_per_host << " MPI procs per host and " << nthreads << " threads per MPI proc set automatically." << std::endl;

    // Probably sub optimal but will figure something out later for non-linux systems
    nthreads = 1;

    return nthreads;
}

