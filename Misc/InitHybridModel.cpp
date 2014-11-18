#include <unistd.h>
#include <stdlib.h>
#include <limits.h>
#include <cstring>
#include <mpi.h>
#include <omp.h>
#include <transition.h>

void *run_threads(void *v);
static BaseThread *B;

// Determine system information required to to setup optimal threading and local
// MPI communications. We assume that if the user has set OMP_NUM_THREADS manually
// that they know what they are doing so we don't try to override their choices.
// We will issue a warning if // they do something that doesn't make sense though.
void InitHybridModel(int nthreads, int npes, int thispe, MPI_Comm comm)
{

    int omp_num_threads_set = -1;
    MPI_Group parent, localgroup;

    // Determine hardware resources and how many MPI procs there are per host
    int name_len, localrank;
    int stride = MPI_MAX_PROCESSOR_NAME+2;
    char *hnames = new char[stride * npes]();

    int *ranks = new int[npes];
    MPI_Comm_rank(comm, &localrank);
    ranks[thispe] = localrank;

    MPI_Get_processor_name(&hnames[thispe * stride], &name_len);

    MPI_Allgather(MPI_IN_PLACE, stride, MPI_CHAR, hnames, stride, MPI_CHAR, comm);
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, ranks, 1, MPI_INT, comm);

    pct.procs_per_host = 0;
    for(int i = 0;i < npes;i++) {
        if(!std::strcmp(&hnames[thispe * stride], &hnames[i * stride])) {
            pct.procs_per_host++;
        }
        else {
            ranks[i] = -1;
        }
    }


    pct.mpi_local_ranks = new int[pct.procs_per_host];
    int j = 0;
    for(int i = 0;i < npes;i++) {
        if(ranks[i] >= 0) {
            pct.mpi_local_ranks[j] = ranks[i];
            j++;
        }
    }

    // Next create the local group and a communicator to go with it
    MPI_Comm_group(comm, &parent);
    MPI_Group_incl(parent, pct.procs_per_host, pct.mpi_local_ranks, &localgroup);
    MPI_Comm_create(comm, localgroup, &pct.local_comm);

    // And finally determine if we are the master process in this group
    pct.is_local_master = false;
    MPI_Comm_rank(pct.local_comm, &pct.local_rank);
    if(pct.local_rank == 0) pct.is_local_master = true;

    delete [] ranks;
    delete [] hnames; 

    // Sysconf works on linux, hopefully C++11 works elsewhere
#if __linux__
    pct.ncpus = sysconf( _SC_NPROCESSORS_ONLN );
#else
    pct.ncpus = std::thread::hardware_concurrency();
#endif


    // Check if OMP_NUM_THREADS was set?
    char *tptr = getenv("OMP_NUM_THREADS");
    if(tptr) omp_num_threads_set = atoi(tptr);

    // If user has set nthreads manually then we don't try to adjust it
    if(nthreads > 0) {

        if(pct.worldrank == 0)
            std::cout << "Running with " << nthreads << " threads set manually via input file." << std::endl;

        // User set threads in input file but did not set OMP_NUM_THREADS so use input file value
        if(omp_num_threads_set < 0) { 

            omp_set_num_threads(nthreads);
            if(pct.worldrank == 0)
                std::cout << "OMP_NUM_THREADS environment variable was not set so using input file value for threads_per_node." << std::endl;

        }
        else {

            // Both input file and OMP_NUM_THREADS set so check if they are equal and issue warning if not
            if(nthreads != omp_num_threads_set) {
                if(pct.worldrank == 0)
                    std::cout << "Warning: OMP_NUM_THREADS = " << omp_num_threads_set << " and threads_per_node = " << nthreads << "." << std::endl;
                if(pct.worldrank == 0)
                    std::cout << "This may be OK but just checking if that is what you intended." << std::endl;
            }

        }

    }
    else {

        if(omp_num_threads_set > 0) {

            // set nthreads to OMP_NUM_THREADS
            nthreads = omp_num_threads_set;
            if(pct.worldrank == 0)
                std::cout << "Running with " << nthreads << " threads set via OMP_NUM_THREADS." << std::endl;

        }
        else {

            if(pct.ncpus >= pct.procs_per_host) {
                nthreads = pct.ncpus / pct.procs_per_host;
            }
            else {
                nthreads = 1;
            }

            omp_set_num_threads(nthreads);
            if(pct.worldrank == 0)
                std::cout << "Running with " << pct.procs_per_host << " MPI procs per host and " << nthreads << " threads per MPI proc set automatically." << std::endl;
            if(pct.worldrank == 0)
                std::cout << "OMP_NUM_THREADS environment variable was not set so using automatically determined value of " << nthreads << "." << std::endl;

        }

    }


    ct.THREADS_PER_NODE = nthreads;
    B = BaseThread::getBaseThread(nthreads);
    B->RegisterThreadFunction(run_threads);

}

