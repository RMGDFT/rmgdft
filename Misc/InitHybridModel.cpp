/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "portability.h"
#include <stdlib.h>
#include <limits.h>
#include <cstring>
#include <mpi.h>
#include <omp.h>
#include <transition.h>

#ifdef USE_NUMA
    #include <numa.h>
    #include <numaif.h>
#endif


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

    // Determine if we are the master process in this group
    pct.is_local_master = false;
    MPI_Comm_rank(pct.local_comm, &pct.local_rank);
    if(pct.local_rank == 0) pct.is_local_master = true;

    // Create a communicator consisting of only the local masters on each node
    int color = MPI_UNDEFINED;
    if(pct.is_local_master) color = 1;
    MPI_Comm_split(comm, color, localrank, &pct.local_master_comm);

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
            std::cout << "RMG running with " << nthreads << " threads set manually via input file." << std::endl;

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
                std::cout << "RMG running with " << nthreads << " threads set via OMP_NUM_THREADS." << std::endl;

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
                std::cout << "RMG running with " << pct.procs_per_host << " MPI procs per host and " << nthreads << " threads per MPI proc set automatically." << std::endl;
            if(pct.worldrank == 0)
                std::cout << "OMP_NUM_THREADS environment variable was not set so using automatically determined value of " << nthreads << "." << std::endl;

        }

    }

#ifdef USE_NUMA
    // Determine if this is a numa system and setup up policies correctly if that is the case.
    // We need to optimize the interaction between MPI procs, threads and numa nodes.
    //
    // Case 1: A single MPI process per host.
    //         Main thread data structures prefer node interleaved memory allocation
    //         Worker thread data structures prefer node local memory allocation
    // Case 2: 1 MPI process per each core on a host
    //         Both main and worker threads prefer node local memory allocation 
    // Case 3: Number of MPI procs per host is equal to the number of numa nodes
    //         Both main and worker threads prefer node local memory allocation
    //         but we want to lock each proc and it's threads to separate nodes.
    // Case 4: Number of MPI_procs > numa_nodes
    //         Each MPI proc and it's threads should be locked to a specific node
    //         and distributed evenly.
    // Case 5: Number of MPI_procs < numa_nodes but not equal to 1
    //         In this case we want to restrict each process to a subset of the
    //         numa nodes based on distance. So for example with a two socket system
    //         with 2 nodes per socket and 2 MPI_procs we would want to run each proc
    //         on a separate socket with two numa nodes.
    
    if(ct.use_numa && (numa_available() < 0))
    {
        ct.use_numa = false;
    }

    if(ct.use_numa)
    {
        // Obviously assuming all CPU's active which is reasonable for HPC
        pct.numa_nodes_per_host = numa_max_node() + 1;
        pct.cpumask = numa_allocate_cpumask();
        pct.nodemask = numa_allocate_nodemask();
    }
    if(pct.numa_nodes_per_host == 1) ct.use_numa = false;

    if(ct.use_numa)
    {
        if(pct.procs_per_host == 1) {
            // Case 1
            bitmask *tmask = numa_allocate_cpumask();
            numa_bitmask_clearall(tmask);
            numa_bitmask_clearall(pct.cpumask);
            for(int nid = 0;nid < pct.numa_nodes_per_host;nid++)
            {
                numa_bitmask_setbit(pct.nodemask, nid);
                numa_node_to_cpus(nid, tmask);
                for(unsigned int idx=0;idx < tmask->size;idx++)
                {
                    if(numa_bitmask_isbitset(tmask, idx)) numa_bitmask_setbit(pct.cpumask, idx); 
                }
            }
            numa_set_interleave_mask(pct.nodemask);
            numa_bind(pct.nodemask);

            //printf("set_mempolicy ret=%d   %d\n",ret,pct.numa_nodes_per_host);
            if(pct.gridpe==0)
                printf("C1: Numa aware allocation with 1 MPI proc, %d cores and %d numa nodes per host.\n", pct.ncpus, pct.numa_nodes_per_host);

            numa_free_cpumask(tmask);

        }
        else if(pct.ncpus == pct.procs_per_host) {

            // Case 2
            unsigned int nid = pct.local_rank / (pct.ncpus / pct.numa_nodes_per_host);
            unsigned int procs_per_node = pct.procs_per_host / pct.numa_nodes_per_host;
            unsigned int proc_rank_in_node = pct.local_rank % procs_per_node;

            numa_node_to_cpus(nid, pct.cpumask);
            numa_bitmask_setbit(pct.nodemask, nid);
            numa_bind(pct.nodemask);
            numa_migrate_pages(getpid(), numa_all_nodes_ptr, pct.nodemask);
            // We have one MPI proc per cpu so set affinity based on a one to one mapping from low to high
            unsigned int offset = 0;
            for(unsigned int idx=0;idx<pct.cpumask->size;idx++) 
            {
                if(numa_bitmask_isbitset(pct.cpumask, idx))
                {
                    if(offset == proc_rank_in_node)
                    {
                        numa_bitmask_clearall(pct.cpumask);
                        numa_bitmask_setbit(pct.cpumask, idx); 
                        numa_sched_setaffinity(0, pct.cpumask);
                        if(ct.verbose) printf("Binding rank %d  to cpu %d\n", pct.local_rank, idx);
                        break;
                    }
                    offset++;
                }
            }
            if(pct.gridpe==0)
                printf("C2: Numa aware allocation with %d MPI procs, %d cores and %d numa nodes per host.\n", pct.procs_per_host, pct.ncpus, pct.numa_nodes_per_host);
        }
        else if(pct.procs_per_host == pct.numa_nodes_per_host) {

            // Case 3
            unsigned int nid = pct.local_rank;
            numa_node_to_cpus(nid, pct.cpumask);

            numa_bitmask_setbit(pct.nodemask, nid); 
            numa_bind(pct.nodemask);
            numa_migrate_pages(getpid(), numa_all_nodes_ptr, pct.nodemask);

            if(pct.gridpe==0)
                printf("C3: Numa aware allocation with %d MPI procs, %d cores and %d numa nodes per host.\n", pct.procs_per_host, pct.ncpus, pct.numa_nodes_per_host);
        }
        else if((pct.procs_per_host > pct.numa_nodes_per_host) && (pct.ncpus != pct.procs_per_host)) 
        {
            // Case 4
            int nid = pct.local_rank / (pct.procs_per_host / pct.numa_nodes_per_host);
            numa_bitmask_setbit(pct.nodemask, nid);
            numa_bind(pct.nodemask);
            numa_node_to_cpus(nid, pct.cpumask);
            numa_migrate_pages(getpid(), numa_all_nodes_ptr, pct.nodemask);
            if(pct.gridpe == 0)
                printf("C4: Numa aware allocation with %d MPI procs, %d cores and %d numa nodes per host.\n", pct.procs_per_host, pct.ncpus, pct.numa_nodes_per_host);
        }
    }
#endif

    ct.THREADS_PER_NODE = nthreads;
    ct.MG_THREADS_PER_NODE = nthreads;

    // Check if RMG_NUM_THREADS was set?
    tptr = getenv("RMG_NUM_THREADS");
    if(tptr) ct.MG_THREADS_PER_NODE = atoi(tptr);
    if((ct.MG_THREADS_PER_NODE < 1) || (ct.MG_THREADS_PER_NODE > ct.THREADS_PER_NODE))
    {
        printf("Warning: RMG_NUM_THREADS set to invalid value. Resetting to default.\n");
        ct.MG_THREADS_PER_NODE = ct.THREADS_PER_NODE;
    }


    B = BaseThread::getBaseThread(nthreads);
    B->RegisterThreadFunction(run_threads, pct.grid_comm);

}

