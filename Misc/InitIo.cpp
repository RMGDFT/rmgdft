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
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#if (defined(_WIN32) || defined(_WIN64))
    #include <io.h>
#else
    #include <unistd.h>
#endif



#include "const.h"
#include "params.h"
#include "rmg_alloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "GlobalSums.h"
#include "RmgException.h"
#include "InputKey.h"
#include "InputOpts.h"
#include "Functional.h"

#if GPU_ENABLED
    #if MAGMA_LIBS
        #include <magma.h>
    #endif

    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>

#endif

#if LINUX
void get_topology(void)
{
  int rank,nproc;
  int core;
  cpu_set_t my_set;
  char hostname[256];

  setbuf(stdout, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  core=sched_getcpu();
  gethostname(hostname,255);

  CPU_ZERO(&my_set);
  sched_getaffinity(0, sizeof(my_set), &my_set);
  char str[1024];
  strcpy(str," ");
  int count = 0;
  int j;
  for (j = 0; j < CPU_SETSIZE; ++j)
  {
    if (CPU_ISSET(j, &my_set))
    {
      ++count;
      char cpunum[2];
      sprintf(cpunum, "%d ", j);
      strcat(str, cpunum);
    }
  }
  printf("Rank: %4d; NumRanks: %d; RankCore: %3d; Hostname: %s; Affinity: %3d\n", rank, nproc, core, hostname, count);
}
#endif


void InitIo (int argc, char **argv, std::unordered_map<std::string, InputKey *>& ControlMap)
{

    int npes, worldpe, status, provided=0;

    MPI_Init_thread(&argc, &argv, ct.mpi_threadlevel, &provided);

    /* get this cores mpi rank */
    MPI_Comm_rank (MPI_COMM_WORLD, &worldpe);
    pct.worldrank = worldpe;

    // Set error handler to only print to rank 0
    RmgErrorSetPrint(pct.worldrank == 0);

    /* get total mpi core count */
    MPI_Comm_size (MPI_COMM_WORLD, &npes);
    pct.total_npes = npes;

    if(argc == 2)
    {
        std::unordered_map<std::string, InputKey *> SetupMap;
        ReadInit(argv[1], ct, pct, SetupMap);

    }
    else {

        std::unordered_map<std::string, InputKey *> SetupMap;
        ReadInit("ctrl_init.dat", ct, pct, SetupMap);

    }

    InitPe4image();

    snprintf (ct.cfile, MAX_PATH, "%s%s", pct.image_path[pct.thisimg], pct.image_input[pct.thisimg]);
    snprintf (ct.shortname, MAX_PATH, "%s%s", pct.image_path[pct.thisimg], pct.image_input[pct.thisimg]);

    ReadCommon(argc, argv, ct.cfile, ct, pct, ControlMap);

    if(Verify("start_mode", "Restart From File", ControlMap) 
            || Verify("calculation_mode", "Band Structure Only", ControlMap) ) {
        std::string dynfile(ct.infile);
        dynfile = dynfile + ".restart";
        ReadDynamics((char *)dynfile.c_str(), ct, ControlMap);

        // Always use absolute coords in restart file
        InputKey *Ik = ControlMap["atomic_coordinate_type"];
        static std::string AbsoluteCoords("Absolute");
        Ik->Readstr = AbsoluteCoords;
        ct.crd_flag = 1;
    }
    else {
        ReadDynamics(ct.cfile, ct, ControlMap);
    }
    ReadPseudo(ct.num_species, ct, ControlMap);

    // If fine/coarse grid ratio is not set then autoset it. By default we
    // use 1 for norm conserving pseudopotentials and 2 for ultrasoft.
    ct.norm_conserving_pp = true;
    for(int isp = 0;isp < ct.num_species;isp++) {
        SPECIES *sp = &ct.sp[isp];
        if(!sp->is_norm_conserving) ct.norm_conserving_pp = false;
    }
    if(!ct.FG_RATIO)
    {
        ct.FG_RATIO = 2;
        if(ct.norm_conserving_pp) ct.FG_RATIO = 1;
    }

    if(ct.forceflag == BAND_STRUCTURE)
    {
        int special_klines = ReadKpointsBandstructure(ct.cfile, ct, ControlMap);
        if(special_klines == 0) 
            ReadKpoints(ct.cfile, ct, ControlMap);
    }
    else if((ct.kpoint_mesh[0] < 1) || (ct.kpoint_mesh[1] < 1) || (ct.kpoint_mesh[2] < 1) ) 
    {
        ReadKpoints(ct.cfile, ct, ControlMap);
    }
    else
    {
        init_kpoints(ct.kpoint_mesh, ct.kpoint_is_shift);
    }

    InitPe4kpspin();

    AutoSet(ct, pct, ControlMap);
    Rmg_G->set_rank(pct.gridpe, pct.grid_comm);
    SetLaplacian();

    InitHybridModel(ct.OMP_THREADS_PER_NODE, ct.MG_THREADS_PER_NODE, pct.grid_npes, pct.gridpe, pct.grid_comm);

    /* Next address grid coalescence. Grids are only coalesced in the x-coordinate. For example if the
       global grid is (96,96,96) and there are 512 MPI process's then a non-coalesced arrangment would
       be to have an (8,8,8) processor grid with the wavefunctions defined on each MPI process in a
       (12,12,12) grid. This breakdown though limits the number of multigrid levels and suffers from
       high ratios of communication/computation in the multigrid solver. An alternate coalesced grid
       setup would be a processor grid of 4*(8,4,4) where the coalesced wavefunction grids are now
       (12,24,24) but each MPI process only handles a subset of the wavefunctions. In particlar.
 
         PE_X = 0,4  handles orbitals 0,4,8,12,16 ...
         PE_X = 1,5  handles orbitals 1,5,9,13,17 ...
         PE_X = 2,6  handles orbitals 2,6,10,14,18 ...
         PE_X = 3,7  handles orbitals 3,7,11,15,19 ...
    
       Threads introduce an additional complication. With non-coalesced grids individual orbitals are assigned
       to a separate thread with each orbital having a thread assigned on each PE. Since orbitals are only assigned
       to a subset of the PE's using coalesced grids some convention needs to be followed. The chosen convention is
       illustrated below where there are 8 PEs in an (8,1,1) arrangement with 4 threads/PE and a coalesce factor of 2.
       Threads increment vertically and PE's horizontally in the diagram.

       PE/Thread|  0  |  1  |  2  |  3  |  4  |  5  |  6  |  7
       ----------------------------------------------------------
          0     | S0  | S4  | S0  | S4  | S0  | S4  | S0  | S4
          1     | S1  | S5  | S1  | S5  | S1  | S5  | S1  | S5
          2     | S2  | S6  | S2  | S6  | S2  | S6  | S2  | S6
          3     | S3  | S7  | S3  | S7  | S3  | S7  | S3  | S7

       The current implementation only handles coalesce factors up to 16. While it is possible to coalesce
       in more than one coordinate dimension that would require repacking the orbitals so for now this
       implementation is limited to the x-direction.

       Implementation details.
         The standard global TradeImages object can be used. The only difference is that the offset of the
         neighboring PE's in the x-direction is multiplied by the coalesce factor.

         We have to perform an MPI reduction in MgEigState in order to compute the eigenvalue for a given
         orbital so we need another communicator that only spans the required subset of MPI procs.

         The Gather and Scatter grid functions also need to pass data back and forth between MPI procs
         with the same y and z proc coordinates so we create a local coalesced communicator for that.

         The number of states must be an integral multiple of the coalesce factor.

         The numer of PE's in the x-direction must be evenly divisible by the coalesce factor.

*/

    // processors in x must be an integral multiple of coalesce factor and grid points in x must be evenly
    // divisible by processors in x. And for now the number of orbitals must be evenly divisible by the
    // number of coalesced grids so we adjust the orbital count later if required.
    int rem1 = pct.pe_x % pct.coalesce_factor;
    int rem2 = Rmg_G->get_NX_GRID(1) % pct.pe_x;
    if(ct.coalesce_states && !rem1 && !rem2)
    {
        if(pct.gridpe == 0)
            std::cout << "Notice: Coalescing states in X with factor " << pct.coalesce_factor << "." << std::endl;

        if((ct.MG_THREADS_PER_NODE > 1) && !ct.mpi_queue_mode)
        {
            ct.mpi_queue_mode = true;
            if(pct.gridpe == 0)
                std::cout << "Notice: Coalescing states require mpi_queue_mode=\"true\" if (nthreads>1). " << std::endl;
        }
        // Set up our coalesced grid communicator next
        int px, py, pz;
        Rmg_G->pe2xyz(pct.gridpe, &px, &py, &pz);
        int base_px = (px / pct.coalesce_factor) * pct.coalesce_factor;
        px = px % pct.coalesce_factor;
        MPI_Comm_split(pct.grid_comm, px + 1, pct.gridpe, &pct.coalesced_grid_comm);

        // And now our coalesced local communicator
        int color = Rmg_G->xyz2pe(0, py, pz);
        color += (base_px + 1)*100000;
        MPI_Comm_split(Rmg_G->comm, color, pct.gridpe, &pct.coalesced_local_comm);
    }
    else
    {
        if((pct.gridpe == 0) && ct.coalesce_states)
        {
            std::cout << "Warning: Coalesced grids cannot be used with these global and processor grids. " << std::endl;
        }
        ct.coalesce_states = false;
        pct.coalesce_factor = 1;
        pct.coalesced_grid_comm = pct.grid_comm;
    }

    // Now that coalescing is sorted out we need to check valid MG levels on the PE level (as opposed to
    // on the global level which was done in Autoset.cpp
    // Next check if the PE grids allow these levels
    int PX0_GRID = pct.coalesce_factor * Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);
    for(int checklevel = 1;checklevel < ct.eig_parm.levels;checklevel++)
    {
        int eig_level_err = false;
        if ((PX0_GRID / (1 << checklevel)) < 1) eig_level_err = true;
        if ((PY0_GRID / (1 << checklevel)) < 1) eig_level_err = true;
        if ((PZ0_GRID / (1 << checklevel)) < 1) eig_level_err = true;
        MPI_Allreduce(MPI_IN_PLACE, &eig_level_err, 1, MPI_INT, MPI_SUM, pct.grid_comm);
        if (eig_level_err) {
            ct.eig_parm.levels = checklevel - 1;
            if(pct.imgpe == 0) std::cout << "Too many eigenvalue multigrid levels specified. Resetting to " << ct.eig_parm.levels << std::endl;
            break;
        }
    }


    // We do not (currently) coalesce outside of the multigrid solver so the number of grid points on a PE in any
    // coordinate direction (non-coalesced) must be greater than or equal to the number of points used in the global FD routines.
    int fd_check_err = false;
    PX0_GRID = Rmg_G->get_PX0_GRID(1);
    if(PX0_GRID < ct.kohn_sham_fd_order/2) fd_check_err = true;
    if(PY0_GRID < ct.kohn_sham_fd_order/2) fd_check_err = true;
    if(PZ0_GRID < ct.kohn_sham_fd_order/2) fd_check_err = true;
    MPI_Allreduce(MPI_IN_PLACE, &fd_check_err, 1, MPI_INT, MPI_SUM, pct.grid_comm);
    if(fd_check_err) 
        rmg_error_handler (__FILE__, __LINE__, "The Number of grid points per PE must be >= kohn_sham_fd_order/2.\n");

    




    /* if logname exists, increment until unique filename found */
    if (pct.imgpe == 0)
    {
        int name_incr;
        name_incr = FilenameIncrement(ct.shortname);
        snprintf (ct.basename, sizeof(ct.basename) - 1, "%s.%02d", ct.shortname, name_incr);
        snprintf (ct.logname, sizeof(ct.logname) - 1, "%s.%02d.log", ct.shortname, name_incr);

        /* open and save logfile handle, printf is stdout before here */
        ct.logfile = fopen(ct.logname, "w");
    }
    else {
#if (defined(_WIN32) || defined(_WIN64))
        ct.logfile = fopen("NUL:", "w");
#else
        ct.logfile = fopen("/dev/null", "w");
#endif
    }

    MPI_Bcast(ct.logname, MAX_PATH, MPI_CHAR, 0, pct.img_comm);
    MPI_Comm_size (pct.img_comm, &status);
    rmg_printf ("RMG initialization ...");
    rmg_printf (" %d image(s) total, %d per node.", pct.images, ct.images_per_node);
    rmg_printf (" %d MPI processes/image. ", status);


#if GPU_ENABLED
    CUdevice dev;
    if(pct.local_rank == 0)
    {
        cudaDeviceReset();
        cuInit(0);
        cudaSetDeviceFlags(cudaDeviceScheduleAuto);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(pct.local_rank != 0)
    {
        cuInit(0);
    }
    cuDeviceGetCount( &ct.num_gpu_devices);
    size_t deviceMem;
    int clock;
    char name[1024];

    // Get device list and memory capacities. While this is not ideal under all circumstances
    // we will find the device with the largest memory and use all devices that have just as
    // much memory
    std::vector<size_t> device_mem;
    rmg_printf("\n");
    for(int idevice = 0; idevice < ct.num_gpu_devices; idevice++ ) {
        cuDeviceGet( &dev, idevice );
        cuDeviceGetName( name, sizeof(name), dev );
        cuDeviceTotalMem( &deviceMem, dev );
        ct.gpu_mem[idevice] = deviceMem;
        cuDeviceGetAttribute( &clock, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev );
        rmg_printf( "device %d: %s, %.1f MHz clock, %.1f MB memory\n", idevice, name, clock/1000.f, deviceMem/1024.f/1024.f );
        device_mem.push_back(deviceMem/1024.0/1024.0);
    }

    size_t max_mem = *std::max_element(device_mem.begin(), device_mem.end());
    ct.num_usable_gpu_devices = 0;
    int idevice = 0;
    int does_managed;
    ct.gpus_support_managed_memory = true;
    for(auto it = device_mem.begin();it != device_mem.end();++it) {
        if(*it >= 0.9*max_mem) {
            ct.gpu_device_ids[ct.num_usable_gpu_devices] = idevice;
            cuDeviceGet( &ct.cu_devices[ct.num_usable_gpu_devices], idevice);
            cuDeviceGetAttribute( &does_managed, CU_DEVICE_ATTRIBUTE_MANAGED_MEMORY, ct.cu_devices[ct.num_usable_gpu_devices]);
            if(!does_managed) ct.gpus_support_managed_memory = false;
            ct.num_usable_gpu_devices++;
        }
        idevice++;
    }

    if(ct.num_usable_gpu_devices == 0)
        rmg_error_handler (__FILE__, __LINE__, "No usable GPU devices were found on the system. Exiting.\n");

    // Now we have to decide how to allocate the GPU's to MPI procs if we have more than
    // one GPU/node.
    if((ct.num_usable_gpu_devices == pct.numa_nodes_per_host) && (pct.numa_nodes_per_host == pct.procs_per_host))
    {
        cudaSetDevice(ct.gpu_device_ids[pct.local_rank]);
        if( CUBLAS_STATUS_SUCCESS != cublasCreate(&ct.cublas_handle) ) {
            rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
        }
    }
    else if(pct.procs_per_host > ct.num_usable_gpu_devices)
    {
        // Round robin
        int next_gpu = 0;
        for(int rank = 0;rank < pct.procs_per_host;rank++)
        {
            if(rank == pct.local_rank)
            {
                cudaSetDevice(ct.gpu_device_ids[next_gpu]);
                if( CUBLAS_STATUS_SUCCESS != cublasCreate(&ct.cublas_handle) ) {
                    rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
                }
            }
            next_gpu++;
            next_gpu = next_gpu % ct.num_usable_gpu_devices;
        }
        
    }
    else
    {
        if( CUDA_SUCCESS != cuDeviceGet( &ct.cu_dev, 0 ) ) {
            rmg_error_handler (__FILE__, __LINE__, "CUDA: Cannot get the device\n");
        }
        cudaSetDevice(ct.cu_dev);
        if( CUBLAS_STATUS_SUCCESS != cublasCreate(&ct.cublas_handle) ) {
            rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
        }
    }

    cusolverStatus_t cusolver_status = cusolverDnCreate(&ct.cusolver_handle);
    if(cusolver_status != CUSOLVER_STATUS_SUCCESS)
    {
        rmg_error_handler (__FILE__, __LINE__, "cusolver initialization failed.\n");
    }
    cudaStreamCreateWithFlags(&ct.cusolver_stream, cudaStreamNonBlocking);
    cusolver_status = cusolverDnSetStream(ct.cusolver_handle, ct.cusolver_stream);
    if(cusolver_status != CUSOLVER_STATUS_SUCCESS)
    {
        rmg_error_handler (__FILE__, __LINE__, "cusolver stream initialization failed.\n");
    }


#if MAGMA_LIBS
    magma_init();
#endif

#endif

    // This is placed down here since the IO is not setup yet when provided is obtained above.
    if(provided < ct.mpi_threadlevel) {

        rmg_printf("Thread support requested = %d but only %d provided. Terminating.\n", ct.mpi_threadlevel, provided);
        MPI_Finalize();
        exit(0);

    }
    rmg_printf("Thread level %d.\n", provided);
    fflush(NULL);

    // Allocate storage for trade_images and global sums routines
    size_t elem_len = sizeof(std::complex<double>);
    if(ct.is_gamma) elem_len = sizeof(double);
    Rmg_Q = NULL;
    if((ct.MG_THREADS_PER_NODE > 1) && ct.mpi_queue_mode)
    {
#if USE_NUMA && USE_HWLOC
        Rmg_Q = new MpiQueue(64, ct.MG_THREADS_PER_NODE, pct.manager_cpumask, &pct.topology, ct.spin_manager_thread, ct.spin_worker_threads);
#elif USE_NUMA
        Rmg_Q = new MpiQueue(64, ct.MG_THREADS_PER_NODE, pct.manager_cpumask, NULL, ct.spin_manager_thread, ct.spin_worker_threads);
#elif USE_HWLOC
        Rmg_Q = new MpiQueue(64, ct.MG_THREADS_PER_NODE, NULL, &pct.topology, ct.spin_manager_thread, ct.spin_worker_threads);
#else
        Rmg_Q = new MpiQueue(64, ct.MG_THREADS_PER_NODE, NULL, NULL, ct.spin_manager_thread, ct.spin_worker_threads);
#endif
    }
    else
    {
        ct.mpi_queue_mode = false;
    }

    Rmg_T = new TradeImages(Rmg_G, elem_len, ct.mpi_queue_mode, Rmg_Q, pct.coalesce_factor);
    if(ct.verbose) Rmg_T->set_timer_mode(true);
    Rmg_T->set_MPI_comm(pct.grid_comm);

    GlobalSumsInit();

    // Check individual node sizes on all levels for poisson mg solver
    Mgrid MG(&Rmg_L, Rmg_T);
    int dx[MAX_MG_LEVELS];dx[0] = Rmg_G->get_PX0_GRID(ct.FG_RATIO);
    int dy[MAX_MG_LEVELS];dy[0] = Rmg_G->get_PY0_GRID(ct.FG_RATIO);
    int dz[MAX_MG_LEVELS];dz[0] = Rmg_G->get_PZ0_GRID(ct.FG_RATIO);
    int dimx = dx[0], dimy = dy[0], dimz = dz[0];
    int ixoff, iyoff, izoff, dx2, dy2, dz2;

    int maxlevel = ct.poi_parm.levels;
    for(int level=1;level <= ct.poi_parm.levels;level++)
    {
        dx2 = MG.MG_SIZE (dx[level-1], level-1, Rmg_G->get_NX_GRID(ct.FG_RATIO), Rmg_G->get_PX_OFFSET(ct.FG_RATIO), dimx, &ixoff, ct.boundaryflag);
        dy2 = MG.MG_SIZE (dy[level-1], level-1, Rmg_G->get_NY_GRID(ct.FG_RATIO), Rmg_G->get_PY_OFFSET(ct.FG_RATIO), dimy, &iyoff, ct.boundaryflag);
        dz2 = MG.MG_SIZE (dz[level-1], level-1, Rmg_G->get_NZ_GRID(ct.FG_RATIO), Rmg_G->get_PZ_OFFSET(ct.FG_RATIO), dimz, &izoff, ct.boundaryflag);

        if((dx2 < 2) || (dy2 < 2) || (dz2 < 2))
        {
            maxlevel = level - 1;
            break;
        }

        dx[level] = dx2;
        dy[level] = dy2;
        dz[level] = dz2;

    }

    MPI_Allreduce (MPI_IN_PLACE, &maxlevel, 1, MPI_INT, MPI_MIN, pct.grid_comm);
    ct.poi_parm.levels = maxlevel;

    // Write a copy of the options file
    if(pct.imgpe == 0) {

        // Write out options file
        std::string OptionsFile(ct.basename);
        OptionsFile = OptionsFile + ".options";

        FILE *fhand = fopen(OptionsFile.c_str(), "w");
        if (!fhand)
            throw RmgFatalException() <<  "Unable to write file in " << __FILE__ << " at line " << __LINE__ << "\n";

        for(auto it = ControlMap.begin();it != ControlMap.end(); ++it) {

            std::string KeyName = it->first;
            InputKey *ik = it->second;
            std::string KeyVal = ik->Print();
            fprintf(fhand, "%s = \"%s\"\n", KeyName.c_str(), KeyVal.c_str());fflush(NULL);
        }

        fclose(fhand);
    }

    // Set up exchange correlation type
    InputKey *ik = ControlMap["exchange_correlation_type"];
    if(*ik->Readintval == AUTO_XC) {
        // Type set from pp files
        Functional F( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
        F.set_dft_from_name_rmg(reordered_xc_type[ct.xctype]);
    }
    else {
        // Type set explicitly in input file
        std::string xc_type = reordered_xc_type[*ik->Readintval];
        Functional F( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
        F.set_dft_from_name_rmg(xc_type.c_str());
    }

}


