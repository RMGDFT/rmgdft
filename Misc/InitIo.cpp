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
#include "GpuAlloc.h"
#include "Gpufuncs.h"

#if CUDA_ENABLED
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cublas_v2.h>
#endif

#if HIP_ENABLED
    #include <hip/hip_runtime_api.h>
    #include <hip/hip_vector_types.h>
    #include <rocfft.h>
    #include <rocsolver.h>
    #if MAGMA_LIBS
        #include <magma_v2.h>
    #endif
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
      //char cpunum[2];
      //sprintf(cpunum, "%d ", j);
      //strcat(str, cpunum);
    }
  }
  printf("Rank: %4d; NumRanks: %d; RankCore: %3d; Hostname: %s; Affinity: %3d\n", rank, nproc, core, hostname, count);
}
#endif



namespace Ri = RmgInput;


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

    snprintf (ct.cfile, sizeof(ct.cfile), "%s%s", pct.image_path[pct.thisimg], pct.image_input[pct.thisimg]);
    snprintf (ct.shortname, sizeof(ct.shortname), "%s%s", pct.image_path[pct.thisimg], pct.image_input[pct.thisimg]);

    ReadCommon(ct.cfile, ct, pct, ControlMap);

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

    /* Get the crystal or cartesian coordinates of the ions */
    init_pos ();

    ReadPseudo(ct.num_species, ct, ControlMap);

    if(ct.spinorbit)
    {
        bool pp_has_so = false;
        for(int isp = 0;isp < ct.num_species;isp++)
        {
            SPECIES *sp = &Species[isp];
            if(sp->is_spinorb) pp_has_so = true; 
        }
        if(!pp_has_so)
        {
            rmg_error_handler (__FILE__, __LINE__, "no pseudopotential has spin-orbit.\n");
        }
    }
    if(ct.noncoll) ct.is_ddd_non_diagonal = true;
    // If fine/coarse grid ratio is not set then autoset it. By default we
    // use 1 for norm conserving pseudopotentials and 2 for ultrasoft.

    ct.norm_conserving_pp = true;
    int nc_count = 0;
    int us_count = 0;
    for(int isp = 0;isp < ct.num_species;isp++)
    {
        SPECIES *sp = &Species[isp];
        if(sp->is_norm_conserving)
        {
            nc_count++;
        }
        else
        {
            us_count++;
            ct.norm_conserving_pp = false;
        }
    }
    if(nc_count && us_count)
    {
        rmg_error_handler (__FILE__, __LINE__, "Mixing norm conserving and ultrasoft pseudopotentials is not supported. Check your input files.\n");
    }

    if(!ct.FG_RATIO) ct.FG_RATIO = 2;

    // For USPP force a minimum of 2
    if(!ct.norm_conserving_pp) ct.FG_RATIO = std::max(2, ct.FG_RATIO);

    static Ri::ReadVector<int> WavefunctionGrid;
    double *celldm = Rmg_L.celldm;
    static double grid_spacing;

    InputKey *ik;

    ik = ControlMap["wavefunction_grid"];
    WavefunctionGrid = ik->Vint;

    int NX_GRID = WavefunctionGrid.vals.at(0);
    int NY_GRID = WavefunctionGrid.vals.at(1);
    int NZ_GRID = WavefunctionGrid.vals.at(2);

    if(NX_GRID * NY_GRID * NZ_GRID == 1)
    {
        ik = ControlMap["grid_spacing"];
        grid_spacing = ik->Readdoubleval[0];
        SetupWavefunctionGrid(NX_GRID, NY_GRID, NZ_GRID, celldm, grid_spacing);

        WavefunctionGrid.vals[0] = NX_GRID;
        WavefunctionGrid.vals[1] = NY_GRID;
        WavefunctionGrid.vals[2] = NZ_GRID;
        ik = ControlMap["wavefunction_grid"];
        ik->Vint = WavefunctionGrid;
        ik = ControlMap["grid_spacing"];
        ik->Readdoubleval = &grid_spacing;

    }

    /* Initialize symmetry stuff */

    ct.is_gamma = true;
    ct.is_gamma = ct.is_gamma && (ct.kpoint_mesh[0] == 1);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_mesh[1] == 1);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_mesh[2] == 1);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_is_shift[0] == 0);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_is_shift[1] == 0);
    ct.is_gamma = ct.is_gamma && (ct.kpoint_is_shift[2] == 0);
    ct.is_gamma = ct.is_gamma && (!ct.noncoll);
    ct.is_use_symmetry = ct.is_use_symmetry && (!ct.is_gamma);
    if(ct.AFM) ct.is_use_symmetry = 1;
    if(ct.is_use_symmetry) Rmg_Symm = new Symmetry(Rmg_L, NX_GRID, NY_GRID, NZ_GRID, ct.FG_RATIO);


    if(ct.forceflag == BAND_STRUCTURE)
    {
        if(ct.wannier90) {
            init_kpoints(ct.kpoint_mesh, ct.kpoint_is_shift);
        }
        else {
            int special_klines = ReadKpointsBandstructure(ct.cfile, ct, ControlMap);
            if(special_klines == 0) 
                ReadKpoints(ct.cfile, ct, ControlMap);
        }
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
    Rmg_halfgrid->set_rank(pct.gridpe, pct.grid_comm);


    if(Rmg_Symm) Rmg_Symm->setgrid(*Rmg_G, ct.FG_RATIO);


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
    if(!ct.coalesce_states) pct.coalesce_factor = 1;
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


    // Set up the Laplacian and gradient coefficients
    SetLaplacian();

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
        if (!ct.logfile)
            throw RmgFatalException() <<  "Unable to open logfile in " << __FILE__ << " at line " << __LINE__ << ". Possible disk full or permissions error?\n";
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


#if CUDA_ENABLED || HIP_ENABLED
    size_t deviceMem;
    int clock;
    char name[1024];

#if CUDA_ENABLED
    CUdevice cudev;
    if(pct.local_rank == 0)
    {
        gpuDeviceReset();
        cuInit(0);
        gpuSetDeviceFlags(gpuDeviceScheduleAuto);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(pct.local_rank != 0)
    {
        cuInit(0);
    }
    gpuGetDeviceCount( &ct.num_gpu_devices);

    // Get cuda version
    cudaError_t cuerr =  cudaDriverGetVersion ( &ct.cuda_version );
    rmg_printf ("\nCUDA version %d detected.\n", ct.cuda_version);
#endif

#if HIP_ENABLED
    hipDevice_t hipdev;
    if(pct.local_rank == 0)
    {
        gpuDeviceReset();
        hipInit(0);
        gpuSetDeviceFlags(gpuDeviceScheduleAuto);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(pct.local_rank != 0)
    {
        hipInit(0);
    }
    gpuGetDeviceCount( &ct.num_gpu_devices);

    // Get hip version
    hipError_t hiperr =  hipDriverGetVersion ( &ct.hip_version );
    if(hiperr == hipSuccess) {
        rmg_printf ("\nHIP version %d detected.\n", ct.hip_version);
    }
    else
    {
        rmg_printf ("\nHIP version NOT detected.\n");
    }
#endif

    // Get device list and memory capacities. While this is not ideal under all circumstances
    // we will find the device with the largest memory and use all devices that have just as
    // much memory
    std::vector<size_t> device_mem;
    rmg_printf("\n");
    for(int idevice = 0; idevice < ct.num_gpu_devices; idevice++ ) {
#if CUDA_ENABLED
        cuDeviceGet( &cudev, idevice );
        cuDeviceGetName( name, sizeof(name), cudev );
        cuDeviceTotalMem( &deviceMem, cudev );
        ct.gpu_mem[idevice] = deviceMem;
        cuDeviceGetAttribute( &clock, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, cudev );
        rmg_printf( "device %d: %s, %.1f MHz clock, %.1f MB memory\n", idevice, name, clock/1000.f, deviceMem/1024.f/1024.f );
        device_mem.push_back(deviceMem/1024.0/1024.0);
#endif

#if HIP_ENABLED
        hipDeviceGet( &hipdev, idevice );
        hipDeviceGetName( name, sizeof(name), hipdev );
        hipDeviceTotalMem( &deviceMem, hipdev );
        ct.gpu_mem[idevice] = deviceMem;
        hipDeviceGetAttribute( &clock, hipDeviceAttributeClockRate, hipdev );
        rmg_printf( "device %d: %s, %.1f MHz clock, %.1f MB memory\n", idevice, name, clock/1000.f, deviceMem/1024.f/1024.f );
        device_mem.push_back(deviceMem/1024.0/1024.0);
#endif

    }

    size_t max_mem = *std::max_element(device_mem.begin(), device_mem.end());
    ct.num_usable_gpu_devices = 0;
    int idevice = 0;
    int does_managed;
    ct.gpus_support_managed_memory = true;
    for(auto it = device_mem.begin();it != device_mem.end();++it) {
        if(*it >= 0.9*max_mem) {
            ct.gpu_device_ids[ct.num_usable_gpu_devices] = idevice;
#if CUDA_ENABLED
            cuDeviceGet( &ct.cu_devices[ct.num_usable_gpu_devices], idevice);
            cuDeviceGetAttribute( &does_managed, CU_DEVICE_ATTRIBUTE_MANAGED_MEMORY, ct.cu_devices[ct.num_usable_gpu_devices]);
#endif
#if HIP_ENABLED
            hipDeviceGet( &ct.hip_devices[ct.num_usable_gpu_devices], idevice);
            hipDeviceGetAttribute( &does_managed, hipDeviceAttributeManagedMemory, ct.hip_devices[ct.num_usable_gpu_devices]);
#endif
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
        gpuSetDevice(ct.gpu_device_ids[pct.local_rank]);
#if CUDA_ENABLED
        if( CUBLAS_STATUS_SUCCESS != cublasCreate(&ct.cublas_handle) ) {
            rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
        }
        if( CUBLAS_STATUS_SUCCESS != cublasXtCreate(&ct.cublasxt_handle) ) {
            rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
        }
        cublasXtDeviceSelect(ct.cublasxt_handle, 1, &ct.gpu_device_ids[pct.local_rank]);
        cublasXtSetBlockDim(ct.cublasxt_handle, 2048);
        ct.gpublas_handle = ct.cublas_handle;
#endif
#if HIP_ENABLED
        ct.hip_dev = ct.gpu_device_ids[pct.local_rank];
        hipDeviceReset();
        hipSetDeviceFlags(hipDeviceScheduleAuto);
        if( HIPBLAS_STATUS_SUCCESS != hipblasCreate(&ct.hipblas_handle) ) {
            rmg_error_handler (__FILE__, __LINE__, "HIPBLAS: Handle not created\n");
        }
        ct.gpublas_handle = ct.hipblas_handle;
#endif

    }
    else if(pct.procs_per_host >= ct.num_usable_gpu_devices)
    {
        // Round robin
        int next_gpu = 0;
        for(int rank = 0;rank < pct.procs_per_host;rank++)
        {
            if(rank == pct.local_rank)
            {
                gpuSetDevice(ct.gpu_device_ids[next_gpu]);
#if CUDA_ENABLED
                if( CUBLAS_STATUS_SUCCESS != cublasCreate(&ct.cublas_handle) ) {
                    rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
                }
                if( CUBLAS_STATUS_SUCCESS != cublasXtCreate(&ct.cublasxt_handle) ) {
                    rmg_error_handler (__FILE__, __LINE__, "CUBLAS: Handle not created\n");
                }
                cublasXtDeviceSelect(ct.cublasxt_handle, 1, &ct.gpu_device_ids[next_gpu]);
                cublasXtSetBlockDim(ct.cublasxt_handle, 2048);
                ct.gpublas_handle = ct.cublas_handle;
#endif
#if HIP_ENABLED
                ct.hip_dev = ct.gpu_device_ids[next_gpu];
                hipDeviceReset();
                hipSetDeviceFlags(hipDeviceScheduleAuto);
                if( HIPBLAS_STATUS_SUCCESS != hipblasCreate(&ct.hipblas_handle) ) {
                    rmg_error_handler (__FILE__, __LINE__, "HIPBLAS: Handle not created\n");
                }
                ct.gpublas_handle = ct.hipblas_handle;
#endif
            }
            next_gpu++;
            next_gpu = next_gpu % ct.num_usable_gpu_devices;
        }

    }
    else
    {
#if CUDA_ENABLED
        if( CUDA_SUCCESS != cuDeviceGet( &ct.cu_dev, 0 ) ) {
            rmg_error_handler (__FILE__, __LINE__, "CUDA: Cannot get the device\n");
        }
        gpuSetDevice(ct.cu_dev);
#endif
#if HIP_ENABLED
        // Still need to figure out a way to handle this case
        gpuSetDevice(ct.hip_dev);
        hipDeviceReset();
        hipSetDeviceFlags(hipDeviceScheduleAuto);
        if( HIPBLAS_STATUS_SUCCESS != hipblasCreate(&ct.hipblas_handle) ) {
            rmg_error_handler (__FILE__, __LINE__, "HIPBLAS: Handle not created\n");
        }
        ct.gpublas_handle = ct.hipblas_handle;
#endif
    }

#if CUDA_ENABLED
    cusolverStatus_t cusolver_status = cusolverDnCreate(&ct.cusolver_handle);
    if(cusolver_status != CUSOLVER_STATUS_SUCCESS)
    {
        rmg_error_handler (__FILE__, __LINE__, "cusolver initialization failed.\n");
    }
    gpuStreamCreateWithFlags(&ct.cusolver_stream, gpuStreamNonBlocking);
    cusolver_status = cusolverDnSetStream(ct.cusolver_handle, ct.cusolver_stream);
    if(cusolver_status != CUSOLVER_STATUS_SUCCESS)
    {
        rmg_error_handler (__FILE__, __LINE__, "cusolver stream initialization failed.\n");
    }
#endif

#if HIP_ENABLED
    rocfft_setup();
    rocsolver_create_handle(&ct.roc_handle);
#if MAGMA_LIBS
    magma_init();
    magma_setdevice(ct.hip_dev);
#endif
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
        int max_threads = std::max(ct.MG_THREADS_PER_NODE, ct.OMP_THREADS_PER_NODE);
#if USE_NUMA && USE_HWLOC
        Rmg_Q = new MpiQueue(64, max_threads, pct.manager_cpumask, &pct.topology, ct.spin_manager_thread, ct.spin_worker_threads);
#elif USE_NUMA
        Rmg_Q = new MpiQueue(64, max_threads, pct.manager_cpumask, NULL, ct.spin_manager_thread, ct.spin_worker_threads);
#elif USE_HWLOC
        Rmg_Q = new MpiQueue(64, max_threads, NULL, &pct.topology, ct.spin_manager_thread, ct.spin_worker_threads);
#else
        Rmg_Q = new MpiQueue(64, max_threads, NULL, NULL, ct.spin_manager_thread, ct.spin_worker_threads);
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

    // Now we have to determine the max multigrid level where the Kohn-Sham solver can still
    // use the offset routines. If 
    dx[0] = pct.coalesce_factor * Rmg_G->get_PX0_GRID(1);
    dy[0] = Rmg_G->get_PY0_GRID(1);
    dz[0] = Rmg_G->get_PZ0_GRID(1);
    maxlevel = ct.eig_parm.levels;
    int minsize, maxsize;
#if 0
    ct.mg_offset_level = 0;
    for(int level=1;level <= ct.eig_parm.levels;level++)
    {
        dx2 = MG.MG_SIZE (dx[level-1], level-1, Rmg_G->get_NX_GRID(1), Rmg_G->get_PX_OFFSET(1), dimx, &ixoff, ct.boundaryflag);
        dy2 = MG.MG_SIZE (dy[level-1], level-1, Rmg_G->get_NY_GRID(1), Rmg_G->get_PY_OFFSET(1), dimy, &iyoff, ct.boundaryflag);
        dz2 = MG.MG_SIZE (dz[level-1], level-1, Rmg_G->get_NZ_GRID(1), Rmg_G->get_PZ_OFFSET(1), dimz, &izoff, ct.boundaryflag);

        if((dx2 < 2) || (dy2 < 2) || (dz2 < 2))
        {
            maxlevel = level - 1;
            break;
        }
        dx[level] = dx2;
        dy[level] = dy2;
        dz[level] = dz2;

        minsize=dx2, maxsize=dx2;
        MPI_Allreduce(MPI_IN_PLACE, &minsize, 1, MPI_INT, MPI_MIN, pct.grid_comm);
        MPI_Allreduce(MPI_IN_PLACE, &maxsize, 1, MPI_INT, MPI_MAX, pct.grid_comm);
        if(minsize != maxsize && minsize < 4) break;
        minsize=dy2, maxsize=dy2;
        MPI_Allreduce(MPI_IN_PLACE, &minsize, 1, MPI_INT, MPI_MIN, pct.grid_comm);
        MPI_Allreduce(MPI_IN_PLACE, &maxsize, 1, MPI_INT, MPI_MAX, pct.grid_comm);
        if(minsize != maxsize && minsize < 4) break;
        minsize=dz2, maxsize=dz2;
        MPI_Allreduce(MPI_IN_PLACE, &minsize, 1, MPI_INT, MPI_MIN, pct.grid_comm);
        MPI_Allreduce(MPI_IN_PLACE, &maxsize, 1, MPI_INT, MPI_MAX, pct.grid_comm);
        if(minsize != maxsize && minsize < 4) break;
        ct.mg_offset_level++;
    }
#endif

    // Write a copy of the options file
    if((pct.imgpe == 0)  && (ct.rmg_branch != RMG_TO_QMC)) {

        // Write out options file
        std::string OptionsFile(ct.basename);
        OptionsFile = OptionsFile + ".options";

        FILE *fhand = fopen(OptionsFile.c_str(), "w");
        if (!fhand)
            throw RmgFatalException() <<  "Unable to write file in " << __FILE__ << " at line " << __LINE__ << "\n";

        std::map<InputKey *, InputKey *, keycompare> SortedMap;
        for(auto it = ControlMap.begin();it != ControlMap.end(); ++it)
        {
            std::string KeyName = it->first;
            InputKey *ik = it->second;
            SortedMap.insert(std::make_pair(ik, ik));
        }

        for(auto it = SortedMap.begin();it != SortedMap.end(); ++it) {

            InputKey *ik = it->first;
            std::string KeyName = ik->KeyName;
            std::string KeyVal = ik->Print();
            fprintf(fhand, "%s = \"%s\"\n", KeyName.c_str(), KeyVal.c_str());
        }
        fclose(fhand);
    }

    // Set up exchange correlation type
    Functional F( *Rmg_G, Rmg_L, *Rmg_T, ct.is_gamma);
    ik = ControlMap["exchange_correlation_type"];
    if(*ik->Readintval == AUTO_XC) {
        // Type set from pp files
        F.set_dft_from_name_rmg(reordered_xc_type[ct.xctype]);
    }
    else {
        // Type set explicitly in input file
        std::string xc_type = reordered_xc_type[*ik->Readintval];
        F.set_dft_from_name_rmg(xc_type);
    }
    ct.xc_is_hybrid = F.dft_is_hybrid_rmg();
    if(ct.exx_fraction < 0.0)
        ct.exx_fraction = F.get_exx_fraction_rmg();
    else
        F.set_exx_fraction_rmg(ct.exx_fraction);


#if HIP_ENABLED || CUDA_ENABLED
    size_t factor = 2;
    if(ct.is_gamma) factor = 1;
    int images = ct.kohn_sham_fd_order / 2;
    size_t bufsize = factor * pct.coalesce_factor *
                     (Rmg_G->get_PX0_GRID(1) + 2*images) * (Rmg_G->get_PY0_GRID(1) + 2*images) * (Rmg_G->get_PZ0_GRID(1) + 2*images)*sizeof(double);
#if HIP_ENABLED
    init_hip_fd(ct.MG_THREADS_PER_NODE, bufsize);
#endif
#if CUDA_ENABLED
    init_cuda_fd(ct.MG_THREADS_PER_NODE, bufsize);
#endif

#endif
}


