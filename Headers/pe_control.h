#ifndef RMG_pe_control_H
#define RMG_pe_control_H 1

#ifdef USE_NUMA
    #include <numa.h>
#endif
#ifdef USE_HWLOC
    #include <hwloc.h>
#endif
#include "Mgrid.h"
#include "Scalapack.h"


/** @name PE_CONTROL
  *
  * @memo Processor control structure.
  * 
  * This is a global structure declared as extern PE_CONTROL pct.
  * 
 */
typedef struct
{

    /** Number (rank in MPI terminology) of this processor in this image grid */
    int gridpe, imgpe, thisimg, spinpe;

    /** Number of grids (typically 1) per image to be run simultaneously **/
    int images, grids;
    char image_path[MAX_IMGS][MAX_PATH];
    char image_input[MAX_IMGS][MAX_PATH];
    int image_npes[MAX_IMGS];
    int worldrank;
    int total_npes;
    int grid_npes;

    /* MPI communicators for each code grid (grid_comm) and one (rmg_comm)
     * for all group rank 0 pe's. The later effectively replaces MPI_COMM_WORLD
     * unless you really need all-to-all, even across grids, communication. */
    MPI_Comm rmg_comm, img_topo_comm, grid_topo_comm, grid_comm, img_comm, spin_comm, scalapack_comm;
    MPI_Comm kpsub_comm, allkp_comm, my_comm;

    MPI_Comm img_cross_comm;
    // Coalesced grid MPI communicator
    MPI_Comm coalesced_grid_comm;

    // Coalesced local MPI communicator
    MPI_Comm coalesced_local_comm;

    // coalesce factor
    int coalesce_factor;

    // Number of physical hosts
    int total_hosts;

    // Number of MPI procs per physical host
    int procs_per_host;

    // Number of cpu cores per physical host
    int ncpus;

    // rank for image communication in NEB relax
    int left_img_rank, right_img_rank;

#if defined(USE_NUMA) || defined(USE_HWLOC)
    // Numa nodes per host
    int numa_nodes_per_host;

    // cpumask for mpi process
    struct bitmask *cpumask;

    // nodemask for mpi process
    struct bitmask *nodemask;

    // cpumask for manager thread
    struct bitmask *manager_cpumask;

#endif

#ifdef USE_HWLOC
    hwloc_topology_t topology;
    hwloc_const_bitmap_t pu_set;
    hwloc_const_bitmap_t nn_set;
#endif

    // Local rank of this proc
    int local_rank;

    // MPI rank of each proc local to this host (dimension procs_per_host)
    int *mpi_local_ranks;

    // Communicator associated with the local ranks.
    MPI_Comm local_comm;

    // local master
    int is_local_master;

    // Communicator associated with the master ranks on each node
    MPI_Comm local_master_comm;

    /* scalapack variables */
    int desca[DLEN];
    int ictxt;

    /*Whether pe participates in scalapack calculations*/
    int scalapack_pe;
    int scalapack_npes;

    /* Row or Column that given processor handles when using scalapack*/
    int scalapack_myrow;
    int scalapack_mycol;

    /*Processor distribution for scalapack*/
    int scalapack_nprow;
    int scalapack_npcol;

    /* desca for all nodes that participate */
    int *scalapack_desca;

    /* Max dist matrix size for any scalapack PE */
    int scalapack_max_dist_size;

    /** Processor x-coordinate for domain decomposition */
    int pe_x;
    /** Processor y-coordinate for domain decomposition */
    int pe_y;
    /** Processor z-coordinate for domain decomposition */
    int pe_z;

    /** Number of points in the circle of local projector for each pocessor*/
    int *lptrlen;

    /** points to start of DnmI function storage for this ion*/
//    double **dnmI;
//    double **dnmI_x;
//    double **dnmI_y;
//    double **dnmI_z;

    /** points to start of qqq storage for this ion*/
//    double **qqq;

    int num_owned_ions;
    int num_loc_ions;
    int *loc_ions_list;

    double *localpp;
    double *localrhoc;
    double *localrhonlcc;
    double *localatomicrho;

    int instances;
    /** Neighboring processors in three-dimensional space */

    /** Processor kpoint- coordinate for domain decomposition */
    /*  paralleled for kpoint */
    int pe_kpoint;

    /* kpoint index for start and end for a subdomain processors */
    int kstart;
    int kend;

    /*  processor coordinates in COMM_KP communicator */
    int coords[2];

    /* Number of ions centered on this processor */
    int n_ion_center;
    int n_ion_center_loc;

    /* Projectors per ion in a given region */
    int *prj_per_ion;

    /* Indices of the ions within non-local range */
    int *ionidx;
    int *ionidx_loc;

    /* The size of local orbitals on this PE */
    int psi_size;

    /* pointer to former step solution, used in pulay and KAIN mixing  */
    int descb[DLEN];

    double *psi1, *psi2;
    int num_local_orbit;

} PE_CONTROL;


/* Extern declaration for the processor control structure */
#if __cplusplus
extern "C" PE_CONTROL pct;
#else
extern PE_CONTROL pct;
#endif


#endif
