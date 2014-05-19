
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

	/* MPI communicators for each code grid (grid_comm) and one (rmg_comm)
	 * for all group rank 0 pe's. The later effectively replaces MPI_COMM_WORLD
	 * unless you really need all-to-all, even across grids, communication. */
	MPI_Comm rmg_comm, img_topo_comm, grid_topo_comm, grid_comm, img_comm, spin_comm, scalapack_comm;



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

    /* scalapack mpi rank for all nodes that participate */
    int scalapack_mpi_rank[MAX_PES];

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


    /** Points to start of projector storage for this ion in projector space */
    rmg_double_t *weight;
    rmg_double_t *Bweight;

#if FDIFF_BETA
    /*These are used for non-local force */
    rmg_double_t **weight_derx;
    rmg_double_t **weight_dery;
    rmg_double_t **weight_derz;
#endif


    /** An index array which maps the projectors onto the 3-d grid associated
        with each processor.
    */
    int **nlindex;
    int **Qindex;

    /** An index array which indicate whether the grid map on the current pocessor*/
    int **idxflag;
    int **Qdvec;

    /** Number of points in the nlindex array for each ion */
    int *idxptrlen;
    int *Qidxptrlen;

    /** Number of points in the circle of local projector for each pocessor*/
    int *lptrlen;

    /** Phase shifts for the non-local operators */
    rmg_double_t **phaseptr;

    /** Points to start of storage for theaugument function*/
    rmg_double_t **augfunc;

    /** points to start of DnmI function storage for this ion*/
    rmg_double_t **dnmI;
    rmg_double_t **dnmI_x;
    rmg_double_t **dnmI_y;
    rmg_double_t **dnmI_z;

    /** points to start of qqq storage for this ion*/
    rmg_double_t **qqq;


    int num_owned_ions;
    int owned_ions_list[MAX_NONLOC_IONS];
    
    int num_nonloc_ions;
    int nonloc_ions_list[MAX_NONLOC_IONS];
    int nonloc_ion_ownflag[MAX_NONLOC_IONS];

    int num_nonloc_pes;
    int nonloc_pe_list[MAX_NONLOC_PROCS];
    int nonloc_pe_list_ions[MAX_NONLOC_PROCS][MAX_NONLOC_IONS];
    int nonloc_pe_num_ions[MAX_NONLOC_PROCS];
    
    
    /*For ions owned by current PE*/
    /* Number of cores to cummunicate with about owned ions*/
    int num_owned_pe;
    /*List of ranks of cores to comunicate with about owned ions*/
    int owned_pe_list[MAX_NONLOC_PROCS];
    /* Number of owned ions to communicate about for cores from owned_pe_list */
    int num_owned_ions_per_pe[MAX_NONLOC_PROCS];
    /*List of ion indices to communicate about for core from owned_pe_list 
     * These are indices within nonloc ions, not absolute ones*/ 
    int list_owned_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS];
    
    
    /*For ions NOT owned by current PE*/
    /* Number of cores to cummunicate with about non-owned ions*/
    int num_owners;
    /*Indices of cores to cummunicate with about non-owned ions*/
    int owners_list[MAX_NONLOC_PROCS];
    /* Number of non-owned ions to communicate about for cores from owners_list  */
    int num_nonowned_ions_per_pe[MAX_NONLOC_PROCS];
    /*List of ion indices to communicate about for cores from owners_list
     * These are indices within nonloc ions, not absolute ones*/ 
    int list_ions_per_owner[MAX_NONLOC_PROCS][MAX_NONLOC_IONS];
    
    rmg_double_t *oldsintR_local;
    rmg_double_t *oldsintI_local;
    rmg_double_t *newsintR_local;
    rmg_double_t *newsintI_local;

    // Holds non-local and S operators acting on psi
    rmg_double_t *nv;
    rmg_double_t *ns;
    rmg_double_t *Bns;
    int num_tot_proj;
    rmg_double_t *M_dnm;
    rmg_double_t *M_qqq;


 
    /* Grid sizes on each PE */
//    int PX0_GRID;
//    int PY0_GRID;
//    int PZ0_GRID;

    /* Grid offsets on each PE */
//    int PX_OFFSET;
//    int PY_OFFSET;
//    int PZ_OFFSET;

    /* Basis size on each PE */
//    int P0_BASIS;

    /* Fine grid sizes on each PE */
//    int FPX0_GRID;
//    int FPY0_GRID;
//    int FPZ0_GRID;

    /* Fine Grid offsets on each PE */
//    int FPX_OFFSET;
//    int FPY_OFFSET;
//    int FPZ_OFFSET;

    /* Fine grid basis size on each PE */
//    int FP0_BASIS;

    int instances;
    /** Neighboring processors in three-dimensional space */
    int neighbors[6];

    /** Processor kpoint- coordinate for domain decomposition */
    /*  paralleled for kpoint */
    int pe_kpoint;

 
    /** local grid size for x,y,z **/
//    int nx_grid;
//    int ny_grid;
//    int nz_grid;

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

    /* Points to start of projectors for this ion in projector space */
    /* All projectors are stored consecutively.                      */
    int *prj_ptr;

    /* Pointers to the index array for the non-local potentials */
    int *idxptr;




    /* The size of local orbitals on this PE */
    int psi_size;
    /* pointer to former step solution, used in pulay and KAIN mixing  */
    int descb[DLEN];
    int mycol;
    int myrow;
    int nprow;
    int npcol;

    int num_local_orbit;
    rmg_double_t *psi1, *psi2;

} PE_CONTROL;


/* Extern declaration for the processor control structure */
extern PE_CONTROL pct;



