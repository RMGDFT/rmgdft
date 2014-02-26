
/* multigrid-parameter structure */
typedef struct
{

    /* number of global-grid pre/post smoothings and timestep */
    double gl_step;
    int gl_pre;
    int gl_pst;

    /* timestep for the subiteration */
    double sb_step;

    /* timestep for the Richardson-Iteration */
    double ri_step;

    /* lowest MG level */
    int levels;

    /* Number of Mu-cycles to use */
    int mucycles;

    /* Number of Smoother iterations on the coarsest level */
    int coarsest_steps;

} MG_PARM;

/** @name CONTROL
  @memo Main control structure
 
  This is a global structure declared as extern CONTROL ct
 
 */
typedef struct
{

    // Special variables when running with multiple images per node
    // The number of images stacked on a single node
    int images_per_node;

    // Image node index ranging from 0 to (images_per_node-1)
    int image_node_id;

    //pid_t main_thread_pid;
    int THREADS_PER_NODE;

    // requested mpi thread level
    int mpi_threadlevel;

    /** Description of the run. */
    char description[MAX_CHAR];

    // Precision level for mg_eig_state. Either single (4) or double (8) May be switched dynamically
    int mg_eig_precision;

    /* time at which run started */
    rmg_double_t time0;
    
    /* determine if this image is processing spin up or spin down. */
    int spin_flag;

    /* determine whether to initialize up and down density equally or not */
    int init_equal_density_flag; 

    /* determine whether to get_dos or not */
    int pdos_flag; 


    /** Name of the input control file. Passed as a command line argument
     *
     *  Example:
     *  bash$  md in.diamond8
     */
    char cfile[MAX_PATH];

    /* Basename for control file and other related files (beta and q function plots, etc.) */
    char basename[MAX_PATH];

    /** HAndle of the output log file. Constructed from command line argument */
    FILE *logfile;

    /** Input file name to read wavefunctions from when doing a restart */
    char infile[MAX_PATH];

    /** Input file name to write wavefunctions to */
    /* Output file name */
    char outfile[MAX_PATH];

    /** File to read the pseudopotentials from */
    /*  char pspfile[MAX_PATH]; */

    /** Initial run flag. Read from the input file. 0=initial run otherwise a restart */
    int runflag;

    /* output z-average of states */
    int zaverage;

    /* number of state to output */
    int plot_state;

    /* Exchage-Correlation flag */
    int xctype;

    /** Boundary condition flag. Read from the input file. 0=periodic, 1=cluster, 2=surface */
    int boundaryflag;

    /* Coordinate input flag: crystal or cartesian */
    int crd_flag;

    /* Maximum number of MD steps */
    int max_md_steps;

    /* Maximum number of rmg meta loops (NEB, ARTS, etc.) */
    int max_rmg_steps;

    /* MD steps iterator */
    int md_steps;

    /* Emin when get_dos */
    rmg_double_t Emin;

    /* Emax when get_dos */
    rmg_double_t Emax;

    /* number of energy points when get_dos */
    int E_POINTS;

    /* Maximum number of SCF steps in a MD steps */
    int max_scf_steps;

    /* Total number of SCF steps done */
    int total_scf_steps;

    /* SCF steps iterator */
    int scf_steps;


    /* convergence criterion */
    rmg_double_t thr_rms;

    /* force convergence criterion */
    rmg_double_t thr_frc;

    /* Number of steps after which to perform checkpointing */
    int checkpoint;

    /** Sorting flag for wavefunctions. Read from input file. 0=no sort, 1=sort */
    int sortflag;

    /** Number of states */
    int num_states;


    /*Number of states for spin up and down used for initialization*/
    int num_states_up, num_states_down;

    /** Number of unoccupied states above Fermi level */
    int num_unocc_states;

    /** string to store repeat count occupations */
    char occupation_str[256];

    /*string to store repeat count occupations for spin up*/
    char occupation_str_spin_up[256];

    /*string to store repeat count occupations for spin down*/
    char occupation_str_spin_down[256]; 
    

    /** Number of ions */
    int num_ions;

    /** Ion structure */
    ION *ions;

    /** Number of species */
    int num_species;

    /* Cutoff parameter */
    rmg_double_t cparm;
    rmg_double_t betacparm;
    rmg_double_t qcparm;

    /** Total conpensating charge density */
    rmg_double_t crho;

    /** Total charge in supercell */
    rmg_double_t tcharge;

    /** Norm conserving pseudo potential flag */
    int norm_conserving_pp;

    /** Species structure 
     * @see SPECIES */
    SPECIES *sp;

    /** the fine grid size on each coarse grid cell */
    int nxfgrid;
    int nyfgrid;
    int nzfgrid;

    /* Kohn-sham finite difference order */
    int kohn_sham_fd_order;

    /* This is the max of nldim for any species cubed */
    int max_nlpoints;
    int max_lpoints;
    int max_nlfpoints;
    int max_Qpoints;

    /** Maximum grid spacing in any coordinate direction */
    rmg_double_t hmaxgrid;


    /** Minimum grid spacing in any coordinate direction */
    rmg_double_t hmingrid;


    /** Physical grid basis size */
    int nbasis;

    /*Flag used to determine whether the actual ionic masses or equal masses set to 
     * that of carbon are to be used in fastrelax algorithm.*/
    int relax_mass;

    /** Density mixing parameter. Typical values range from 0.2 to 0.9, while larger values provide faster convergence as long as they are stable. */
    rmg_double_t mix;

    /*Order of Pulay mixing for charge density*/
    int charge_pulay_order;

    /*How often to refresh Pulay history*/
    int charge_pulay_refresh;

    /*Scale parameter for residuals in Pulay mixing*/
    rmg_double_t charge_pulay_scale;

    /*Flag to test whether or not the modified metrics should be used in Pulay mixing*/
    int charge_pulay_special_metrics;

    /*Weight for Pulay special metrics*/
    rmg_double_t charge_pulay_special_metrics_weight;

    /* Projector mixing parameter */
    rmg_double_t prjmix;

    /* Global uniform grid corner */
    rmg_double_t xcstart;
    rmg_double_t ycstart;
    rmg_double_t zcstart;


    /* Hartree potential offset from wavefunction grid */
    int vh_xoffset;
    int vh_yoffset;
    int vh_zoffset;


    /* Hartree potential grid sizes */
    int vh_nxgrid;
    int vh_nygrid;
    int vh_nzgrid;


    /* Hartree potential grid sizes per domain */
    int vh_pxgrid;
    int vh_pygrid;
    int vh_pzgrid;


    /* Total points in hartree potential per domain */
    int vh_pbasis;


    /* Wavefunction grid sizes */
    int psi_nxgrid;
    int psi_fnxgrid;
    int psi_nygrid;
    int psi_fnygrid;
    int psi_nzgrid;
    int psi_fnzgrid;

    /* Total points for wavefunctions */
    int psi_nbasis;
    int psi_fnbasis;

    /* Decoupled hartree potential */
    rmg_double_t *vh_ext;


    /* Mean min, and max wavefunction residuals for occupied space */
    rmg_double_t meanres;
    rmg_double_t minres;
    rmg_double_t maxres;

    /* total ionic charge */
    rmg_double_t ionic_charge;

    /* Variable occupation stuff */
    rmg_double_t nel;

    int occ_flag;

    rmg_double_t occ_width;

    rmg_double_t occ_mix;

    /** total background smearing charge -- for charged supercells */
    rmg_double_t background_charge;


    /** Multigrid parameters for the eigenvalue solver */
    MG_PARM eig_parm;

    /** Multigrid parameters for the poisson solver */
    MG_PARM poi_parm;


    /** Nose paramters */
    FINITE_T_PARM nose;

    /* force pointer array */
    int fpt[4];

    /* temperature control */
    int tcontrol;

    /* md integrate level */
    int mdorder;


    /* Milliken population flags. */
    int domilliken;
    int milliken;

    /* Diagonalization driver type */
    int subdiag_driver;

    /* Diagonalization flag and period */
    int initdiag;
    int diag;
    int end_diag;

    /* Folded spectrum width */
    rmg_double_t folded_spectrum_width;

    /* Flag indicating whether to use MPI_Allreduce operations or point to point in subdiag */
    int scalapack_global_sums;

    /* scalapack block size */
    int scalapack_block_factor;

    /* How many steps between writeout of eigenvalues*/
    int write_eigvals_period;

    /* Diagonalizations during the md step */
    int mddiag1;
    int mddiag2;

    /** Force flag. 0=don't compute forces, 1=compute forces */
    int forceflag;

    /*  relax_method: 0 = FIRE, fast relax, 1= LBFGS*/
    int relax_method;

    /* Should we process constraints on forces flag */
    int constrainforces;

    /* Whether to write full memory usage report at the end of calculation */
    int write_memory_report;

    /** Ionic motion timestep */
    rmg_double_t iondt;

    /*** Maximum ionic motion timestep */
    rmg_double_t iondt_max;

    /*** Factor by which iondt is increased */
    rmg_double_t iondt_inc;

    /*** Factor by which iondt is decreased */
    rmg_double_t iondt_dec;

    /*Number of steps after which iondt is increased */
    int relax_steps_delay;
    
    /*Number of steps since F * v was negative, used for dynamic timestep in structure optimization*/
    int relax_steps_counter;


    /** Ionic motion energy */
    rmg_double_t ionke;


    /* Total energies */
    rmg_double_t ES;
    rmg_double_t NUC;
    rmg_double_t KE;
    rmg_double_t XC;
    rmg_double_t NL;
    rmg_double_t II;
    rmg_double_t TOTAL;

    /* fermi energy */
    rmg_double_t efermi;

    /** Total number of k-points being used in the calculation */
    int num_kpts;


    /** K-point control structure */
    KPOINT *kp;

    /** The maximum number of projectors for any species */
    int max_nl;
    
    /** The maximum l quantum number for any species */
    int max_l;

    /*Maximum value of nldim for any species */
    int max_nldim;

    /*This keeps track whether ct.fftw_wisdom_setup was setup or not so that
     * we know whether to release wisdom memory at the end or not*/
    int fftw_wisdom_setup;


    /*Interpolation flags */
    int interp_flag;

    /*Order of B-spline */
    int interp_order;

    /*Order of trade_image used in interpolation */
    int interp_trade;

    /* the external electric field */
    rmg_double_t e_field;

    rmg_double_t x_field_0;

    rmg_double_t y_field_0;

    rmg_double_t z_field_0;

    rmg_double_t neb_spring_constant;

    /*Current RMS value*/
    rmg_double_t rms;

    /* Max number of sweeps in get_vh*/
    int hartree_max_sweeps;
    
    /* Min number of sweeps in get_vh*/
    int hartree_min_sweeps;

    /*Ratio between target RMS for get_vh and RMS total potential*/
    rmg_double_t hartree_rms_ratio;

    /*Boolean flag for using mask function filtering*/
    int mask_function;

    /* Potential acceleration constant step factor */
    rmg_double_t potential_acceleration_constant_step;

    /* Potential acceleration constant step factor */
    rmg_double_t potential_acceleration_poisson_step;

    // Some GPU information. Currently we use at most one device per MPI process
#if GPU_ENABLED

    // Support for gpu direct collectives
    int gpu_direct_collectives;

    // Integer device id
    int gpu_device_id;

    // Cuda device
    CUdevice  cu_dev;

    // Cuda device context
    CUcontext cu_context;

    // CUBLAS library handle
    cublasHandle_t cublas_handle;

    // cuda stream
    cudaStream_t cuda_stream;

    // GPU storage space for wavefunctions
    rmg_double_t *gpu_states;

    // GPU temporary storage space for wavefunctions
    rmg_double_t *gpu_temp;

    // GPU temporary storage space for weights
    rmg_double_t *gpu_weight;
    rmg_double_t *gpu_Bweight;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    rmg_double_t *gpu_work1;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    rmg_double_t *gpu_work2;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    rmg_double_t *gpu_work3;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    rmg_double_t *gpu_work4;

    // GPU storage space for matrix dimensioned (ct.num_states, ct.num_states)
    rmg_double_t *gpu_global_matrix;

    // Pinned host memory for finite difference routines. Allocation is slow so it
    // needs to be done once at initializatio time for each thread.
    rmg_double_t *gpu_host_temp1;
    rmg_double_t *gpu_host_temp2;
    rmg_double_t *gpu_host_temp3;
    rmg_double_t *gpu_host_temp4;
    rmg_double_t *gpu_host_fdbuf1;
    rmg_double_t *gpu_host_fdbuf2;
    rmg_double_t *gpu_host_work;

#endif    

    /* Compression level for some trade images routines */
    int trade_compression_level;

    


    
    /** File to read the pseudopotentials from */
    char pspfile[MAX_PATH];

    
    /*  spin polarize flag, 0 = no spin, 1 = spin polarized */
    int spin;

    
    int freeze_orbital_step;


    /* override occupations */
    int override_occ;

    /* Override of positions (during restart) */
    int override_atoms;

    /* Number of steps after which to output results */
    int outcount;



 
    /** Global uniform grid spacing in x */
    //rmg_double_t hxgrid;

    /** Global uniform grid spacing in y */
    //rmg_double_t hygrid;

    /** Global uniform grid spacing in z */
    //rmg_double_t hzgrid;

 

    /** bravais lattice type */
    //int ibrav;

    /* lengths of the sides of the supercell */
    //rmg_double_t xside;
    //rmg_double_t yside;
    //rmg_double_t zside;

 
    /** Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. 
      A value larger than 1.05 can lead to convergence problems. */
    //rmg_double_t anisotropy;


 


 
 
 

    /* Total points for potential */
    int vh_nbasis;


    rmg_double_t bg_begin;
    rmg_double_t bg_end;
    rmg_double_t BT;


   
    /* Total number of electrons */
    rmg_double_t num_el;

   
    /* movie flags */
    int rmvmovie, chmovie, xbsmovie;


 

    /* Desired vector for constrained dynamics */
    rmg_double_t cd_vector[3];

    /* number of velocity in constrained dynamics */
    rmg_double_t cd_velocity;



    rmg_double_t Evxcold_rho;
    rmg_double_t Evhold_rho;
    rmg_double_t Evh_rho;
    rmg_double_t Evh_rhoc;

    rmg_double_t TOTAL_former;
    rmg_double_t dE;

    rmg_double_t *energies;

    int restart_mix;

    int move_centers_at_this_step;



    int num_storage_st;


    /* Index telling if an orbital has been allocated */
    char allocate_st[MAX_STATES];

    /* Number of orbital in each global function */
    int nb_func[MAX_STATES];

    /* Initial storage of each orbital
       (to restart with same function) */
    int alloc_ini[MAX_STATES];

    /* "Center" of global function */
    rmg_double_t center[MAX_STATES][3];

    rmg_double_t Bc;
    rmg_double_t Bx;
    rmg_double_t By;
    rmg_double_t Bz;
    rmg_double_t Ac;
    rmg_double_t Ax;
    rmg_double_t Ay;
    rmg_double_t Az;
    rmg_double_t Axy;
    rmg_double_t Axz;
    rmg_double_t Ayz;

    int state_per_proc;
    int state_begin;
    int state_end;

    int max_orbit_size;
    int max_orbit_nx;
    int max_orbit_ny;
    int max_orbit_nz;

    int movingCenter;
    int movingSteps;
    int mg_method;
    int mg_steps;
    STATE *states;


    /* number of waves to plot in ON calculation */
    int num_waves;
   

    char file_atomic_orbit[MAX_SPECIES][MAX_PATH];

    /* output information for GW calculations.  --Xiaohong */
    int flag_gw;

    /*  for negf  */


    int metal;
    int num_blocks;
    int block_dim[MAX_BLOCKS];

    int num_cond_curve;
    int *cond_probe1;
    int *cond_probe2;

    int ion_begin;
    int ion_end;


 /* parameter to define the gate bias */
    double gbias_begin;
    double gbias_end;
    double gate_bias;

    /* parameter to define compensating potential */
    int vcomp_Lbegin;
    int vcomp_Lend;
    int vcomp_Rbegin;
    int vcomp_Rend;

    /* tag to determine whether to do auto 3Ddos printing for
 * transmission peaks */
    int auto_3Ddos;


    int plane[5];
    int simpson_depth;
    double simpson_tol;


} CONTROL;


/* Extern declaration for the main control structure */
extern CONTROL ct;


