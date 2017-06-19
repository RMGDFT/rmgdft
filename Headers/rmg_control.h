#ifndef RMG_control_H
#define RMG_control_H 1

#include <stdbool.h>
#ifdef __cplusplus
    #include <complex>
#else
    #include <complex.h>
#endif


/* multigrid-parameter structure */
typedef struct
{

    /* number of global-grid pre/post smoothings and timestep */
    double gl_step;
    int gl_pre;
    int gl_pst;

    /* timestep for multigrid correction */
    double mg_timestep;

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

    // Energy output units 0=Hartrees, 1=Rydbergs
    int energy_output_units;
    double energy_output_conversion[2];
    char *energy_output_string[2];

    // Discretization type flag
    int discretization;

    // Special variables when running with multiple images per node
    // The number of images stacked on a single node
    int images_per_node;

    // Image node index ranging from 0 to (images_per_node-1)
    int image_node_id;

    // OpenMP threads per node
    int THREADS_PER_NODE;

    // Multigrid threads per node. Normally the same as above but may be tunred for performance.
    int MG_THREADS_PER_NODE;

    // requested mpi thread level
    int mpi_threadlevel;

    // mpi queued mode
    bool mpi_queue_mode;

    /** Description of the run. */
    char description[MAX_CHAR];

    /* time at which run started */
    double time0;
    
    /* Flag indicating if this is a spin polarized calculation */
    bool spin_polarization;

    /* determine if this image is processing spin up or spin down. */
    bool spin_flag;

    /* determine whether to initialize up and down density equally or not */
    bool init_equal_density_flag; 

    /* determine whether to get_dos or not */
    bool pdos_flag; 


    /** Name of the input control file. Passed as a command line argument
     *
     *  Example:
     *  bash$  md in.diamond8
     */
    char cfile[2*MAX_PATH];

    /* Part of the log file name without numerical increment and ".log" added */
    char shortname[MAX_PATH];

    /* Basename for control file and other related files (beta and q function plots, etc.)  This is full name except the ".log" extension*/
    char basename[2*MAX_PATH];

    /** Handle of the output log file. Constructed from command line argument */
    FILE *logfile;

    /** Actual full name to the log file */
    char logname[2*MAX_PATH];

    /** Input file name to read wavefunctions from when doing a restart */
    char infile[2*MAX_PATH];

    /** Input file name to write wavefunctions to */
    /* Output file name */
    char outfile[2*MAX_PATH];

    char infile_tddft[2*MAX_PATH];
    char outfile_tddft[2*MAX_PATH];
    bool restart_tddft;

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
    double Emin;

    /* Emax when get_dos */
    double Emax;

    /* number of energy points when get_dos */
    int E_POINTS;

    /* Maximum number of SCF steps in a MD steps */
    int max_scf_steps;
    int tddft_steps;

    /* Total number of SCF steps done */
    int total_scf_steps;

    /* SCF steps iterator */
    int scf_steps;

    /* RMS[dV] convergence criterion */
    double thr_rms;

    /* estimated total energy  convergence criterion */
    double thr_energy;

    /* preconditioner single/double precision switch threshold */
    double preconditioner_thr;

    /* force convergence criterion */
    double thr_frc;

    /* GW convergence criterion */
    double gw_threshold;

    /* Fraction of full spectrum that gw_threshold is applied to */
    double gw_residual_fraction;

    /* Number of steps after which to perform checkpointing */
    int checkpoint;

    /** Sorting flag for wavefunctions. Read from input file. 0=no sort, 1=sort */
    bool sortflag;
    
    /** Flag for verbose output. Read from input file. 0=short output, 1=verbose */
    bool verbose;

    /** Number of states. May switch between init_states and run_states */
    int num_states;

    /** Number of initialization states */
    int init_states;

    /** Extra random states used for LCAO start above and beyond the atomic orbitals. */
    int extra_random_lcao_states;

    /** Number of run states */
    int run_states;

    /** Maximum number of states, used for davidson */
    int max_states;

    /** Number of states to allocate memory for */
    int alloc_states;
    int state_block_size;

    /** total number of atomic orbitals including the m dependency */
    int total_atomic_orbitals;

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
    
    /* Coalesce states */
    bool coalesce_states;

    /* Renormalize forcecs */
    bool renormalize_forces;

    /** Number of ions */
    int num_ions;
    
    /** Number of TF ions */
    /** Used for simple solvent model */
    int num_tfions;

    /** Ion structure */
    ION *ions;
    
    /** TF Ion structure */
    /** Used for simple solvent model */
    TF_ION *tf_ions;

    /** Number of species */
    int num_species;

    // potential grid refinement
    int FG_RATIO;

    /* Cutoff parameter */
    double cparm;        // Multiplicative factor
    double betacparm;    // For beta functions (non-local projectors)
    double rhocparm;     // For qfunctions and density

    /** Total conpensating charge density */
    double crho;

    /** Total charge in supercell */
    double tcharge;

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
    int force_grad_order;

    /* This is the max of nldim for any species cubed */
    int max_nlpoints;
    int max_lpoints;
    int max_nlfpoints;
    int max_Qpoints;

    /** Maximum grid spacing in any coordinate direction */
    double hmaxgrid;


    /** Minimum grid spacing in any coordinate direction */
    double hmingrid;


    /** Physical grid basis size */
    int nbasis;

    /*Flag used to determine whether the actual ionic masses or equal masses set to 
     * that of carbon are to be used in fastrelax algorithm.*/
    int relax_mass;

    /** Density mixing parameter. Typical values range from 0.2 to 0.9, while larger values provide faster convergence as long as they are stable. */
    double mix;

    /*Order of Pulay mixing for charge density*/
    int charge_pulay_order;

    /*How often to refresh Pulay history*/
    int charge_pulay_refresh;

    /*Scale parameter for residuals in Pulay mixing*/
    double charge_pulay_scale;

    /*Flag to test whether or not the modified metrics should be used in Pulay mixing*/
    bool charge_pulay_special_metrics;

    /*Weight for Pulay special metrics*/
    double charge_pulay_special_metrics_weight;

    /*Order of Broyden mixing for charge density*/
    int charge_broyden_order;

    /* Scale factor for Broyden mixing of charge density*/
    double charge_broyden_scale;

    /* For Davidson diagonalization. Factor to multipy desired tolerance for unoccupied states by. */
    double unoccupied_tol_factor;

    /* Type of charge analysis*/
    int charge_analysis_type;
    
    /* How often (in terms of SCF steps) to perform and write out charge analysis*/
    int charge_analysis_period;

    /* Projector mixing parameter */
    double prjmix;

    /* Global uniform grid corner */
    double xcstart;
    double ycstart;
    double zcstart;


    /* Hartree potential grid sizes per domain */
    int vh_pxgrid;
    int vh_pygrid;
    int vh_pzgrid;


    /* Total points in hartree potential per domain */
    int vh_pbasis;


    /* Total points for wavefunctions */
    int psi_nbasis;
    int psi_fnbasis;

    /* Decoupled hartree potential */
    double *vh_ext;


    /* Mean min, and max wavefunction residuals for occupied space */
    double meanres;
    double minres;
    double maxres;

    /* total ionic charge */
    double ionic_charge;

    /* Variable occupation stuff */
    double nel;

    /* spin up and down electronic charge */
    double nel_up;
    double nel_down;

    int occ_flag;

    double occ_width;

    double occ_mix;
    int mp_order;

    /** total background smearing charge -- for charged supercells */
    double background_charge;


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

    /* Kohn sham solver type */
    int kohn_sham_solver;

    /*Charge Mixing type*/
    int charge_mixing_type;

    /* Diagonalization flag and period */
    bool initdiag;
    int diag;
    int end_diag;

    /* Use folded spectrum if possible */
    bool use_folded_spectrum;

    /* Folded spectrum width */
    double folded_spectrum_width;

    /* Expansion factor for non-local projectors */
    double projector_expansion_factor;

    /* Flag indicating whether to use MPI_Allreduce operations or point to point in subdiag */
    bool scalapack_global_sums;

    /* scalapack block size */
    int scalapack_block_factor;

    /* cublasxt block size */
    int cublasxt_block_size;

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

    /** Ionic motion timestep */
    double iondt;

    /*** Maximum ionic motion timestep */
    double iondt_max;

    /*** Factor by which iondt is increased */
    double iondt_inc;

    /*** Factor by which iondt is decreased */
    double iondt_dec;

    /*Number of steps after which iondt is increased */
    int relax_steps_delay;
    
    /*Number of steps since F * v was negative, used for dynamic timestep in structure optimization*/
    int relax_steps_counter;


    /** Ionic motion energy */
    double ionke;


    /* Total energies */
    double ES;
    double NUC;
    double KE;
    double XC;
    double NL;
    double II;
    double TOTAL;

    /* fermi energy */
    double efermi;

    /** Total number of k-points being used in the calculation */
    int num_kpts;
    int num_kpts_pe;



    /** K-point control structure */
    KPOINT *kp;
    int kpoint_mesh[3];
    int kpoint_is_shift[3];

    /** The maximum number of projectors for any species */
    int max_nl;
    
    /** The maximum l quantum number for any species */
    int max_l;

    /*Maximum value of nldim for any species */
    int max_nldim;

    /* Whether or not Bweight is required (not needed for Central FD operator with NCPP */
    bool need_Bweight;

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
    double e_field;

    double x_field_0;

    double y_field_0;

    double z_field_0;

    double neb_spring_constant;

    /*Current RMS value*/
    double rms;

    /* 2nd order convergence term  (Vh_out - Vh_in)*(rho_out - rho_in) */
    double scf_accuracy;

    /* variational correction term for potential convergence */
    double scf_correction;

    /* Max number of sweeps in get_vh*/
    int hartree_max_sweeps;
    
    /* Min number of sweeps in get_vh*/
    int hartree_min_sweeps;

    /*Ratio between target RMS for get_vh and RMS total potential*/
    double hartree_rms_ratio;

    /* Potential acceleration constant step factor */
    double potential_acceleration_constant_step;

    /* Potential acceleration constant step factor */
    double potential_acceleration_poisson_step;

    // Some GPU information. Currently we use at most one device per MPI process
#if GPU_ENABLED

    // Total number of gpu devices present in the node
    int num_gpu_devices;

    // Total number of usable gpu devices present in the node
    int num_usable_gpu_devices;

    // Device ids for the usable gpus
    int gpu_device_ids[MAX_GPU_DEVICES];

    // GPU memory for the usable devices
    size_t gpu_mem[MAX_GPU_DEVICES];

    // Cuda device
    CUdevice  cu_dev;

    // Cuda device context
    CUcontext cu_context;

    // CUBLAS library handles
    cublasHandle_t cublas_handle;
    cublasXtHandle_t cublasXt_handle;

    // cuda stream
    cudaStream_t cuda_stream;

    // GPU storage space for wavefunctions
    double *gpu_states;

    // GPU temporary storage space for wavefunctions

    // GPU temporary storage space for weights

    // Pinned host memory for finite difference routines. Allocation is slow so it
    // needs to be done once at initializatio time for each thread.
    double *gpu_host_work;

    cuDoubleComplex *gpu_Htri, *gpu_Gtri, *gpu_Grow;
    cuDoubleComplex *gpu_GdiagBlocks;
    cuDoubleComplex *gpu_Imatrix, *gpu_Hii,  *gpu_temp, *gpu_Gii;
    cuDoubleComplex *gpu_Gcol;


    int *gpu_ipiv;

#else
    double *gpu_Htri, *gpu_Gtri, *gpu_GdiagBlocks;
    double *gpu_Grow;
    double *gpu_Gcol;
    double *gpu_Imatrix, *gpu_Hii,  *gpu_temp, *gpu_Gii;
#endif

    
   /* RMG2BGW options */
   bool rmg2bgw;
    double ecutrho, ecutwfc;

   
   int vxc_diag_nmin;
   int vxc_diag_nmax;

    
    /** File to read the pseudopotentials from */
    char pspfile[MAX_PATH];
    
    int freeze_orbital_step;

    /* override occupations */
    int override_occ;

    /* Override of positions (during restart) */
    int override_atoms;

    /* Number of steps after which to output results */
    int outcount;

    /* Total points for potential */
    int vh_nbasis;

    double BT;
   
    /* Total number of electrons */
    double num_el;
   
    /* movie flags */
    int rmvmovie, chmovie, xbsmovie;

    /* Desired vector for constrained dynamics */
    double cd_vector[3];

    /* number of velocity in constrained dynamics */
    double cd_velocity;



    double Evxcold_rho;
    double Evhold_rho;
    double Evh_rho;
    double Evh_rhoc;

    double TOTAL_former;
    double dE;

    double *energies;
    double pulay[10];

    int restart_mix;

    int move_centers_at_this_step;



    int num_storage_st;


    /* "Center" of global function */

    double Bc;
    double Bx;
    double By;
    double Bz;
    double Ac;
    double Ax;
    double Ay;
    double Az;
    double Axy;
    double Axz;
    double Ayz;

    int state_per_proc;
    int state_begin;
    int state_end;

    int max_orbit_size;
    int max_orbit_nx;
    int max_orbit_ny;
    int max_orbit_nz;

    bool movingCenter;
    bool bandwidthreduction;
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


    bool metal;
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

    int is_gamma;
    bool is_use_symmetry;

    bool freeze_occupied;

    // Maximum number of valence electrons for any atomic species
    double max_zvalence;

    // Non-local block size
    int non_local_block_size;
    int poisson_solver;
    int dipole_corr[3];

    // Flag to use fine grid for vdf-df
    bool use_vdwdf_finegrid;

    // Flag indicating whether or not to localize the non-local projectors
    bool localize_projectors;
    bool localize_localpp;

    // Flag indicating whether or not to calculate dipole moment for the entrire cell. Default is off
    bool dipole_moment;
    int nsym, *sym_rotate;
    double *sym_trans;

    // In case system has numa whether or not to use it
    // (may not want to try setting it up internally since the user may want to
    // to do it manually with numactl or aprun
    bool use_numa;

    // Default is false. RMG will still be able to use transparent huge pages but
    // certain special optimizations will be disabled. If you set this to true then
    // RMG assumes that sufficient huge pages are available to meet all memory
    // requirements and bad results may occur if that is not true.
    bool require_huge_pages;

    // Controls how far below the Nyquist frequency potentials are cutoff. Default is 1.0
    double filter_factor;

    // Use square root technique for density interpolation
    bool sqrt_interpolation;

    // Filter density dependent potentials
    bool filter_dpot;

    // Alternate laplacian discretization
    bool alt_laplacian;

} CONTROL;


/* Extern declaration for the main control structure */
#if __cplusplus
extern "C" CONTROL ct;
#else
extern CONTROL ct;
#endif

#endif
