#ifndef RMG_control_H
#define RMG_control_H 1
#include <vector>
#include <complex>

#if CUDA_ENABLED
    #include <cusolverDn.h>
    #include <cublas_v2.h>
    #include <cublasXt.h>
    #if USE_NCCL
        #include <nccl.h>
    #endif
#endif

#if HIP_ENABLED
    #include <hipblas/hipblas.h>
    #include <rocsolver/rocsolver.h>
    #include <hip/hip_runtime_api.h> // for hip functions
    #include <hipsolver/hipsolver.h> // for all the hipsolver C interfaces and type declarations
    #if USE_NCCL
        #include <rccl/rccl.h>
    #endif
#endif

#if SYCL_ENABLED
    #include <CL/sycl.hpp>
    #include <sycl/queue.hpp>
#endif

#include "Klist.h"

/** @name CONTROL
  @memo Main control structure
 
  This is a global structure declared as extern CONTROL ct
 
 */
class CONTROL
{
public:
    void save_args(int argc, char **argvin)
    {
        for (int i = 1; i < argc; i++) this->argv.push_back(argvin[i]);
    }

    // argc and argc
    std::vector<std::string> argv;

    // Branch type
    int rmg_branch;

    // Energy output units 0=Hartrees, 1=Rydbergs
    int energy_output_units;
    double energy_output_conversion[2];
    std::string energy_output_string[2];

    // Special variables when running with multiple images per node
    // The number of images stacked on a single node
    int images_per_node;

    // Image node index ranging from 0 to (images_per_node-1)
    int image_node_id;

    // OpenMP threads per node
    int OMP_THREADS_PER_NODE;

    // Multigrid threads per node. Normally the same as above but may be tuned for performance.
    int MG_THREADS_PER_NODE;

    // requested mpi thread level
    int mpi_threadlevel;

    // mpi queued mode
    bool mpi_queue_mode;
    bool spin_manager_thread;
    bool spin_worker_threads;

    bool write_orbital_overlaps;

    /** Description of the run. */
    std::string description;

    /* time at which run started */
    double time0;
    
    /* Flag indicating if this is a spin polarized calculation */
    bool spin_polarization;
    bool noncoll;
    bool spinorbit;
    int noncoll_factor;
    int nspin;
    bool AFM;

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

    /** whether input and output files are compressed */
    bool compressed_infile;
    bool compressed_outfile;

    /** whether to mmap the weights for the projectors weights, work space and orbitals */
    bool nvme_weights;
    bool nvme_work;
    bool nvme_orbitals;
    int nvme_orbital_fd;
    int nvme_work_fd;
    int nvme_weight_fd;
    std::string nvme_weights_path;
    std::string nvme_work_path;
    std::string nvme_orbitals_path;
    std::string qfunction_path;

    /** Input file name to read wavefunctions from when doing a restart */
    char infile[2*MAX_PATH];

    /** Input file name to write wavefunctions to */
    /* Output file name */
    char outfile[2*MAX_PATH];

    char infile_tddft[2*MAX_PATH];
    char outfile_tddft[2*MAX_PATH];
    bool restart_tddft;
    int tddft_mode;
    double tddft_qpos[3];
    double tddft_qgau;
    

    /** Prepended to pseudopotential name */
    std::string pseudo_dir;

    /** Initial run flag. Read from the input file. 0=initial run otherwise a restart */
    int runflag;

    /* output z-average of states */
    int zaverage;

    /* number of state to output */
    int plot_state;

    /* Exchage-Correlation flag */
    int xctype;

    /* Hybrid EXC flag */
    int xc_is_hybrid;
    double exx_fraction;
    int vdw_corr;
    int dftd3_version;
    double Evdw;
    double ldaU_E;
    double *force_vdw=NULL;
    double stress_vdw[9];

    /* EXX integral file */
    std::string exx_int_file;
    bool exx_int_flag;

    /* Exx computation mode */
    int exx_mode;
    int exxchol_max=8;
    bool ExxIntChol;
    int exxdiv_treatment;
    bool gamma_extrapolation;

    /* Exx convergence multiplier. Threshold checked in the inner (Scf) is multiplied by this */
    double exx_convergence_factor;

    /* EXX potential RMS[dV] */
    double vexx_rms;

    /** Boundary condition flag. Read from the input file. 0=periodic, 1=cluster, 2=surface */
    int boundaryflag;

    /* Coordinate input flag: crystal or cartesian */
    int crd_flag;

    /* Maximum number of MD steps */
    int max_md_steps;

    /* Maximum number of rmg meta loops (NEB, ARTS, etc.) */
    int max_neb_steps;

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
    int max_exx_steps;

    /* Total number of SCF steps done */
    int total_scf_steps;

    /* SCF steps iterator */
    int scf_steps;

    /* EXX steps iterator */
    int exx_steps;

    /* RMS[dV] convergence criterion */
    double thr_rms;

    /* estimated total energy  convergence criterion */
    double thr_energy;  // input value
    double adaptive_thr_energy;  // value checked against in Scf.cpp

    double thr_stress; 
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

    /** Flag for write the pseudopotential  */
    bool write_pp_flag;

    /** Number of states. May switch between init_states and run_states */
    int num_states;
    int tddft_start_state;

    /** Number of initialization states */
    int init_states;

    /** Extra random states used for LCAO start above and beyond the atomic orbitals. */
    int extra_random_lcao_states;

    /** Some pseudopotentials contain unbound atomic orbitals and this flag indicates */
    /** whether or not they should be used for LCAO starts. */
    bool lcao_use_empty_orbitals;

    /** Number of elpa or scalapack groups to use */
    int subdiag_groups;

    /** RMG normally writes parallel restart files. These require that restarts have the */
    /** same processor topology. If write_serial_restart = \"true\" then RMG will also */
    /** write a serial restart file that can be used with a different processor topology. */
    bool read_serial_restart;
  
    /** Directs RMG to read from serial restart files. Normally used when changing */
    /** the sprocessor topology used during a restart run */
    bool write_serial_restart;

    /** If true also implies write_serial_restart */
    bool write_qmcpack_restart;
    bool write_qmcpack_restart_localized;
    int qmc_nband;

    /** Number of run states */
    int run_states;

    /** Maximum number of states, used for davidson */
    int max_states;

    /** Davidson multiplicative factor */
    int davidx;

    /** Davidson maxits */
    int david_max_steps;

    /** Davidson pre multigrid steps */
    int davidson_premg;

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

    /** TF Ion structure */
    /** Used for simple solvent model */
    TF_ION *tf_ions;

    /** Number of species */
    int num_species;

    // potential grid refinement
    int FG_RATIO;

    /* Cutoff parameter */
    double cparm;        // Multiplicative factor
    double rhocparm;     // For qfunctions and density

    /** Total conpensating charge density */
    double crho;

    /** GGA guard value for low density regions */
    double epsg_guard;

    /** Total charge in supercell */
    double tcharge;

    /** Norm conserving pseudo potential flag */
    int norm_conserving_pp;

    /** Semi local pseudo potential flag */
    int semilocal_pp;
    bool use_bessel_projectors;

    /* Kohn-sham finite difference order */
    int kohn_sham_fd_order;
    int force_grad_order;
    bool kohn_sham_ke_fft;

    /* Prolong operator default order */
    int prolong_order;
    /* Prolong operator mixing parameter */
    double cmix = 0.0;
    /* Flag indicating whether or not to use cmix */
    bool use_cmix = true;


    // Flag indicating whether or not to use gpu finite differencing for the hamiltonian
    bool use_gpu_fd;

    /* This is the max of nldim for any species cubed */
    int max_nlpoints;
    int max_lpoints;
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
    double drho_q0;
    int drho_precond_type;

    /** Minimum charge density after mixing */
    double min_rho;

    bool charge_pulay_Gspace;
    bool drho_precond;
    /*Order of Pulay mixing for charge density*/
    int charge_pulay_order;

    /*How often to refresh Pulay history*/
    int charge_pulay_refresh;

    /*Order of Broyden mixing for charge density*/
    int charge_broyden_order;

    /* Resta beta parameter. Should be conventional unit cell A0 */
    double resta_beta;

    int ldau_mixing_type;
    int ldau_pulay_refresh;
    int ldau_pulay_order;
    double ldau_mix;
    double ldau_pulay_scale;
    

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
    int dos_flag;
    double gaus_broad;

    double occ_width;

    double occ_mix;
    int mp_order;

    /** total background smearing charge -- for charged supercells */
    double background_charge;
    double system_charge;


    /** Multigrid parameters for the eigenvalue solver */
    MG_PARM eig_parm;

    /** Multigrid parameters for the poisson solver */
    MG_PARM poi_parm;

    /** maximum multigrid level where offset routines can be used */
    int mg_offset_level;
 
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

    /* Use folded spectrum if possible */
    bool use_folded_spectrum;

    /* Folded spectrum width */
    double folded_spectrum_width;

    /* Folded spectrum iterations */
    int folded_spectrum_iterations;

    /* Expansion factor for non-local projectors */
    double projector_expansion_factor;

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

    /** Ionic motion timestep */
    double iondt;

    //  TDDFT time step in atomic unit
    double tddft_time_step;

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
    double xcstate;
    double vtxc;
    double NL;
    double II;
    double FOCK;
    double TOTAL;

    // EXX convergence measure
    double exx_delta;

    // Exit criteria from EXX loop
    double exx_convergence_criterion;

    // Single/double precision delta fft threshold for EXX Vexx computation
    double vexx_fft_threshold;

    /* fermi energy */
    double efermi;

    /* Tetrahedron method. 0=Bloechl, 1=Linear, 2=Optimized */
    int tetra_method;

    /** Total number of k-points being used in the calculation */
    int num_kpts;
    int num_kpts_pe;
    Klist klist;



    /** minimal K-point structure 
     *  the Kpoint class stores a reference to this so don't change this unless
     *  you know what you are doing.
    */
    std::vector<KSTRUCT> kp;
    int kpoint_mesh[3];
    int qpoint_mesh[3];
    int kpoint_is_shift[3];

    /** The maximum number of projectors for any species */
    int max_nl;
    
    /** The maximum number of atomic orbitals for any species */
    int max_orbitals;

    /** The maximum number of LDA+U orbitals used for projections for any species */
    int max_ldaU_orbitals;

    /** Type of atomic orbitals used for LCAO inits and projections */
    int atomic_orbital_type;

    /** If interntal pseudopotential is wanted which type */
    int internal_pseudo_type;

    /** Gygi all-electron parameter **/
    int all_electron_parm;

    /** Number of semilocal projectors */
    int semilocal_projectors;

    /** The maximum l quantum number for any species */
    int max_l;

    /** Clebsch-Gordon coefficients */
    double_3d_array cg_coeff;

    /*Maximum value of nldim for any species */
    int max_nldim;

    double max_nlradius, min_nlradius, max_qradius, min_qradius;

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
    double ns_occ_rms;

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

    // Some GPU information.
    // Total number of gpu devices present in the node
    int num_gpu_devices;

    // Flag to indicate whether or not to use gpu managed memory
    bool gpu_managed_memory;

    // Flag to use new energy correction terms
    bool use_energy_correction;

#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED


    // Total number of usable gpu devices present in the node
    int num_usable_gpu_devices;

    // Device ids for the usable gpus
    int gpu_device_ids[MAX_GPU_DEVICES];

    // GPU memory for the usable devices
    size_t gpu_mem[MAX_GPU_DEVICES];

    // Default is to use managed memory for non-local weights but if GPU memory
    // is constrained performance is much better using pinned memory.
    bool pin_nonlocal_weights;

    // Flag indicating whether all of the gpu devices we plan on using support managed memory
    bool gpus_support_managed_memory;

#endif

#if CUDA_ENABLED

    // Cuda version
    int cuda_version;

    // Cuda devices
    CUdevice  cu_dev;
    CUdevice  cu_devices[MAX_GPU_DEVICES];

    // Cuda device context
    CUcontext cu_context;

    // CUBLAS library handles
    cublasHandle_t cublas_handle;
    cublasHandle_t gpublas_handle;
    cublasXtHandle_t cublasxt_handle;
    cusolverDnHandle_t cusolver_handle;
    cudaStream_t cusolver_stream;
    bool use_cublasxt;
#if USE_NCCL
    ncclUniqueId nccl_nd_id;  
    ncclComm_t nccl_world_comm;
    ncclComm_t nccl_local_comm;
#endif
#endif

#if HIP_ENABLED

    // Cuda version
    int hip_version;

    // Hip devices
    hipDevice_t hip_dev;
    hipDevice_t hip_devices[MAX_GPU_DEVICES];
    int smemSize[MAX_GPU_DEVICES];

    // Hip device context
    hipCtx_t hip_context;

    // hipblas library handles
    hipblasHandle_t hipblas_handle;
    hipblasHandle_t gpublas_handle;
    hipStream_t rocsolver_stream;
    rocsolver_handle roc_handle;
    hipsolverHandle_t hipsolver_handle;
    bool use_cublasxt;
#if USE_NCCL
    ncclUniqueId nccl_nd_id;  
    ncclComm_t nccl_world_comm;
    ncclComm_t nccl_local_comm;
#endif
#endif

#if SYCL_ENABLED
    std::vector<cl::sycl::device> sycl_devs;
    cl::sycl::queue sycl_Q;
    int host_dev;
#endif

    /* RMG2BGW options */
    bool rmg2bgw;
    double ecutrho, ecutwfc;

   
    int vxc_diag_nmin;
    int vxc_diag_nmax;

    int freeze_orbital_step;
    int exx_freeze_orbital_step;

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

    int restart_mix;

    int move_centers_at_this_step;


    /* "Center" of global function */
    double Ac;

    int state_per_proc;
    int state_begin;
    int state_end;
    int num_orbitals_total; // including the orbitals stored in the process and orbitals overlaped with them.
    int *orbitals_list;  // dimension of ct.num_orbital_total, sorted.
    int max_orbit_size;
    int max_orbit_nx;
    int max_orbit_ny;
    int max_orbit_nz;

    bool movingCenter;
    bool bandwidthreduction;
    int movingSteps;
    STATE *states;


    bool cube_rho;
    bool cube_vh;
    bool cube_pot;
    std::vector<int> cube_states_list;
    std::vector<double> stm_bias_list;
    std::vector<double> stm_height_list;
    
    std::vector<std::string> file_atomic_orbit;

    /* output information for GW calculations.  --Xiaohong */
    int flag_gw;

    /*  for negf  */


    bool metal;
    int num_blocks;
    int block_dim[MAX_BLOCKS];
    std::vector<int> block_dim_phi, block_dim_nl;

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

    bool is_gamma;
    int is_use_symmetry;
    bool time_reversal;
    bool frac_symm;

   bool wannier90;
   int wannier90_scdm;
   double wannier90_scdm_mu;
   double wannier90_scdm_sigma;
   int num_wanniers;


   bool freeze_occupied;

   // Maximum number of valence electrons for any atomic species
   double max_zvalence;

   // Non-local block size
   int non_local_block_size;
   int poisson_solver;
   int dipole_corr[3];

   // Flag to use fine grid for vdf-df
   bool use_vdwdf_finegrid;

   // vdW kernel table file
   std::string vdW_kernel_file;

   // Flag indicating whether or not to localize the non-local projectors
   bool localize_projectors;
   bool localize_localpp;
   bool proj_nophase;

   // Flag indicating whether or not to calculate dipole moment for the entrire cell. Default is off
   bool dipole_moment;
   int nsym;

   // In case system has numa whether or not to use it
   // (may not want to try setting it up internally since the user may want to
   // to do it manually with numactl or aprun
   bool use_numa;

   // In case system has hwloc whether or not to use it
   bool use_hwloc;

   // Some CPU versions of zgemm are very slow because of compiler implementations
   // of std::complex. In that case setting this flag lets you use an alternate implementation.
   bool use_alt_zgemm;

   // Default is false. RMG will still be able to use transparent huge pages but
   // certain special optimizations will be disabled. If you set this to true then
   // RMG assumes that sufficient huge pages are available to meet all memory
   // requirements and bad results may occur if that is not true.
   bool require_huge_pages;

   // Controls how far below the Nyquist frequency potentials are cutoff. Default is 0.25
   double filter_factor;

   // Use square root technique for density interpolation
   bool sqrt_interpolation;

   // Use a faster but less accurate method for generating the charge density
   bool fast_density;

   // Alternate laplacian discretization
   bool alt_laplacian;

   // Value for CFAC read from input file instead of being computed
   double afd_cfac;

   // Flag is true if the ddd is non-diagonal for any atomic species
   bool is_ddd_non_diagonal;

   // Flag controlling whether or not to use asynchronous MPI allreduce operations in certain places.
   // Async reduction is not available on all platforms so there is both a compilation flag and a
   // runtime flag. In particular on a GPU platform that supports asynchronous all reduce operations
   // on CPU's only you want to set this to false.
   bool use_async_allreduce;

   bool ON_read_from_RMG;
   char infile_ON_from_RMG[2*MAX_PATH];
   int freeze_rho_steps;
   bool use_cpdgemr2d;

   int orbital_mixing_method;
   /*Order of Pulay mixing for orbital*/
   int orbital_pulay_order;

   /*How often to refresh Pulay history for orbital mixing*/
   int orbital_pulay_refresh;

   /*Scale parameter for residuals in Pulay orbital mixing*/
   double orbital_pulay_scale;

   double orbital_pulay_mixfirst;

   int fd_allocation_limit;

   // LDA+U options
   int ldaU_mode;
   int num_ldaU_ions;
   int max_ldaU_l;  // max ldaU l value for any ion
   int freeze_ldaU_steps;


   bool stress;
   bool cell_relax;
   int cell_movable[9];
   double ldaU_radius;
   double stress_tensor[9];

   // Memory usage options
   size_t q_alloc[3];
   size_t beta_alloc[3];
   size_t psi_alloc[3];
   size_t vexx_alloc[3];
   int LocalizedOrbitalLayout=1;

   std::string input_initial, input_final;
   double totale_initial, totale_final;

   // Whether or not forces and stress have converged or not.
   bool forces_converged;
   bool is_converging;

   // bfgs options
   int bfgs_ndim = 1;
   double trust_radius_max = 0.8;
   double trust_radius_min = 1.0e-4;
   double trust_radius_ini = 0.5;
   double w_1 = 0.01;
   double w_2 = 0.5;

   // Testing options
   double test_energy=NAN;
   double test_energy_tolerance=1.0e-7;
   double test_bond_length=NAN;
   double test_bond_length_tolerance=1.0e-3;
   int test_steps=0;
   int test_steps_tolerance=0;
   int kpoint_units = 0;

};


/* Extern declaration for the main control structure */
extern CONTROL ct;
extern std::vector<ION> Atoms;
extern std::vector<SPECIES> Species;
#endif
