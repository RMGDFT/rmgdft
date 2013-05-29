/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*** QMD-MGDFT/main.h *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */

typedef struct
{
    int count;
} QMD_sem_t;



/** @name PE_CONTROL
  *
  * @memo Processor control structure.
  * 
  * This is a global structure declared as extern PE_CONTROL pct.
  * 
 */

typedef struct
{



  /** Number (rank in MPI terminology) of this processor in this image
 *  * grid */
    int gridpe, imgpe, thisimg, spinpe;

    /** Number of grids (typically 1) per image to be run simultaneously
 *  * **/
    int images, grids;

    /* MPI communicators for each code grid (grid_comm) and one
 *  * (rmg_comm)
 *   *      * for all group rank 0 pe's. The later effectively replaces
 *    *      MPI_COMM_WORLD
 *     *           * unless you really need all-to-all, even across
 *     grids,
 *      *           communication. */
   MPI_Comm rmg_comm, img_topo_comm, grid_topo_comm, grid_comm, img_comm, spin_comm;


    /* Grid sizes on each PE */
    int PX0_GRID;
    int PY0_GRID;
    int PZ0_GRID;

  /* Grid offsets on each PE */
    int PX_OFFSET;
    int PY_OFFSET;
    int PZ_OFFSET;

    /* Basis size on each PE */
    int P0_BASIS;

    /* Fine grid sizes on each PE */
    int FPX0_GRID;
    int FPY0_GRID;
    int FPZ0_GRID;

    /* Fine Grid offsets on each PE */
    int FPX_OFFSET;
    int FPY_OFFSET;
    int FPZ_OFFSET;


    /* Fine grid basis size on each PE */
    int FP0_BASIS;

    /** Number (rank in MPI terminology) of this processor */
    int instances;
  
    /** Neighboring processors in three-dimensional space */
    int neighbors[6];
  
    /** Processor kpoint- coordinate for domain decomposition */
    /*  paralleled for kpoint */
    int pe_kpoint;
  
    /** Processor x-coordinate for domain decomposition */
    int pe_x;
    /** Processor y-coordinate for domain decomposition */
    int pe_y;
    /** Processor z-coordinate for domain decomposition */
    int pe_z;
  
    /** An index array which maps the projectors onto the 3-d grid associated 
     * whith each processor */
    int *Qindex[MAX_IONS];

    /** An index array which indicate whether the grid map on the current processor */
    int *Qdvec[MAX_IONS];

    /** Number of points in the Qindex array for each ion */
    int Qidxptrlen[MAX_IONS];

    /** Number of points in the circle of local projector for eac processor*/
    int lptrlen[MAX_IONS];
  
    /** Point to start of storage for the augument function */
    REAL *augfunc[MAX_IONS];

    /** Point to start of DnmI function storage for this ion */
    REAL *dnmI[MAX_IONS];
    REAL *dnmI_x[MAX_IONS];
    REAL *dnmI_y[MAX_IONS];
    REAL *dnmI_z[MAX_IONS];

    /** Point to start of qqq storage for this ion */
    REAL *qqq[MAX_IONS];

    /** local grid size for x,y,z **/
    int nx_grid; 
    int ny_grid; 
    int nz_grid; 

    /* kpoint index for start and end for a subdomain processors */
    int kstart;
    int kend;  
  
    /*  processor coordinates in COMM_KP communicator */
    int coords[2];

    /* Number of ions centered on this processor */
    int n_ion_center;
    int n_ion_center_loc;

    /* Indices of the ions centered on this PE */
    int ionidxcenter[IONS_PER_PE];
  
    /* Projectors per ion in a given region*/
    int prj_per_ion[MAX_IONS];
  
    /* Indices of the ions within non-local range */
    int ionidx[IONS_PER_PE];
    int ionidx_loc[IONS_PER_PE];
  
    /* Points to start of projectors for this ion in projector space */
    /* All projectors are stored consecutively.                      */
    int prj_ptr[IONS_PER_PE];
  
    /* Pointers to the index array for the non-local potentials */
    int idxptr[IONS_PER_PE];
  
    /* The size of local orbitals on this PE */
    int   psi_size;
    /* pointer to former step solution, used in pulay and KAIN mixing  */
    double *former_psi;
    double *psi;
    double *psi1; 
    double *psi2; 

    double *projectors; 

    int desca[DLEN];
    int descb[DLEN];
    int mycol;
    int myrow;
    int nprow;
    int npcol;
    int num_local_orbit;

} PE_CONTROL;


/**@name STATE
 *
 * @memo Wavefunction storage structure */
typedef struct
{
  
    /** First iteration flag */
    int firstflag;
  
    /** Current estimate of the eigenvalue for this orbital (state). */
    REAL eig;
  
    /** Wavefunction residual error computed by multigrid solver */
    REAL res;
  
    /** Points to the storage area for the real part of the orbital */
    REAL *psiR;
    /** Points to the storage area for the imaginary part of the orbital */
    REAL *psiI;
  
    /** Occupation of the orbital */
    REAL occupation;
    REAL oldeig;
  
    /** Index showing which k-point this orbital is associated with */
    int kidx;
    int istate;

    /* The ion on which this state is localized */
    int inum;
  
    /* index for functions with same localization */
    int loc_index;

    /* Actual Physical coordinates at current time step */
    int  pe;
    REAL crds[3];
    REAL radius;
    int movable;
    int frozen;
    int index;

    int ixmin;
    int ixmax;
    int iymin;
    int iymax;
    int izmin;
    int izmax;
    int xfold;
    int yfold;
    int zfold;
    int ixstart;
    int iystart;
    int izstart;
    int ixend;
    int iyend;
    int izend;
    int orbit_nx, orbit_ny, orbit_nz;
    int size;
    /* Localization mask */
    char *lmask[4];
    int ixmin_old;
    int ixmax_old;
    int iymin_old;
    int iymax_old;
    
    int atomic_orbital_index;
    
} STATE;



/**@name SPECIES
 * @memo Species (pseudopotential) control structure
 * @doc Structure holds data about the pseudopotentials used to
 * represent different types of atomic species. 
*/
typedef struct 
{

    /** Description of the species (e.g Atomic carbon generated using 
     * hamann's code with  rcs=0.80 rcp=0.85 bohr
     */

    /* symbol read from control file */
    char pseudo_symbol[32];

    /* pseudopotential filename */
    char pseudo_filename[MAX_PATH];

            
    char description[256];
  
    /** Atomic number */
    int atomic_number;
  
    /** Atomic symbol */
    char *atomic_symbol;

    /** Atomic mass */
    REAL atomic_mass;
  
    /** Number of valence electrons */
    REAL zvalence;
  
    /** Gaussian charge parameter used for compensating the 1/r Coulomb
     * tail of the pseudopotentials
     */
  
    REAL rc;
  
    /* Number of grid points in the local and non-local potential localized grids
     * in each coordinate direction. These used to be L0_LDIM and L0_NLDIM.
     */
    int ldim;
    int nldim; 
    int ldim_coar;
    int nlfdim; 
    int qdim;

  
    /* These are input parameters in the pseudopotential file. They represent the
     * real radii that are used in generating ldim and nldim.
     */
    REAL lradius;
    REAL nlradius;
    REAL qradius;


    /*Radius for milliken analysis*/
    REAL mill_radius;
    /*Radius in number of grid points*/
    int mill_dim;
    /*l-numbers for states for which we have atomic orbitals*/
    int lstate_atomic_wave[5];
    /*Sum of all atomic states (with different l or m numbers*/
    int sum_atomic_waves;



    /*This will store name of atomic wavefunctions, such as s, px, py, pz, dxx etc*/
    char atomic_wave_symbol[20][12];




    /** Number of radial grid points in the pseudopotential file */
    int rg_points;

    /* Log mesh parameter 
       REAL almesh;*/

    /* Log mesh parameter, where aa=exp(-aasf)/Z, bb=1.0/bbsf */
    REAL aa,bb;
    /** Non-linear core correction flag */
    int nlccflag;

    /* Number of potentials */
    int num_potentials;

    /* L-values for the reference states */
    int lval[10];

    /* L-value for local pseudopotential state */
    int local;

    /* Index for local pseudopotential state */
    int localidx;

    /* Radial grids */
    /*REAL r[MAX_L+1][MAX_RGRID];*/

    /* Reference states */
    REAL psi[MAX_L+1][MAX_RGRID];


    /* Pseudopotentials */
    REAL psp[MAX_L+1][MAX_RGRID];

    /* Number of grid points in the beta function */
    int kkbeta;

    /* Matrix ddd0(nbeta,nbeta)*/
    REAL ddd0[18][18];
    REAL ddd[18][18];

    /* Matrix qqq(nbeta,nbeta)*/
    REAL qqq[18][18];

    /* The number of L=|l1-l2|.....|l1+l2|, we limit nlc <=5 */
    int nlc;

    /* The number of component in polynomial of the pseudized Q_I(r) function we limit nqf<=10 */
    int nqf;

    /* L-independent inner cutoff radii rinner for Q_I(r) function */
    REAL rinner[5];

    /* ultrosoft Vanderbilt Qnm_rad(r) function and */
    REAL *qnm;
    REAL *qnmlig;
    REAL *drqnmlig;

    /* the coefficient for pseudosation of Qnm_L(r)*/
    REAL *qfcoef;

    /* logarithmic radial mesh information */
    REAL r[MAX_RGRID];
    REAL rab[MAX_RGRID];

    /* Local Pseudopotentials */
    REAL vloc0[MAX_RGRID];

    /* Core charge radial grids */
    REAL cr[MAX_RGRID];

    /* Pseudo atomic valence density */
    REAL avdens[MAX_RGRID];

    /* Pseudo atomic core density */
    REAL rspsco[MAX_RGRID];

    /* Kleinman-Bylander Projectors on radial grid */
    /*REAL kbp[MAX_L+1][MAX_RGRID];*/

    /* The L-value for the beta function */
    int llbeta[MAX_NB];

    /* Ultrasoft Vanderbilt beta_n(r) function on radial grid */
    REAL beta[MAX_NB][MAX_RGRID];


    /* Total number of projectors */
    int nbeta;

    /* Radial derivative terms for Kleinman-Bylander projectors */
    REAL drkbp[MAX_L+1][MAX_RGRID];


    /* Local difference potential on radial grid */
    REAL difflocal[MAX_RGRID];


    /* Radial derivative of local difference potential */
    /*REAL drdifflocal[MAX_RGRID];*/


    /* Total number of projectors */
    int ipcount;


    /* Normalization coefficients for the KB projectors */
    REAL kbnorm[MAX_L+1];


    /* Milliken projector normalization coefficients */
    REAL mnorm[MAX_L + 1];


    /* Linear interpolation storage for the compensated local potential
     * and for it's radial derivative.
     */
    REAL localig[MAX_LOCAL_LIG];
    REAL drlocalig[MAX_LOCAL_LIG];

    /* Linear interpolation storage for the core charge density */
    REAL rhocorelig[MAX_LOCAL_LIG];

    /* Kleinman-Bylander Projectors on linear interpolation grid */
    REAL kbplig[MAX_L+1][MAX_LOCAL_LIG];

    /* Radial derivatives of the Kleinman-Bylander Projectors on linear interpolation grid */
    REAL drkbplig[MAX_L+1][MAX_LOCAL_LIG];

    /* Reference states on linear interpolation grid */
    REAL psilig[MAX_L+1][MAX_LOCAL_LIG];

    /* Ultrasoft PP Projectors on linear interpolation grid */
    REAL betalig[MAX_NB][MAX_LOCAL_LIG];

    /* Radial derivatives of the Ultrasoft PP Projectors on linear interpolation grid */
    REAL drbetalig[MAX_NB][MAX_LOCAL_LIG];

    /* Local potential linear interpolation grid spacing */
    REAL drlig;


    /* Non-local linear interpolation grid spacing */
    REAL drnlig;

    /* Qfunction linear interpolation grid spacing */
    REAL drqlig;


    /* Pseudopotential filtering parameters */
    REAL lrcut;		/* Real space local cutoff */
    REAL nlrcut[5]; 		/* Real space non-local cutoff */
    REAL rwidth;            /* Real-space width parameter */
    REAL gwidth;            /* G-space width parameter */

    /* The total number of projectors for the nonlocal part of the pseudopotential*/
    int num_projectors;

    fftw_complex *forward_beta;
    fftw_complex *forward_derbeta_x;
    fftw_complex *forward_derbeta_y;
    fftw_complex *forward_derbeta_z;

    char *backward_wisdom;

    /*Some parameters for Q function*/
    int indv[18];
    int nhtol[18];
    int nhtom[18];
    int nh;
    
    /*Filtering parameters for atomic wavefunctions and charge density*/
    REAL acut; 
    REAL aradius; 
    REAL agwidth;
    REAL arwidth;
    
    
    /* radius of atomic wavefunctions and charge in terms of number of grid points*/
    int adim_rho;
    int adim_wave;

    /*Number of radial atomic wave functions - these depend on l only, not on m*/
    int num_atomic_waves;
    /*l-numbers for states for which we have atomic orbitals*/
    int atomic_wave_l[5];
    REAL atomic_wave_oc[5];
    
    char atomic_wave_label[5][3];

    REAL *atomic_rho;
    
    /* Pseudo atomic valence density read from PP file in log grid*/
    REAL **atomic_wave;
    /* Pseudo atomic valence density on linear grid*/
    REAL **awave_lig;

} SPECIES;



/* Ion structure */
typedef struct 
{

    /* Initial physical coordinates at start of run */
    REAL icrds[3];

    /* Actual Physical coordinates at current time step */
    REAL crds[3];

    /* Positions at the previous time step */
    REAL ocrds[3];

    /* Initial crystal coordinates at start of run */
    REAL ixtal[3];

    /* Actual crystal coordinates at current time step */
    REAL xtal[3];

    /* Crystal coordinates  at the previous time step */
    REAL oxtal[3];

    /* Coordinates of the corner of the grid that the local */
    /* difference potential is nonzero on.                  */
    REAL lxcstart;
    REAL lycstart;
    REAL lzcstart;


    /* Coordinates of the corner of the grid that the non-local */
    /* potential is nonzero on.                                 */
    REAL nlxcstart;
    REAL nlycstart;
    REAL nlzcstart;

    REAL xcstart_loc;
    REAL ycstart_loc;
    REAL zcstart_loc;

    /* Coordinates of the corner of the grid that the non-local */
    /* potential is nonzero on.                                 */
    REAL Qxcstart;
    REAL Qycstart;
    REAL Qzcstart;


    /* Integer species type when using a raw pseudopotential */
    int species;


    /* Forces on the ion */
    REAL force[4][3];

    /* Current velocity of the ion */
    REAL velocity[3];

    /* Kleinman-Bylander normalization coefficients */
    REAL pd[(MAX_L + 1) * (MAX_L + 1)];

    /* Milliken normalization coefficients */
    REAL mnorm[(MAX_L + 1) * (MAX_L + 1)];

    /* Total number of projectors */
    int prjcount;

    /* Movable flag */
    int movable;

    int frozen;

    /*  number of local orbitals on the ion */
    int n_loc_states;

    /* Localization mask */
    char *lmask[4];

    int ixstart;
    int iystart;
    int izstart;
    int ixend;
    int iyend;
    int izend;

    int ixstart_loc;
    int iystart_loc;
    int izstart_loc;
    int ixend_loc;
    int iyend_loc;
    int izend_loc;


    int first_state_index;



 /* Force modifier parameters */
     struct {
         REAL setA_weight;
         REAL setA_coord[3];
         REAL setB_weight;
         REAL setB_coord[3];
         double forcemask[3];
     } constraint;


} ION;






/** @name KPOINT
 * @memo Holds data specific to individual k-points.
 */                                    
typedef struct 
{

    /** The index of the k-point for backreferencing */
    int kidx;

    /** The k-point */
    REAL kpt[3];

    /** The corresponding vector */
    REAL kvec[3];

    /** The weight associated with the k-point */
    REAL kweight;

    /** The magnitude of the k-vector */
    REAL kmag;

    STATE *kstate;						     

    /* Mean min, and max wavefunction residuals for occupied space */
    REAL meanres;
    REAL minres;
    REAL maxres;

} KPOINT;


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

    pid_t main_thread_pid;
    int THREADS_PER_NODE;

    /** Description of the run. */
    char description[200];

    int mg_eig_precision;
    int mpi_threadlevel;
    /** Name of the input control file. Passed as a command line argument
     *
     *  Example:
     *  bash$  md in.diamond8
     */
    double time0;
    int spin_flag;
 /* determine whether to initialize up and down density equally or not
 * */
    int init_equal_density_flag;
    /* determine whether to get_dos or not */
    int pdos_flag;


    char cfile[MAX_PATH];
    char basename[MAX_PATH];

    FILE *logfile;

    /** Input file name to read wavefunctions  and potential and rho from when doing a restart */
    char infile[MAX_PATH];
    char infile_rho[MAX_PATH];

    /** Input file name to write wavefunctions to */
    /* Output file name */
    char outfile[MAX_PATH];

    /** File to read the pseudopotentials from */
    char pspfile[MAX_PATH];

    /** Initial run flag. Read from the input file.
      0=initial run otherwise a restart */
    int runflag;

    /*  spin polarize flag, 0 = no spin, 1 = spin polarized */
    int spin;

    /* output z-average of states */
    int zaverage;

    /* number of state to output */
    int plot_state;

    /* Exchage-Correlation flag */
    int xctype;

    /** Boundary condition flag. Read from the input file. 
      0=periodic, 1=cluster, 2=surface */
    int boundaryflag;

    /* Coordinate input flag: crystal or cartesian */
    int crd_flag;


    /* Maximum number of MD steps */
    int max_md_steps;

    /* MD steps iterator */
    int md_steps;

    /* Maxium number of SCF steps in an MD step */
    int max_scf_steps;

    /* Actual number of steps done */
    int scf_steps;

    /* Total number of SCF steps done */
    int total_scf_steps;


    /* Maximum number of steps to do */
    int num_steps;

    /* Actual number of steps done */
    int steps;

    /* parameter to define the gate bias */
    double gbias_begin;
    double gbias_end;
    double BT;
    double gate_bias;
    double bg_begin;
    double bg_end;

    /* override occupations */
    int override_occ;

    /* Override of positions (during restart) */
    int override_atoms;

    char occupation_str[256];

    /* convergence criterion */
    REAL thr_rms;

    /* force convergence criterion */
    REAL thr_frc;

    /* Number of steps after which to perform checkpointing */
    int checkpoint;

    /* Number of steps after which to output results */
    int outcount;

    /** Sorting flag for wavefunctions. Read from input file. 0=no sort, 1=sort */
    int sortflag;

    /** Number of states */
    int num_states;

    /** Number of ions */
    int num_ions;

   /*string to store repeat count occupations for spin up*/
    char occupation_str_spin_up[256];

    /*string to store repeat count occupations for spin down*/
    char occupation_str_spin_down[256];


    /** Ion structure */
    ION *ions;

    /** Number of species */
    int num_species;

    /* Cutoff parameter */
    REAL cparm;
    REAL betacparm;
    REAL qcparm;

    /** Total conpensating charge density */
    REAL crho;

    /** Total charge in supercell */
    REAL tcharge;

    /** Species structure 
     * @see SPECIES */
    SPECIES *sp;

    /** Global uniform grid spacing in x */
    REAL hxgrid;

    /** Global uniform grid spacing in y */
    REAL hygrid;

    /** Global uniform grid spacing in z */
    REAL hzgrid;

    /* The fine uniform grid spacing in x */
    REAL hxxgrid;

    /* The fine uniform grid spacing in y */
    REAL hyygrid;

    /* The fine uniform grid spacing in z */
    REAL hzzgrid;

    /* Kohn-sham finite difference order */
    int kohn_sham_fd_order;

    /** bravais lattice type */
    int ibrav;

    /** Lattice information */
    REAL celldm[6];

    /* lattice vectors */
    REAL a0[3];
    REAL a1[3];
    REAL a2[3];

    /** Total cell volume */
    REAL omega;

    /* lengths of the sides of the supercell */
    REAL xside;
    REAL yside;
    REAL zside;

    /* This is the max of nldim for any species cubed */
    int max_nlpoints;
    int max_lpoints;
    int max_Qpoints;

    /** Maximum grid spacing in any coordinate direction */
    REAL hmaxgrid;


    /** Minimum grid spacing in any coordinate direction */
    REAL hmingrid;


    /** Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. 
      A value larger than 1.05 can lead to convergence problems. */
    REAL anisotropy;


    /** Volume element associated with each grid point */
    REAL vel;
    REAL vel_f;


    /** Physical grid basis size */
    int nbasis;


    /** Density mixing parameter. Typical values range from 0.2 to 0.9, while larger values provide faster convergence as long as they are stable. */
    REAL mix;

    /*Order of Pulay mixing for charge density*/
    int charge_pulay_order;

    /*How often to refresh Pulay history*/
    int charge_pulay_refresh;

    /*Scale parameter for residuals in Pulay mixing*/
    REAL charge_pulay_scale;

    /*Flag to test whether or not the modified metrics should be used in
 * Pulay mixing*/
    int charge_pulay_special_metrics;

    /*Weight for Pulay special metrics*/
    REAL charge_pulay_special_metrics_weight;



    /* Projector mixing parameter */
    REAL prjmix;


    /* Global uniform grid corner */
    REAL xcstart;
    REAL ycstart;
    REAL zcstart;


    /* Hartree potential grid sizes per domain */
    int vh_pxgrid;
    int vh_pygrid;
    int vh_pzgrid;

    /* Potential grid sizes */
    int vh_nxgrid;
    int vh_nygrid;
    int vh_nzgrid;


    /* Total points in hartree potential per domain */
    int vh_pbasis;


    /* Wavefunction grid sizes */
    int psi_nxgrid;
    int psi_nygrid;
    int psi_nzgrid;

    int psi_fnxgrid;
    int psi_fnygrid;
    int psi_fnzgrid;

    /* Total points for wavefunctions */
    int psi_nbasis;
    int psi_fnbasis;

    /* Total points for potential */
    int vh_nbasis;

    /* Decoupled hartree potential */
    REAL *vh_ext;


    /* Mean min, and max wavefunction residuals for occupied space */
    REAL meanres;
    REAL minres;
    REAL maxres;

    /* total ionic charge */
    REAL ionic_charge;

    /* Variable occupation stuff */
    REAL nel;

    int occ_flag;

    REAL occ_width;

    REAL occ_mix;

    /** total background smearing charge -- for charged supercells */
    REAL background_charge;


    /** Multigrid parameters for the eigenvalue solver */
    MG_PARM eig_parm;

    /** Multigrid parameters for the poisson solver */
    MG_PARM poi_parm;


    /* Total number of electrons */
    REAL num_el;


    /* force pointer array */
    int fpt[4];

    /* temperature control */
    int tcontrol;

    /* md integrate level */
    int mdorder;

    /* movie flags */
    int rmvmovie,chmovie;

    int domilliken, milliken;

    /* Diagonalization flag and period */
    int initdiag;
    int diag;
    int end_diag;

  /* scalapack block size */
    int scalapack_block_factor;

    /* How many steps between writeout of eigenvalues*/
    int write_eigvals_period;


    /** Force flag. 0=don't compute forces, 1=compute forces */
    int forceflag;


    /* Whether to write full memory usage report at the end of
 * calculation */
    int write_memory_report;

    /* Number of scf steps per md step */

    /* Desired vector for constrained dynamics */
    REAL cd_vector[3];

    /* number of velocity in constrained dynamics */
    REAL cd_velocity;


    /** Ionic motion timestep */
    REAL iondt;


    /** Ionic motion energy */
    REAL ionke;


    /* Total energies */
    REAL ES;
    REAL NUC;
    REAL KE;
    REAL XC;
    REAL NL;
    REAL II;
    REAL TOTAL;
    REAL Evxcold_rho;
    REAL Evhold_rho;
    REAL Evh_rho;
    REAL Evh_rhoc;

    REAL TOTAL_former;
    REAL dE;

    REAL *energies;

    int  restart_mix;

    int   move_centers_at_this_step;


    /* fermi energy */
    REAL efermi;

    /** Total number of k-points being used in the calculation */
    int num_kpts;


    /** K-point control structure */
    KPOINT *kp;                                                     

    /** The maximum number of projectors for any species */
    int max_nl;
    
    /** The maximum l quantum number for any species */
    int max_l;

    REAL cut_radius;
    int num_storage_st;


    /* Index telling if an orbital has been allocated */
    char allocate_st[MAX_STATES];

    /* Number of orbital in each global function */
    int nb_func[MAX_STATES];

    /* Initial storage of each orbital
       (to restart with same function) */
    int alloc_ini[MAX_STATES];

    /* "Center" of global function */
    REAL center[MAX_STATES][3];

    int state_per_proc;
    int state_begin;
    int state_end;
    int ion_begin;
    int ion_end;

    int max_orbit_size;
    int max_orbit_nx;
    int max_orbit_ny;
    int max_orbit_nz;

    int     movingCenter;
    int     movingSteps;
    int     mg_method;
    int     mg_steps;

    STATE *states;

    int metal;
    int num_blocks;
    int block_dim[MAX_BLOCKS];

    int num_cond_curve;
    int *cond_probe1;
    int *cond_probe2;

    int nxfgrid;
    int nyfgrid;
    int nzfgrid;

    char file_atomic_orbit[MAX_SPECIES][MAX_PATH];

    int plane[5];

    /* the external electric field */
    REAL e_field;

    REAL x_field_0;

    REAL y_field_0;

    REAL z_field_0;

    /* Should we process constraints on forces flag */
    int constrainforces;
 
    REAL neb_spring_constant;

  /*Current RMS value*/
    REAL rms;

    /* Max number of sweeps in get_vh*/
    int hartree_max_sweeps;

    /* Min number of sweeps in get_vh*/
    int hartree_min_sweeps;

    /*Ratio between target RMS for get_vh and RMS total potential*/
    REAL hartree_rms_ratio;
   
    int ncpus;

   /* Potential acceleration constant step factor */
    REAL potential_acceleration_constant_step;

    /* Potential acceleration constant step factor */
    REAL potential_acceleration_poisson_step;



    int mask_function;
    int norm_conserving_pp;
    int images_per_node;
    int simpson_depth;
    double simpson_tol;


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

    cuDoubleComplex *gpu_Htri, *gpu_Gtri;
    cuDoubleComplex *gpu_Gtem; 
    cuDoubleComplex *gpu_Gii; 
    cuDoubleComplex *gpu_temp; 
    cuDoubleComplex *gpu_Hii; 
    cuDoubleComplex *gpu_Imatrix; 
    int *gpu_ipiv;

    // GPU storage space for wavefunctions
    REAL *gpu_states;

    // GPU temporary storage space for weights
    REAL *gpu_weight;
    REAL *gpu_Bweight;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    REAL *gpu_work1;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    REAL *gpu_work2;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    REAL *gpu_work3;

    // GPU work space. Sized as sbasis*THREADS_PER_NODE
    REAL *gpu_work4;

    // GPU storage space for matrix dimensioned (ct.num_states, ct.num_states)
    REAL *gpu_global_matrix;

    // Pinned host memory for finite difference routines. Allocation is slow so it
    // needs to be done once at initializatio time for each thread.
    REAL *gpu_host_temp1;
    REAL *gpu_host_temp2;
    REAL *gpu_host_temp3;
    REAL *gpu_host_temp4;
    REAL *gpu_host_fdbuf1;
    REAL *gpu_host_fdbuf2;
    REAL *gpu_host_work;

#endif


} CONTROL;


/* Extern declaration for the main control structure */
extern CONTROL ct;
STATE states[MAX_STATES];
STATE states1[MAX_STATES];
STATE states_res[MAX_STATES];
STATE states_res1[MAX_STATES];
STATE states_tem[MAX_STATES];

/* Extern declaration for the processor control structure */
extern PE_CONTROL pct;



/* Thread control structures */
typedef struct
{

    /* Thread ID number assigned by us */
    int tid;

    /* MPI communicator for use by this thread */
    MPI_Comm grid_comm;

#if GPU_ENABLED
    // Cuda device context
    cudaStream_t cstream;
    REAL *gpu_host_temp1;
    REAL *gpu_host_temp2;
#endif

    /* Thread identifier from pthread_self. Needed to send signals */
    pthread_t pthread_tid;

    /* Assigned job */
    int job;

    /* Synchronization semaphore */
    sem_t sync;

    /* These volatiles are used as synchronization variables for the threads */
    volatile int start;

    /* With the complex option this lets the threads know which k-point is
     * currently being worked on in ortho and subdiag. */
    int kidx;

    /* Pointer to current state assigned to the thread when used in sections that process a single state */
    STATE *sp;

    /* Pointer to state array used by each thread */
    STATE *my_states;

    /* Local variable -- summed to obtain total charge for all orbitals */
    REAL tcharge;

    /* Spacial offset for the thread */
    int offset;

    /* Points to base of distributed storage array for this thread */
    REAL *base_mem;

    /* Points to base of distributed scratch array for this thread */
    REAL *scratch1;

    /* Number of points per wavefunction in the distributed storage array */
    int numpt;

    /* leading dimension of the distributed wave function storage array */
    int lda;

    /* Local copies of eigenvalues and occupations for this thread */
    REAL *eigs;
    REAL *occs;

    /* Force contributions computed by this thread */
    REAL force[MAX_IONS][3];

    /* Pointer to dynamically allocated arrays of size ct.num_states*ct.num_states */
    /* that is required in ortho. Each thread has it's own copy */
    REAL *darr;
    REAL *barr;


    /* The same array as referenced by darr but this copy is 
     *allocated in the main program rather than in one of the threads.
     */
    REAL *farr;


    REAL *rho;
    REAL *rhocore;
    REAL *vtot;
    REAL *vnuc;

    /* Pointers to the non-local potential index list 
     *and to the projectors themselves */
    int *nlindex;
    REAL *projectors;

    // Pointers to special args
    void *p1;
    void *p2;
    REAL *trade_buf;// Used by trade_images
    int ion;        // Used for threaded beta_xpsi
    int nion;       // Used for threaded beta_xpsi
    REAL *sintR;    // Used for threaded beta_xpsi
    REAL *sintI;    // Used for threaded beta_xpsi
    REAL *weiptr;   // Used for threaded beta_xpsi
    int kpt;    // Used for threaded beta_xpsi
} SCF_THREAD_CONTROL;

/* Extern declarations for thread control structures */
extern SCF_THREAD_CONTROL thread_control[];

