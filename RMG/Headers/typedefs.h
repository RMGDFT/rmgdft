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

#ifndef TYPEDEFS_H
#define TYPEDEFS_H 1

#if GPU_ENABLED
#include <cuda.h>
#include <cublas_v2.h>
#endif

#include "fftw.h"
#include "mpi.h"
#include "my_scalapack.h"

#if MPI
typedef struct
{
    int count;
} QMD_sem_t;
#endif


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
    //rmg_double_t *weight[MAX_IONS];
    rmg_double_t *weight;
    rmg_double_t *Bweight;

#if FDIFF_BETA
    /*These are used for non-local force */
    //rmg_double_t *weight_derx[MAX_IONS];
    rmg_double_t **weight_derx;
    //rmg_double_t *weight_dery[MAX_IONS];
    rmg_double_t **weight_dery;
    //rmg_double_t *weight_derz[MAX_IONS];
    rmg_double_t **weight_derz;
#endif


    /** An index array which maps the projectors onto the 3-d grid associated
        with each processor.
    */
    //int *nlindex[MAX_IONS];
    int **nlindex;
    //int *Qindex[MAX_IONS];
    int **Qindex;

    /** An index array which indicate whether the grid map on the current pocessor*/
    //int *idxflag[MAX_IONS];
    int **idxflag;
    //int *Qdvec[MAX_IONS];
    int **Qdvec;

    /** Number of points in the nlindex array for each ion */
    //int idxptrlen[MAX_IONS];
    int *idxptrlen;
    //int Qidxptrlen[MAX_IONS];
    int *Qidxptrlen;

    /** Number of points in the circle of local projector for each pocessor*/
    //int lptrlen[MAX_IONS];
    int *lptrlen;

    /** Phase shifts for the non-local operators */
    //rmg_double_t *phaseptr[MAX_IONS];
    rmg_double_t **phaseptr;

    /** Points to start of storage for theaugument function*/
    //rmg_double_t *augfunc[MAX_IONS];
    rmg_double_t **augfunc;

    /** points to start of DnmI function storage for this ion*/
    //rmg_double_t *dnmI[MAX_IONS];
    rmg_double_t **dnmI;

    /** points to start of qqq storage for this ion*/
    //rmg_double_t *qqq[MAX_IONS];
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

} PE_CONTROL;



/**@name STATE
 *
 * @memo Wavefunction storage structure */
typedef struct
{

    /** First iteration flag */
    int firstflag;

    /** Current estimate of the eigenvalue for this orbital (state). */
    rmg_double_t eig[2];

    /** Previous estimate */
    rmg_double_t oldeig[2];

    /** Wavefunction residual error computed by multigrid solver */
    rmg_double_t res;

    /** Points to the storage area for the real part of the orbital */
    rmg_double_t *psiR;
    /** Points to the storage area for the imaginary part of the orbital */
    rmg_double_t *psiI;


    /** Nuclear potential */
    rmg_double_t *vnuc;
    /** Hartree potential */
    rmg_double_t *vh;
    /** Exchange correlation potential */
    rmg_double_t *vxc;
    /** Total potential */
    rmg_double_t *vtot;

    /** dvhxc */
    rmg_double_t *dvhxc;

    /** Core charge for non-linear core corrections */
    rmg_double_t *rhocore;

    /** Grid dimension in the x-coordinate direction on this processor */
    int dimx;
    /** Grid dimension in the y-coordinate direction on this processor */
    int dimy;
    /** Grid dimension in the z-coordinate direction on this processor */
    int dimz;


    /** Grid spacings */
    rmg_double_t hxgrid;
    rmg_double_t hygrid;
    rmg_double_t hzgrid;


    /** Total basis size on each processor (dimx*dimy*dimz) */
    int pbasis;

    /* Total basis size in a smoothing grid on each processor (dimx+2)*(dimy+2)*(dimz+2) */
    int sbasis;

    /*8 Index of the orbital */
    int istate;


    /** Volume element associated with each real space grid point */
    rmg_double_t vel;


    /** Occupation of the orbital */
    rmg_double_t occupation[2];
//    rmg_double_t oldeig;

    /* The offsets and the sizes of the grid that the orbital
     * is defined on relative to the global grid. These will
     * be used in the future for cluster boundary condition or
     * localized orbitals in an Order(N) formulation.
     */
    int xoff, yoff, zoff;
    int xsize, ysize, zsize;

    /** Index showing which k-point this orbital is associated with */
    int kidx;


} STATE;


/**@name SPECIES
 * @memo Species (pseudopotential) control structure
 * @doc Structure holds data about the pseudopotentials used to
 * represent different types of atomic species. 
*/
typedef struct
{

    /* symbol read from control file */
    char pseudo_symbol[32];

    /* pseudopotential filename */
    char pseudo_filename[MAX_PATH];

    /** Description of the species (e.g Atomic carbon generated using 
     * hamann's code with  rcs=0.80 rcp=0.85 bohr
     */
    char description[MAX_CHAR];

    /** Atomic number */
    int atomic_number;

    /** Atomic symbol */
    char *atomic_symbol;

    /** Atomic mass */
    rmg_double_t atomic_mass;

    /** Number of valence electrons */
    rmg_double_t zvalence;

    /** Gaussian charge parameter used for compensating the 1/r Coulomb
     * tail of the pseudopotentials
     */

    rmg_double_t rc;

    /* Number of grid points in the local in each coordinate direction. 
     * These used to be L0_LDIM and L0_NLDIM.
     */
    int ldim;
    int nldim;
    int nlfdim;
    int qdim;




    /* These are input parameters in the pseudopotential file. They represent the
     * real radii that are used in generating ldim and nldim.
     */
    rmg_double_t lradius;
    rmg_double_t nlradius;
    rmg_double_t qradius;

    /*Radius for milliken analysis*/
    rmg_double_t mill_radius;
    /*Radius in number of grid points*/
    int mill_dim;
    /*Number of radial atomic wave functions - these depend on l only, not on m*/
    int num_atomic_waves;
    /*l-numbers for states for which we have atomic orbitals*/
    int atomic_wave_l[5];
    rmg_double_t atomic_wave_oc[5];
    
    char atomic_wave_label[5][3];

    rmg_double_t *atomic_rho;
    
    /* Pseudo atomic valence density read from PP file in log grid*/
    rmg_double_t **atomic_wave;
    /* Pseudo atomic valence density on linear grid*/
    rmg_double_t **awave_lig;


    /*Sum of all atomic states (with different l or m numbers*/
    //int sum_atomic_waves;

    /*This will store name of atomic wavefunctions, such as s, px, py, pz, dxx etc*/
    //char atomic_wave_symbol[20][12];


    /** Number of radial grid points in the pseudopotential file */
    int rg_points;

    /* Log mesh parameter, where aa=exp(-aasf)/Z, bb=1.0/bbsf */
    rmg_double_t aa, bb;

    /** Non-linear core correction flag */
    int nlccflag;

    /* Number of potentials */
    int num_potentials;

    /* L-values for the reference states */
    int lval[10];

    /* L-value for local pseudopotential state */
    int local;

    /*Number of grid points in the beta function */
    int kkbeta;

    /*matrix ddd0(nbeta,nbeta) */
    rmg_double_t ddd0[18][18];
    rmg_double_t ddd[18][18];

    /*matrix qqq(nbeta,nbeta) */
    rmg_double_t qqq[18][18];

    /*the number of L=|l1-l2|.....|l1+l2|, we limit nlc <=5 */
    int nlc;

    /*the number of component in polynomial of the pseudized Q_I(r) function we limit nqf<=10 */
    int nqf;

    /*L-independent inner coutoff radii rinner for Q_I(r) function */
    rmg_double_t rinner[5];

    /* ultrosoft Vanderbilt Qnm_rad(r) function and */
    rmg_double_t *qnm;
    rmg_double_t *qnmlig;
    rmg_double_t *drqnmlig;

    /* the coefficient for pseudosation of Qnm_L(r) */
    rmg_double_t *qfcoef;

    /* Logarithmic radial mesh information */
    rmg_double_t r[MAX_RGRID];
    rmg_double_t rab[MAX_RGRID];


    /* Local Pseudopotentials */
    rmg_double_t vloc0[MAX_RGRID];

    /* Core charge radial grids */
    rmg_double_t cr[MAX_RGRID];



    /* Pseudo atomic core density */
    rmg_double_t rspsco[MAX_RGRID];

    /*the L-value for the beta function */
    int llbeta[MAX_NB];

    /*utrosoft Vanderbilt beta_n(r) function on radial grid */
    rmg_double_t beta[MAX_NB][MAX_RGRID];


    /* Total number of projectors */
    int nbeta;


    /* Linear interpolation storage for the compensated local potential
     * and for it's radial derivative.
     */
    rmg_double_t localig[MAX_LOCAL_LIG];
    rmg_double_t drlocalig[MAX_LOCAL_LIG];

    /* Linear interpolation storage for the core charge density */
    rmg_double_t rhocorelig[MAX_LOCAL_LIG];

    /* Utrosoft Vandbelit Projectors on linear interpolation grid */
    rmg_double_t betalig[MAX_NB][MAX_LOCAL_LIG];

    /* Radial derivatives of the Utrosoft Vandbelit Projectors on linear interpolation grid */
    rmg_double_t drbetalig[MAX_NB][MAX_LOCAL_LIG];

    /* Local potential linear interpolation grid spacing */
    rmg_double_t drlig;

    /* Non-local linear interpolation grid spacing */
    rmg_double_t drnlig;

    /* Qfunction linear interpolation grid spacing */
    rmg_double_t drqlig;

    /*Grid spacing for atomic charge density on linear grid*/
    rmg_double_t drlig_arho;
    
    /*Grid spacing for atomic wave functions on linear grid*/
    rmg_double_t drlig_awave;


    /* Pseudopotential filtering parameters */
    rmg_double_t lrcut;                 /* Real space local cutoff */
    rmg_double_t nlrcut[4];             /*Real space nonlocal cutoff */
    rmg_double_t rwidth;                /* Real-space width parameter */
    rmg_double_t gwidth;                /* G-space width parameter */

    /*Filtering parameters for atomic wavefunctions and charge density*/
    rmg_double_t acut; 
    rmg_double_t aradius; 
    rmg_double_t agwidth;
    rmg_double_t arwidth;

    /* radius of atomic wavefunctions and charge in terms of number of grid points*/
    int adim_rho;
    int adim_wave;


    /*Total number (of what exactly ???) */
    int num_projectors;


    /*This will store results of forward fourier transform on the coarse grid */
    fftw_complex *forward_beta;

#if !FDIFF_BETA
    /*This will store results of forward fourier transform for derivatives of beta on the coarse grid */
    fftw_complex *forward_derbeta_x;
    fftw_complex *forward_derbeta_y;
    fftw_complex *forward_derbeta_z;
#endif

    /*Backwards wisdom for fftw */
    char *backward_wisdom;

    /*Some parameters for Q function*/
    int indv[18];
    int nhtol[18];
    int nhtom[18];
    int nh;

    /*Atomic charge density on linear grid*/
    rmg_double_t arho_lig[MAX_LOCAL_LIG];

} SPECIES;


/* Structure for storing species information for internal pseudopotentials */
typedef struct
{
    char name[4];
    rmg_double_t valence;
    rmg_double_t mass;
    rmg_double_t rc;
    int nlccflag;
    int maxl;
    int local;
} ISPECIES;



/* multigrid-parameter structure */
typedef struct
{

    /* number of global-grid pre/post smoothings and timestep */
    rmg_double_t gl_step;
    int gl_pre;
    int gl_pst;

    /* timestep for the subiteration */
    rmg_double_t sb_step;

    /* timestep for the Richardson-Iteration */
    rmg_double_t ri_step;

    /* lowest MG level */
    int levels;

    /* Number of Mu-cycles to use */
    int mucycles;

    /* Number of Smoother iterations on the coarsest level */
    int coarsest_steps;

} MG_PARM;

/* Nose control structure */
typedef struct
{

    /* number of atoms allowed to move */
    int N;

    /* ionic target temperature in Kelvin */
    rmg_double_t temp;

    /* ionic target kinetic energy */
    rmg_double_t k0;

    /* randomize velocity flag */
    int randomvel;

    /* Nose oscillation frequency */
    rmg_double_t fNose;

    /* number of thermostats used */
    int m;

    /* thermostat positions,velocities,masses and forces */
    rmg_double_t xx[10];
    rmg_double_t xv[10];
    rmg_double_t xq[10];
    rmg_double_t xf[4][10];

} FINITE_T_PARM;


/** @name KPOINT
 * @memo Holds data specific to individual k-points.
 */
typedef struct
{

    /** The index of the k-point for backreferencing */
    int kidx;

    /** The k-point */
    rmg_double_t kpt[3];

    /** The corresponding vector */
    rmg_double_t kvec[3];

    /** The weight associated with the k-point */
    rmg_double_t kweight;

    /** The magnitude of the k-vector */
    rmg_double_t kmag;

    /* The orbital structure for this k-point */
    STATE *kstate;


    /* Mean min, and max wavefunction residuals for occupied space */
    rmg_double_t meanres;
    rmg_double_t minres;
    rmg_double_t maxres;

    /* Total energies */
    rmg_double_t ES;
    rmg_double_t NUC;
    rmg_double_t KE;
    rmg_double_t XC;
    rmg_double_t NL;
    rmg_double_t II;
    rmg_double_t TOTAL;

} KPOINT;




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

    /** The fine uniform grid spacing in x */
    rmg_double_t hxxgrid;

    /** The fine uniform grid spacing in y */
    rmg_double_t hyygrid;

    /** The fine uniform grid spacing in z */
    rmg_double_t hzzgrid;

    /* Kohn-sham finite difference order */
    int kohn_sham_fd_order;

    /** Lattice information */
    rmg_double_t celldm[6];

    /* lattice vectors */
    rmg_double_t a0[3];
    rmg_double_t a1[3];
    rmg_double_t a2[3];

    /** Total cell volume */
    rmg_double_t omega;

    /* This is the max of nldim for any species cubed */
    int max_nlpoints;
    int max_nlfpoints;
    int max_Qpoints;

    /** Maximum grid spacing in any coordinate direction */
    rmg_double_t hmaxgrid;


    /** Minimum grid spacing in any coordinate direction */
    rmg_double_t hmingrid;


    /** Volume element associated with each grid point */
    rmg_double_t vel;
    rmg_double_t vel_f;


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

    /* Number of CPU's in system */
    int ncpus;

    /* Potential acceleration constant step factor */
    rmg_double_t potential_acceleration_constant_step;

    /* Potential acceleration constant step factor */
    rmg_double_t potential_acceleration_poisson_step;

#if PAPI_PERFMON

    // Global event set for serial portion of code
    int EventSet;

    // Event set for Open MP threads
    int OpenMpEventSet[MAX_SCF_THREADS];

    // Flop counts for OpenMp threads
    long long OpenMpFlopCount[MAX_SCF_THREADS];

    pthread_t OpenMpPthreadId[MAX_SCF_THREADS];

    // Flop counts for pthreads
    long long PthreadFlopCount[MAX_SCF_THREADS];
    
#endif

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

} CONTROL;


/* Extern declaration for the main control structure */
extern CONTROL ct;


/* Extern declaration for the processor control structure */
extern PE_CONTROL pct;

#endif
