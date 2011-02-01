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

/*

                        main.h

    Header file for molecular dynamics code.



*/

#define REAL_SPACE 1


/* stdio for file handle argument type */
#include <stdio.h>
#include <stdlib.h>

/* Use predefined booleans */
#include <stdbool.h>

/* Size of floating point variables used in QMD */
#define     REAL    double

/* Version information */

#include    "version.h"

/* Compile time parameters */
#include    "params.h"

/* Constants and symbolic definitions */
#include    "const.h"

/* include the mpi wrapper */
#include    "my_mpi.h"

/* fourier transformation structure definitions*/
#include    "fftw.h"

#ifdef SMP
#  include <pthread.h>
#  include <semaphore.h>

typedef struct
{
    volatile int count;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} QMD_thread_barrier_struct;

void QMD_thread_barrier (QMD_thread_barrier_struct * bs);

#  if (AIX || LINUX || XT3)
typedef struct
{
    pthread_mutex_t lock;
    pthread_cond_t cond;
    int count;
} QMD_sem_t;
#  else
typedef struct
{
    sem_t sem;
} QMD_sem_t;
#  endif
#endif


#if MPI
typedef struct
{
    int count;
} QMD_sem_t;
#endif

/* Processor grid storage on finest level */
typedef struct
{

    REAL b[PX0_GRID][PY0_GRID][PZ0_GRID];

} P0_GRID_S;


/* Here we use a union of P0_GRID_S and a single long real array so that */
/* we can access things in that manner.                              */
typedef union
{

    P0_GRID_S s1;
    REAL s2[PX0_GRID * PY0_GRID * PZ0_GRID];

} P0_GRID;

/* Processor grid storage on finest level on the Fine Grid*/
typedef struct
{

    REAL b[FPX0_GRID][FPY0_GRID][FPZ0_GRID];

} FP0_GRID_S;


/* Here we use a union of FP0_GRID_S and a single long real array so that */
/* we can access things in that manner.                              */
typedef union
{
    FP0_GRID_S s1;
    REAL s2[FPX0_GRID * FPY0_GRID * FPZ0_GRID];

} FP0_GRID;


/* Smoothing grid storage on finest level */
typedef struct
{

    REAL b[PX0_GRID + 2][PY0_GRID + 2][PZ0_GRID + 2];

} S0_GRID_S;


/* Here we use a union of S0_GRID_S and a single long real array so that */
/* we can access things in that manner.                              */
typedef union
{

    S0_GRID_S s1;
    REAL s2[(PX0_GRID + 2) * (PY0_GRID + 2) * (PZ0_GRID + 2)];

} S0_GRID;

/* Smoothing grid storage on finest level on the Fine Grid*/
typedef struct
{

    REAL b[FPX0_GRID + 2][FPY0_GRID + 2][FPZ0_GRID + 2];

} FS0_GRID_S;

/* Here we use a union of FS0_GRID_S and a single long real array so that */
/* we can access things in that manner.                              */
typedef union
{

    FS0_GRID_S s1;
    REAL s2[(FPX0_GRID + 2) * (FPY0_GRID + 2) * (FPZ0_GRID + 2)];

} FS0_GRID;


/* For applying higher order finite difference operators */
typedef struct
{

    REAL b[PX0_GRID + 4][PY0_GRID + 4][PZ0_GRID + 4];

} SS0_GRID;

typedef struct
{

    REAL b[PX0_GRID + 6][PY0_GRID + 6][PZ0_GRID + 6];

} S30_GRID;

typedef struct
{

    REAL b[PX0_GRID + 10][PY0_GRID + 10][PZ0_GRID + 10];

} S50_GRID;


typedef struct
{

    REAL b[FPX0_GRID + 4][FPY0_GRID + 4][FPZ0_GRID + 4];

} FSS0_GRID;


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
    int thispe, imgpe, thisimg, thisspin, thisgrid;

	/** Number of grids (typically 1) per image to be run simultaneously **/
	int images, grids;

	/* MPI communicators for each code grid (grid_comm) and one (rmg_comm)
	 * for all group rank 0 pe's. The later effectively replaces MPI_COMM_WORLD
	 * unless you really need all-to-all, even across grids, communication. */
	MPI_Comm rmg_comm, img_topo_comm, grid_topo_comm, grid_comm, img_comm, spin_comm;

    /* determine if this image is processing spin up or spin down. */
    int spin_flag;

    /* determine whether to initialize up and down density equally or not */
    int init_equal_density_flag;

    /*Whether pe participates in scalapack calculations*/
    int scalapack_pe;

    /* Row that given processor handles when using scalapack*/
    int scalapack_myrow;
    
    /* Column that given processor handles when using scalapack*/
    int scalapack_mycol;

    /*Processor distribution for scalapack*/
    int scalapack_nprow;
    int scalapack_npcol;

    /** Neighboring processors in three-dimensional space */
    int neighbors[6];


    /** Processor x-coordinate for domain decomposition */
    int pe_x;
    /** Processor y-coordinate for domain decomposition */
    int pe_y;
    /** Processor z-coordinate for domain decomposition */
    int pe_z;


    /** Points to start of projector storage for this ion in projector space */
    //REAL *weight[MAX_IONS];
    REAL **weight;

#if FDIFF_BETA
    /*These are used for non-local force */
    //REAL *weight_derx[MAX_IONS];
    REAL **weight_derx;
    //REAL *weight_dery[MAX_IONS];
    REAL **weight_dery;
    //REAL *weight_derz[MAX_IONS];
    REAL **weight_derz;
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
    //REAL *phaseptr[MAX_IONS];
    REAL **phaseptr;

    /** Number of projectors associated with each ion. */
    //int prj_per_ion[MAX_IONS];
    int *prj_per_ion;

    /** Points to start of storage for theaugument function*/
    //REAL *augfunc[MAX_IONS];
    REAL **augfunc;

    /** points to start of DnmI function storage for this ion*/
    //REAL *dnmI[MAX_IONS];
    REAL **dnmI;

    /** points to start of qqq storage for this ion*/
    //REAL *qqq[MAX_IONS];
    REAL **qqq;
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


    /** Nuclear potential */
    REAL *vnuc;
    /** Hartree potential */
    REAL *vh;
    /** Exchange correlation potential */
    REAL *vxc;
    /** Total potential */
    REAL *vtot;

    /** Core charge for non-linear core corrections */
    REAL *rhocore;

    /** Grid dimension in the x-coordinate direction on this processor */
    int dimx;
    /** Grid dimension in the y-coordinate direction on this processor */
    int dimy;
    /** Grid dimension in the z-coordinate direction on this processor */
    int dimz;


    /** Grid spacings */
    REAL hxgrid;
    REAL hygrid;
    REAL hzgrid;


    /** Total basis size on each processor (dimx*dimy*dimz) */
    int pbasis;

    /* Total basis size in a smoothing grid on each processor (dimx+2)*(dimy+2)*(dimz+2) */
    int sbasis;

    /*8 Index of the orbital */
    int istate;


    /** Volume element associated with each real space grid point */
    REAL vel;


    /** Occupation of the orbital */
    REAL occupation;

    REAL oldeig;

    /* The offsets and the sizes of the grid that the orbital
     * is defined on relative to the global grid. These will
     * be used in the future for cluster boundary condition or
     * localized orbitals in an Order(N) formulation.
     */
    int xoff, yoff, zoff;
    int xsize, ysize, zsize;

    /** Index showing which k-point this orbital is associated with */
    int kidx;

    /* eigenvalue of the opposite spin for the corresponding orbital (state)*/
    REAL eig_oppo;

    /* Hold occupation of the opposite spins orbital*/
    REAL occupation_oppo;

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
    REAL atomic_mass;

    /** Number of valence electrons */
    REAL zvalence;

    /** Gaussian charge parameter used for compensating the 1/r Coulomb
     * tail of the pseudopotentials
     */

    REAL rc;

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
    REAL lradius;
    REAL nlradius;
    REAL qradius;

    /*Radius for milliken analysis*/
    REAL mill_radius;
    /*Radius in number of grid points*/
    int mill_dim;
    /*Number of radial atomic wave functions - these depend on l only, not on m*/
    int num_atomic_waves;
    /*l-numbers for states for which we have atomic orbitals*/
    int lstate_atomic_wave[5];
    /*Sum of all atomic states (with different l or m numbers*/
    int sum_atomic_waves;

    /*This will store name of atomic wavefunctions, such as s, px, py, pz, dxx etc*/
    char atomic_wave_symbol[20][12];


    /** Number of radial grid points in the pseudopotential file */
    int rg_points;

    /* Log mesh parameter, where aa=exp(-aasf)/Z, bb=1.0/bbsf */
    REAL aa, bb;

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

    /*Number of grid points in the beta function */
    int kkbeta;

    /*matrix ddd0(nbeta,nbeta) */
    REAL ddd0[18][18];
    REAL ddd[18][18];

    /*matrix qqq(nbeta,nbeta) */
    REAL qqq[18][18];

    /*the number of L=|l1-l2|.....|l1+l2|, we limit nlc <=5 */
    int nlc;

    /*the number of component in polynomial of the pseudized Q_I(r) function we limit nqf<=10 */
    int nqf;

    /*L-independent inner coutoff radii rinner for Q_I(r) function */
    REAL rinner[5];

    /* ultrosoft Vanderbilt Qnm_rad(r) function and */
    REAL *qnm;
    REAL *qnmlig;
    REAL *drqnmlig;

    /* the coefficient for pseudosation of Qnm_L(r) */
    REAL *qfcoef;

    /* Logarithmic radial mesh information */
    REAL r[MAX_RGRID];
    REAL rab[MAX_RGRID];


    /* Local Pseudopotentials */
    REAL vloc0[MAX_RGRID];

    /* Core charge radial grids */
    REAL cr[MAX_RGRID];


    /* Pseudo atomic valence density */
    /*REAL avdens[MAX_RGRID];*/
    REAL **atomic_wave;


    /* Pseudo atomic core density */
    REAL rspsco[MAX_RGRID];

    /*the L-value for the beta function */
    int llbeta[MAX_NB];

    /*utrosoft Vanderbilt beta_n(r) function on radial grid */
    REAL beta[MAX_NB][MAX_RGRID];


    /* Total number of projectors */
    int nbeta;


    /* Linear interpolation storage for the compensated local potential
     * and for it's radial derivative.
     */
    REAL localig[MAX_LOCAL_LIG];
    REAL drlocalig[MAX_LOCAL_LIG];

    /* Linear interpolation storage for the core charge density */
    REAL rhocorelig[MAX_LOCAL_LIG];

    /* Utrosoft Vandbelit Projectors on linear interpolation grid */
    REAL betalig[MAX_NB][MAX_LOCAL_LIG];

    /* Radial derivatives of the Utrosoft Vandbelit Projectors on linear interpolation grid */
    REAL drbetalig[MAX_NB][MAX_LOCAL_LIG];

    /* Local potential linear interpolation grid spacing */
    REAL drlig;

    /* Non-local linear interpolation grid spacing */
    REAL drnlig;

    /* Qfunction linear interpolation grid spacing */
    REAL drqlig;


    /* Pseudopotential filtering parameters */
    REAL lrcut;                 /* Real space local cutoff */
    REAL nlrcut[4];             /*Real space nonlocal cutoff */
    REAL rwidth;                /* Real-space width parameter */
    REAL gwidth;                /* G-space width parameter */


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


} SPECIES;


/* Structure for storing species information for internal pseudopotentials */
typedef struct
{
    char name[4];
    REAL valence;
    REAL mass;
    REAL rc;
    int nlccflag;
    int maxl;
    int local;
} ISPECIES;



/*Structure for storing PDB information
 * Each ion should have it*/
typedef struct
{

/* 1 -  6  Record name*/
char record_name[7];

/* 7 - 11 Atom serial number*/
int serial_num;

/*13 - 16  Atom name*/
char name[5];

/* 17 Alternate location indicator.*/
char altLoc[2];

/* 18 - 20 Residue name*/
char resName[4];

/* 22 Chain identifier*/
char chainID[2];

/* 23 - 26 Residue sequence number*/
int resSeq;

/* 27 Code for insertion of residues*/
char iCode[2];

/* 55 - 60 Occupancy*/
REAL occupancy;

/* 61 - 66 Temperature factor*/
REAL tempFactor;

/* 77 - 78  Element symbol, right-justified. */
char element[3];

/*79 - 80  Charge on the atom.*/
char charge[3];

} PDB_INFO;




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

    /*Position of ion relative to the middle of non-local box around the ion 
     *          * determined in get_nlop, AIget_cindex sets this up*/
    REAL nlcrds[3];


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


    /* Coordinates of the corner of the grid that the Qfunction */
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

	/* Movement constraint: float[3] vector, int constraint_type (0=disabled, 1=in-plane, 2=along vector) */
	REAL constraint[3];
	int constraint_type;

    /* Stored non-local projectors */
    REAL *oldsintR;
    REAL *oldsintI;
    REAL *newsintR;
    REAL *newsintI;

    /* Stores sine and cosine of a phase factor for backwards fourier transform */
    REAL *fftw_phase_sin;
    REAL *fftw_phase_cos;


    /*Stores PDB information*/
    PDB_INFO pdb;


} ION;



/* multigrid-parameter structure */
typedef struct
{

    /* number of global-grid pre/post smoothings and timestep */
    REAL gl_step;
    int gl_pre;
    int gl_pst;

    /* timestep for the subiteration */
    REAL sb_step;

    /* timestep for the Richardson-Iteration */
    REAL ri_step;

    /* lowest MG level */
    int levels;


} MG_PARM;

/* Nose control structure */
typedef struct
{

    /* number of atoms allowed to move */
    int N;

    /* ionic target temperature in Kelvin */
    REAL temp;

    /* ionic target kinetic energy */
    REAL k0;

    /* randomize velocity flag */
    int randomvel;

    /* Nose oscillation frequency */
    REAL fNose;

    /* number of thermostats used */
    int m;

    /* thermostat positions,velocities,masses and forces */
    REAL xx[10];
    REAL xv[10];
    REAL xq[10];
    REAL xf[4][10];

} FINITE_T_PARM;


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

    /* The orbital structure for this k-point */
    STATE *kstate;


    /* Mean min, and max wavefunction residuals for occupied space */
    REAL meanres;
    REAL minres;
    REAL maxres;

    /* Total energies */
    REAL ES;
    REAL NUC;
    REAL KE;
    REAL XC;
    REAL NL;
    REAL II;
    REAL TOTAL;

} KPOINT;




/** @name CONTROL
  @memo Main control structure
 
  This is a global structure declared as extern CONTROL ct
 
 */
typedef struct
{

#ifdef SMP
    /* Number of threads to run concurrently. We always 
     * use one thread per eigenfunction but it's inefficient 
     * to have more threads than the number of CPU's available
     * running concurrently. Defaults to compile time parameter 
     * of 1 but can be overridden by setting the environment 
     * variable QMD_NUM_THREADS to any value up to MAX_THREADS.
     */
    int thread_concurrency;
    int num_threads;

    /* When running on a cluster this is the total number of nodes
     * in the cluster.
     */
    int num_nodes;

    /* While this is the identity of this node */
    int this_node;

    /* And this is the space offset of this node */
    int node_space_offset;

    /* And this is the number of space points handled by this node */
    int node_space_size;
#endif

    /** Description of the run. */
    char description[MAX_CHAR];

    /* time at which run started */
    REAL time0;

    /** Name of the input control file. Passed as a command line argument
     *
     *  Example:
     *  bash$  md in.diamond8
     */
    char cfile[MAX_PATH];

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

    /* Maximum number of fast relax steps */
    int max_rlx_steps;

    /* Maximum number of rmg meta loops (NEB, ARTS, etc.) */
    int max_rmg_steps;

    /* MD steps iterator */
    int md_steps;

    /* Maximum number of SCF steps in a MD steps */
    int max_scf_steps;

    /* Total number of SCF steps done */
    int total_scf_steps;

    /* SCF steps iterator */
    int scf_steps;

    /* override occupations */
    int override_occ;

    /* Override current ionic positions (in control file) with positions from wave file (during restart) */
    int override_current;

    /* Override initial ionic positions (in control file) with positions from wave file (during restart) */
    int override_initial;

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

    /* Number of states for the opposite spin*/
    int num_states_oppo; 

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

    /** the fine grid size on each coarse grid cell */
    int nxfgrid;
    int nyfgrid;
    int nzfgrid;

    /** Global uniform grid spacing in x */
    REAL hxgrid;

    /** Global uniform grid spacing in y */
    REAL hygrid;

    /** Global uniform grid spacing in z */
    REAL hzgrid;

    /** The fine uniform grid spacing in x */
    REAL hxxgrid;

    /** The fine uniform grid spacing in y */
    REAL hyygrid;

    /** The fine uniform grid spacing in z */
    REAL hzzgrid;

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
    int max_nlfpoints;
    int max_Qpoints;

    /** Maximum grid spacing in any coordinate direction */
    REAL hmaxgrid;


    /** Minimum grid spacing in any coordinate direction */
    REAL hmingrid;


    /** Grid anisotropy defined as the ratio of hmaxgrid to hmingrid. A value larger than 1.05 can lead to convergence problems. */
    REAL anisotropy;


    /** Volume element associated with each grid point */
    REAL vel;
    REAL vel_f;


    /** Physical grid basis size */
    int nbasis;


    /** Density mixing parameter. Typical values range from 0.2 to 0.9, while larger values provide faster convergence as long as they are stable. */
    REAL mix;


    /* Projector mixing parameter */
    REAL prjmix;

    /* Global uniform grid corner */
    REAL xcstart;
    REAL ycstart;
    REAL zcstart;


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


    /** Nose paramters */
    FINITE_T_PARM nose;

    /* force pointer array */
    int fpt[4];

    /* temperature control */
    int tcontrol;

    /* md integrate level */
    int mdorder;

    /* movie flags */
    int rmvmovie, chmovie, xbsmovie;

    /* Milliken population flags. */
    int domilliken;
    int milliken;

    /* Diagonalization flag and period */
    int initdiag;
    int diag;
    int end_diag;

    /* How many steps between writeout of eigenvalues*/
    int write_eigvals_period;

    /* Diagonalizations during the md step */
    int mddiag1;
    int mddiag2;

    /** Force flag. 0=don't compute forces, 1=compute forces */
    int forceflag;

    /* Whether to write full memory usage report at the end of calculation */
    int write_memory_report;

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

    /* fermi energy */
    REAL efermi;

    /** Total number of k-points being used in the calculation */
    int num_kpts;


    /** K-point control structure */
    KPOINT *kp;

    /** Throttles data transfer rates when writing wavefunctions to disk
     *
     * On clusters with NFS mounted filesystems having all nodes
     * dump there data at the same time can cause network congestion
     * and hangups so wait_flag can be set in the input file to throttle
     * the total bandwidth being written. */
    int wait_flag;

    /** The maximum number of projectors for any species */
    int max_nl;

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
    REAL e_field;

    REAL x_field_0;

    REAL y_field_0;

    REAL z_field_0;


} CONTROL;


#ifdef SMP

/* Thread control structures */
typedef struct
{

    /* Thread ID number */
    int tid;

    /* These volatiles are used as synchronization variables for the threads */
    volatile int start;

    /* With the complex option this lets the threads know which k-point is
     * currently being worked on in ortho and subdiag. */
    int kidx;

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


} SCF_THREAD_CONTROL;
#endif



/* Extern declaration for the main control structure */
extern CONTROL ct;


/* Extern declaration for the processor control structure */
extern PE_CONTROL pct; 


/* Extern declarations for thread control structures */
#ifdef SMP
extern SCF_THREAD_CONTROL thread_control[];
#endif


/* Header file for blas routines */
#include "blas.h"


/* Function prototypes */
void app_4del2 (S0_GRID *f, P0_GRID *work);
REAL app_del2c (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                REAL gridhx, REAL gridhy, REAL gridhz);
void app6_del2 (REAL *rho, P0_GRID * work);
void app6_del2f (REAL *rho, FP0_GRID * work);
void app_smooth (S0_GRID *f, S0_GRID *work, REAL sfac);
/*void app_smooth1(FS0_GRID *f, FS0_GRID *work);*/
void app_cir_sixth (REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir (REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_ortho (REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_bcc (REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_fcc (REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_hex (REAL *a, REAL *b, int dimx, int dimy, int dimz);
REAL app_cilr (REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
               REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_bcc (REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
                   REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_fcc (REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
                   REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_hex (REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
                   REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_ortho (REAL *a, REAL *b, REAL *c, int dimx, int dimy,
                     int dimz, REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil (REAL *a, REAL *b, int dimx, int dimy, int dimz, REAL gridhx,
              REAL gridhy, REAL gridhz);
REAL app_cil_sixth (REAL *psi, REAL *b, int dimx, int dimy, int dimz,
                    REAL gridhx, REAL gridhy, REAL gridhz);
void app_grad (REAL  * rho, P0_GRID * wx, P0_GRID * wy, P0_GRID * wz);
void app_gradf (REAL * rho, FP0_GRID * wx, FP0_GRID * wy, FP0_GRID * wz);
void constrain(ION *iptr);
void corlyp (REAL *dp, REAL *dm, REAL *dp1, REAL *dm1, REAL *dp2, REAL *dm2, REAL *ec,
             REAL *vcp0, REAL *vcm0, int *ndm);
void cross_product (REAL *a, REAL *b, REAL *c);
void destroy_fftw_wisdom (void);
void eval_residual (REAL *mat, REAL *f_mat, int dimx, int dimy, int dimz,
                    REAL gridhx, REAL gridhy, REAL gridhz, REAL *res);
void solv_pois (REAL *vmat, REAL *fmat, REAL *work,
                int dimx, int dimy, int dimz, REAL gridhx,
                REAL gridhy, REAL gridhz, REAL step);
REAL fill (STATE *states, REAL width, REAL nel, REAL mix,
           int num_st, int occ_flag);
REAL fill_spin(STATE *states_up, REAL * eigval_dn, REAL width, REAL nel, REAL mix, int num_st, int occ_flag);

void find_phase (int nldim, REAL *nlcdrs, REAL *phase_sin,
                 REAL *phase_cos);
void finish_release_mem(STATE *states);
void genvpsi (REAL *psi, REAL *twovpsi, REAL *pvtot, REAL *pvnl,
              REAL *kd, REAL kmag, int dimx, int dimy, int dimz);
void get_nlop (void);
void get_weight (void);
void get_phase (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_smp (int tid);
void get_eig (STATE *states, P0_GRID *vxc, P0_GRID *vh, P0_GRID *vnuc);
char *get_num (char *str);
void get_te (REAL *rho, REAL *rhocore, REAL *rhoc, REAL *vh, REAL *vxc,
             STATE *states);
void get_te_spin (REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc, REAL *vh, REAL *vxc,
             STATE *states);

void get_vxc (REAL *rho, REAL *rhocore, REAL *vxc);
void get_vxc_spin (REAL *rho, REAL * rho_oppo, REAL *rhocore, REAL *vxc);


void get_zdens (STATE *states, int state, REAL *zvec);
void xclda_pz81 (REAL *rho, REAL *vxc);
void exclda_pz81 (REAL *rho, REAL *exc);

/* new lda function incorporating both  1981 and 1994 monte carlo data */
void xclda(REAL *rho, REAL *vxc, REAL *exc);
void slater(REAL rs, REAL *ex, REAL *vx);
void pz(REAL rs, int iflag, REAL *ec, REAL *vc); 

/* lsda exchange correlation functional */
void xclsda_spin(REAL *rho, REAL *rho_oppo, REAL *vxc, REAL *exc);
void slater_spin(REAL arhox, REAL zeta, REAL *ex, REAL *vx, REAL *vx_oppo);
void pz_spin(REAL rs, REAL zeta, REAL *ec, REAL *vc, REAL *vc_oppo);
void pz_polarized(REAL rs, REAL *ec, REAL *vc);
void pw_spin (REAL rs, REAL zeta, REAL *ec, REAL *vcup, REAL *vcdw);
void pw (REAL rs, int iflag, REAL *ec, REAL *vc) ;


double mu_pz (double rho);
double e_pz (double rho);
void xcgga (REAL *rho, REAL *vxc, REAL *exc, int flag);
void xcgga_spin (REAL *rho, REAL *rho_oppo, REAL *vxc, REAL *exc, int flag);


/* exchange correlation functional for PBE */

void gcxcpbe_spin(REAL rho_up, REAL rho_dw,
  REAL grad_up, REAL grad_dw, REAL grad, REAL *enxc,
  REAL *vxc1_up, REAL *vxc1_dw, REAL *vxc2_upup, REAL *vxc2_dwdw,
  REAL *vxc2_updw, REAL *vxc2_dwup);
void pbex (REAL rho, REAL grho, int iflag, REAL *sx, REAL *v1x, REAL *v2x);
void pbec_spin (REAL rho, REAL zeta, REAL grho, int iflag, REAL *sc, REAL *v1cup, REAL *v1cdw, REAL *v2c);

void gcxcpbe (REAL rho, REAL grad, REAL *enxc, REAL *vxc1, REAL *vxc2);

void pbec (REAL rho, REAL grad, int iflag, REAL *sc, REAL *v1c, REAL *v2c);



/* becke88 exchange and perdew86 correlation */

void becke88 ( REAL rho, REAL grho, REAL *sx, REAL *v1x, REAL *v2x );
void perdew86 ( REAL rho, REAL grho, REAL *sc, REAL *v1c, REAL *v2c );
void gcxbcp (REAL rho, REAL grad, REAL *enxc, REAL *vxc1, REAL *vxc2);
void perdew86_spin ( REAL rho, REAL zeta, REAL grho, REAL *sc, REAL *v1cup, REAL *v1cdw, REAL *v2c );
void becke88_spin ( REAL rho, REAL grho, REAL *sx, REAL *v1x, REAL *v2x );
void gcxbcp_spin (REAL rho_up, REAL rho_dw,
  REAL grad_up, REAL grad_dw, REAL grad, REAL *enxc,
  REAL *vxc1_up, REAL *vxc1_dw, REAL *vxc2_upup, REAL *vxc2_dwdw,
  REAL *vxc2_updw, REAL *vxc2_dwup);


/* PW91 exchange correlation */

void ggax(REAL rho, REAL grho, REAL *sx, REAL *v1x, REAL *v2x);
void ggac (REAL rho, REAL grho, REAL *sc, REAL *v1c, REAL *v2c);
void ggac_spin (REAL rho, REAL zeta, REAL grho, REAL *sc, REAL *v1cup, REAL *v1cdw, REAL *v2c);
void gcxcpw91 (REAL rho, REAL grad, REAL *enxc, REAL *vxc1, REAL *vxc2);
void gcxcpw91_spin(REAL rho_up, REAL rho_dw,
  REAL grad_up, REAL grad_dw, REAL grad, REAL *enxc,
  REAL *vxc1_up, REAL *vxc1_dw, REAL *vxc2_upup, REAL *vxc2_dwdw,
  REAL *vxc2_updw, REAL *vxc2_dwup);


/* BLYP exchange correlation */

void lyp ( REAL rs, REAL *ec, REAL *vc );
void glyp ( REAL rho, REAL grho, REAL *sc, REAL *v1c, REAL *v2c );
void lsd_lyp ( REAL rho, REAL zeta, REAL * elyp, REAL * valyp, REAL * vblyp );
void lsd_glyp ( REAL rhoa, REAL rhob, REAL grhoaa, REAL grhoab2, REAL grhobb, REAL *sc, REAL *v1ca, REAL *v2ca, REAL *v1cb, REAL *v2cb, REAL *v2cab );
void gcxcblyp (REAL rho, REAL grad, REAL *enxc, REAL *vxc1, REAL *vxc2);
void gcxcblyp_spin (REAL rho_up, REAL rho_dw,
  REAL grad_up, REAL grad_dw, REAL grad_updw2, REAL *enxc,
  REAL *vxc1_up, REAL *vxc1_dw, REAL *vxc2_upup, REAL *vxc2_dwdw,
  REAL *vxc2_updw, REAL *vxc2_dwup);



void gram (KPOINT *kpoint, REAL h, int numst, int maxst, int numpt,
           int maxpt);
int get_input (FILE *fh, char *id, void *dest, unsigned int flag, char *def);
REAL get_ke (STATE *sp, int tid);
void get_vh (REAL *rho, REAL *rhoc, REAL *vh, int cycles, int maxlevel);
char *get_symbol (int atomic_number);
void global_sums (REAL *vect, int *length);
void global_sums_spin (REAL *vect, int *length);
void init (REAL *vh, REAL *rho, REAL *rhocore, REAL *rhoc, STATE *states,
           REAL *vnuc, REAL *vxc);
void init_spin (REAL *vh, REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc, STATE *states,
           REAL *vnuc, REAL *vxc);
void init_derweight (void);
void init_derweight_s (SPECIES *sp, fftw_complex *rtptr_x,
                       fftw_complex *rtptr_y, fftw_complex *rtptr_z, int ip,
                       fftwnd_plan p1);
void init_derweight_p (SPECIES *sp, fftw_complex *rtptr_x,
                       fftw_complex *rtptr_y, fftw_complex *rtptr_z, int ip,
                       fftwnd_plan p1);
void init_derweight_d (SPECIES *sp, fftw_complex *rtptr_x,
                       fftw_complex *rtptr_y, fftw_complex *rtptr_z, int ip,
                       fftwnd_plan p1);
void init_fftw_wisdom (void);
void init_kbr (void);
void init_IO ( int argc, char **argv );
void init_pe (void);
void init_img_topo ( int dimensionality );
void init_pegrid (void);
STATE *init_states (void);
STATE *init_states_spin(void);
void init_weight (void);
void init_weight_s (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftwnd_plan p1);
void init_weight_p (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftwnd_plan p1);
void init_weight_d (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftwnd_plan p1);
void init_wf (STATE *states);
void init_nuc (REAL *vnuc, REAL *rhoc, REAL *rhocore);
void init_pos (void);
void init_sym (void);
void symmetrize_rho (FP0_GRID *rho);
void symforce (void);
void cforce (P0_GRID *rho, P0_GRID *vh);
void mgrid_solv (REAL *v_mat, REAL *f_mat, REAL *work,
                 int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy,
                 REAL gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, REAL step);
void rmg_timings (int what, REAL time, int tid);
REAL minimage (ION *ip1, ION *ip2, REAL *xtal_r);
REAL my_crtc (void);
/* Local function prototypes */
FILE *open_xbs_movie (char *filename);
int open_wave_file (char *filename);
void xbsmovie (FILE *movie);
void ortho_half (STATE *states);
void ortho_bcc (STATE *states);
void output_eigenvalues( STATE *states, int ikbs, int iscf );
void pe2xyz (int pe, int *x, int *y, int *z);
void pack_ptos (REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void pack_stop (REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (REAL *sg, REAL *pg, REAL alpha, int dimx, int dimy,
                     int dimz);
void pack_ptos_trade (REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void pack_vhstod (REAL *s, REAL *d, int dimx, int dimy, int dimz);
void pack_vhdtos (REAL *s, REAL *d, int dimx, int dimy, int dimz);
double radint (double *f, double *r, int n, double al);
REAL radint1 (REAL *f, REAL *r, REAL *dr_di, int n);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
void read_control(void);
void write_pdb (void);
int read_atom_line(char *species, REAL *crds, int *movable, FILE *fhand, char *tbuf, int index);
int assign_species (CONTROL * c, char *buf);
void read_data (char *name, REAL *vh, REAL *rho, REAL *vxc,
                STATE *states);
void read_data_spin (char *name, REAL *vh, REAL *rho, REAL *rho_oppo, REAL *vxc,
                STATE *states);

void read_pseudo (void);
REAL real_sum_all (REAL x);
REAL real_sum_all_spin (REAL x);
REAL real_min_all (REAL x);


void reset_timers (void);
/*void scf(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
           REAL *rho, REAL *rhocore, REAL *rhoc, int *CONVERGENCE);*/
void sortpsi (STATE *states);
/*void subdiag(STATE *states, REAL *vh, REAL *vnuc, REAL *vxc);*/
void trade_images (REAL *mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_imagesx (REAL *f, REAL *w, int dimx, int dimy, int dimz,
                    int images);
void set_bc (REAL *mat, int dimx, int dimy, int dimz, int images, REAL val);
void set_bcx (REAL *mat, int dimx, int dimy, int dimz, int images, REAL val);
void getpoi_bc (REAL *rho, REAL *vh_bc, int dimx, int dimy, int dimz);
/*void trade_images2(S0_GRID *f, SS0_GRID *w);
void trade_images2f(FS0_GRID *f, FSS0_GRID *w);
void trade_images3(S0_GRID *f, S30_GRID *w);
void trade_images5(S0_GRID *f, S50_GRID *w);*/
void vol_rho (P0_GRID *rho, int step);
void vol_wf (STATE *states, int state, int step);
void write_avgd (REAL *rho);
void write_avgv (REAL *vh, REAL *vnuc);
void write_zstates (STATE *states);
void write_data (char *name, REAL *vh, REAL *rho, REAL *vxc,
                 STATE *states);

void write_data_spin (char *name, REAL *vh, REAL *rho, REAL *rho_oppo, REAL *vxc,
                 STATE *states);
void write_header (void);
void write_occ (STATE *states);
void write_force (void);
void write_timings (void);
void wvfn_residual(STATE *states);
REAL rand0 (long *idum);

void mg_restrict (REAL *full, REAL *half, int dimx, int dimy, int dimz);
void mg_prolong (REAL *full, REAL *half, int dimx, int dimy, int dimz);
void gather_psi (REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid);
void scatter_psi (REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid);
void get_milliken (STATE *states);

#ifdef SMP
void create_threads (STATE *states);
void start_threads (int action);
void wait_for_threads (void);
void thread_spinlock (int *lockptr);
void *thread_scheduler (void *);
void thread_dispatch (SCF_THREAD_CONTROL *s, int job);
void ortho1_smp (SCF_THREAD_CONTROL *s);
void ortho2_smp (SCF_THREAD_CONTROL *s);
void get_rho_smp (SCF_THREAD_CONTROL *s);
void sort_psi_smp (SCF_THREAD_CONTROL *s);
void subdiag1_smp (SCF_THREAD_CONTROL *s);
void subdiag2_smp (SCF_THREAD_CONTROL *s);
#endif

void bandstructure( STATE *states, REAL *vxc, REAL *vh, REAL *vnuc );
void output_wave( STATE *states, int kpt, int fhand );

void QMD_sem_init (QMD_sem_t *sem);
void QMD_sem_destroy (QMD_sem_t *sem);
void QMD_sem_wait (QMD_sem_t *sem);
void QMD_sem_post (QMD_sem_t *sem);

/* Blas wrappers */
void QMD_saxpy (int n, REAL alpha, REAL *x, int incx, REAL *y, int incy);
void QMD_sscal (int n, REAL alpha, REAL *x, int incx);
void QMD_scopy (int n, REAL *x, int incx, REAL *y, int incy);
REAL QMD_sdot (int n, REAL *x, int incx, REAL *y, int incy);



int get_index (ION *iptr, int *Aix, int *Aiy, int *Aiz,
               int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
               int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid,
               REAL *lxcstart, REAL *lycstart, REAL *lzcstart);

REAL linint (REAL *y, REAL rv, REAL invdr);
void my_barrier (void);


#ifdef UNICOS_T3E
#define exit globalexit
#endif

/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, REAL *celldm, REAL *A0I, REAL *A1I, REAL *A2I,
             REAL *OMEGAI, int *flag);
void recips (void);
void to_cartesian (REAL crystal[], REAL cartesian[]);
void to_crystal (REAL crystal[], REAL cartesian[]);
REAL metric (REAL *crystal);

/* Md run types */
void quench (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc, REAL *rho,
             REAL *rhocore, REAL *rhoc);
void quench_spin (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc, REAL *rho,
             REAL *rho_oppo, REAL *rhocore, REAL *rhoc);
void fastrlx (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
              REAL *rho, REAL *rhocore, REAL *rhoc);
void fastrlx_spin (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
              REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc);
void neb_relax (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
              REAL *rho, REAL *rhocore, REAL *rhoc);
void cdfastrlx (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
                REAL *rho, REAL *rhocore, REAL *rhoc);
void moldyn (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
             REAL *rho, REAL *rhoc, REAL *rhocore);
void dx (STATE *states, P0_GRID *vxc, P0_GRID *vh, P0_GRID *vnuc,
         P0_GRID *rho, P0_GRID *rhoc);
void psidx (STATE *states, P0_GRID *vxc, P0_GRID *vh, P0_GRID *vnuc,
            P0_GRID *rho, P0_GRID *rhoc);
void cholesky (REAL *a, int n);


/*the function for softpseudopotential*/
void aainit (int lli, int mix, int lx, int mx, int nlx, double ap[][9][9],
             int lpx[][9], int lpl[][9][9]);
REAL app_cil1 (REAL *a, REAL *b, int dimx, int dimy, int dimz,
               REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil1_bcc (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                   REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil1_fcc (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                   REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil1_hex (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                   REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil1_ortho (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                     REAL gridhx, REAL gridhy, REAL gridhz);
void app_nl_psi (REAL *psiR, REAL *psiI, REAL *workR, REAL *workI,
                 int state, int flag, int kidx, int tid);
void app_nl_eig (REAL *psiR, REAL *psiI, REAL *workR, REAL *workI,
                 int state, int flag, int kidx, int tid);
void app_ns_psi (REAL *psiR, REAL *psiI, REAL *workR, REAL *workI,
                 int state, int kidx, int tid);
void app_ns_eig (REAL *psiR, REAL *psiI, REAL *workR, REAL *workI,
                 int state, int kidx, int tid);
void get_ddd (REAL *veff);
void get_nlop_d (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_p (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_s (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_QI (void);
void get_qqq (void);
void get_rho (STATE * states, REAL * rho, REAL * rhocore);
void init_psp (void);
void init_qfunct (void);
void mg_eig_state (STATE *sp, int tid, REAL *vtot_psi);
void ortho_full (STATE *states);
void ortho (STATE *states, int kpt);
REAL qval (int ih, int jh, REAL r, REAL invdr, REAL *ptpr, int *nhtol,
           int *nhtom, int *indv, REAL *ylm, REAL ap[][9][9], int lpx[][9],
           int lpl[][9][9], SPECIES *sp);
void scf (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
          REAL *rho, REAL *rhocore, REAL *rhoc, int *CONVERGENCE);
void scf_spin (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
          REAL *rho,REAL *rho_oppo, REAL *rhocore, REAL *rhoc, int *CONVERGENCE);

#if GAMMA_PT
void subdiag_gamma (STATE *states, REAL *vh, REAL *vnuc, REAL *vxc);
#else
void subdiag_nongamma (STATE * states, REAL * vh, REAL * vnuc, REAL * vxc);
#endif

void ylmr2 (double *r, double *ylm);
REAL gcutoff (REAL g1, REAL gcut, REAL width);
void rft1 (REAL cparm, REAL *f, REAL *r, REAL *ffil, REAL *rab,
           int rg_points, int lval, REAL dr, REAL width, int lrg_points);
void norm_psi1 (STATE *sp, int istate, int kpt);
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL *cR, REAL *cI);
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL cR, REAL cI);
REAL get_QnmL (int idx, int ltot, REAL r, SPECIES *sp);

/*force for softpseudopotential*/
void force (REAL *rho, REAL *rhoc, REAL *vh, REAL *vxc, REAL *vnuc,
            STATE *states);

void force_spin (REAL *rho, REAL *rho_oppo, REAL *rhoc, REAL *vh, REAL *vxc, REAL *vnuc,
            STATE *states);

void iiforce (void);
void lforce (REAL *rho, REAL *vh);
void nlforce1 (REAL *veff);
void get_gamma (REAL *gammaR, ION *iptr, int nh);
void partial_gamma (int ion, REAL *par_gammaR, REAL *par_omegaR, ION *iptr,
                    int nh, REAL *newsintR_x, REAL *newsintR_y,
                    REAL *newsintR_z, REAL *newsintI_x, REAL *newsintI_y,
                    REAL *newsintI_z);
void partial_betaxpsi (int ion, fftwnd_plan p2, REAL *newsintR_x,
                       REAL *newsintR_y, REAL *newsintR_z,
                       REAL *newsintI_x, REAL *newsintI_y,
                       REAL *newsintI_z, ION *iptr);
void nlforce1_par_Q (REAL *veff, REAL *gamma, int ion, ION *iptr, int nh,
                     REAL *forces);
void nlforce1_par_gamma (REAL *par_gamma, int ion, int nh);
void nlforce1_par_omega (REAL *par_omega, int ion, ION *iptr, int nh);
void partial_QI (int ion, REAL *QI_R, ION *iptr);
void qval_R (int ih, int jh, REAL r, REAL *x, REAL *qlig, REAL *drqlig,
             REAL invdr, int *nhtol, int *nhtom, int *indv, REAL *ylm,
             REAL *ylm_x, REAL *ylm_y, REAL *ylm_z, REAL ap[][9][9],
             int lpx[][9], int lpl[][9][9], REAL *Q_x, REAL *Q_y,
             REAL *Q_z, SPECIES *sp);
void ylmr2_x (double *r, double *ylm_x);
void ylmr2_y (double *r, double *ylm_y);
void ylmr2_z (double *r, double *ylm_z);
void nlccforce (REAL *rho, REAL *vxc);
REAL get_ve_nl (STATE *sta, int istate);
void pack_rho_ctof (REAL *rhoc, REAL *rhof);
void bspline_interp_full (REAL *rho, REAL *rho_f);
void get_vtot_psi (REAL *vtot_psi, REAL *vtot);
void betaxpsi (STATE *states);
void betaxpsi1 (STATE *states, int kpt);
void assign_weight (SPECIES *sp, int ion, fftw_complex *beptr,
                    REAL *rtptr);
void assign_weight2 (int nldim, int ion, REAL *beptr, REAL *rtptr);
void pack_gftoc (SPECIES *sp, fftw_complex *gwptr, fftw_complex *gbptr);
void debug_write_rho_z (REAL *rhoz);
void print_density_z_direction (int grid_x, int grid_y, REAL *density,
                                int px0_grid, int py0_grid, int pz0_grid,
                                REAL zside);
void get_derweight (int ion, REAL *beta_x, REAL *beta_y, REAL *beta_z,
                    ION *iptr, fftwnd_plan p2);
void partial_beta_fdiff (fftw_complex *beptr, int nldim, REAL *beta_x,
                         REAL *beta_y, REAL *beta_z);

void mulliken (STATE *states);
REAL ylm(int l, REAL *r);
int listlen (FILE * fh, char *id);
int del_space (FILE * fh, int *tchr, int isdata);
void norm_psi1_parallel (STATE * sp, int istate, int kidx);
void print_matrix(double *b, int n, int ldb);
void sl_init(int *ictxt, int size);
void sl_exit(int ictxt);
void set_desca(int *desca, int *ictxt, int *size);
void distribute_mat(int *desca, double *bigmat, double *dismat, int *size);
void matinit(int *desca, double *dismat, double *globmat, int size);
void print_distribute_mat(double *dismat, int *desca, int size);
void init_efield (REAL * vnuc);
void pulay(int step, int N, double *xm, double *fm, int NsavedSteps, int preconditioning);


/* Some stuff for timing and performance measurements */
#define TOTAL_TIME (0)
#define ORTHO_TIME (1)
//#define NL_TIME (2)
//#define NS_TIME (3)
#define EIG_TIME (4)
#define IMAGE_TIME (5)
//#define APPCIL_TIME (6)
//#define APPCIR_TIME (7)
#define RESTRICT_TIME (8)
#define EXPAND_TIME (9)
#define PACK_TIME (10)
#define INIT_TIME (11)
#define HARTREE_TIME (12)
#define DIAG_TIME (13)
#define LFORCE_TIME (14)
#define NLFORCE_TIME (15)
#define APPGRAD_TIME (16)
#define GATHER_TIME (17)
#define MG_EIGTIME (18)
#define INTERPOLATION_TIME (19)
#define RHO_TIME (20)
#define FORCE_TIME (21)
#define SCF_TIME (22)
//#define MD_LOOP (23)
#define NLCCFORCE_TIME (24)
#define IIFORCE_TIME (25)
#define MG_EIG_NLS_TIME (26)
#define MG_EIG_APPCIL_TIME (27)
#define MG_EIG_APPCIR_TIME (28)
#define MG_EIG_TRADE_TIME (29)
#define DIAG_NL_TIME (30)
#define DIAG_APPCIL_TIME (31)
#define DIAG_APPCIR_TIME (32)
#define GET_TE_XC_TIME (33)
#define GET_TE_II_TIME (34)
#define GET_TE_TIME (35)
#define SCF_XC_TIME (36)
#define INTERP_SETUP_TIME (37)
#define INTERP_EVAL_TIME (38)
#define DIAG_SUBDIAG1_TIME (39)
#define DIAG_MATRIX_TIME (40)
#define DIAG_WAVEUP_TIME (41)
#define DIAG_SUBDIAG1_LOOP_TIME (42)
#define DIAG_APP_A (43)
#define DIAG_APP_S (44)
#define DIAG_APP_B (45)
#define DIAG_DGEMM (46)
#define DIAG_GENVPSI_TIME (47)
#define DIAG_GLOB_SUMS (48)
#define DIAG_BETAXPSI (49)
#define ALLOC_TIME (50)
#define ORTHO_BETAXPSI (51)
#define ORTHO_NORM_PSI (52)
#define ORTHO_NEW_PSI (53)
#define ORTHO_GET_COEFF (54)
#define ORTHO_GLOB_SUM (55)
#define ORTHO_UPDATE_WAVES (56)
#define DIAG_APPCIR_TIME2 (57)
#define MG_EIG_GENVPSI_TIME (58)
#define MG_EIG_EIGVALUE_TIME (59)
#define MG_EIG_APPSMOOTH_TIME (60)
#define MG_EIG_MGRIDSOLV_TIME (61)
#define MG_EIG_PACK_TIME (62)
#define DIAG_NLS_TIME (63)
#define PREINIT_TIME (64)
#define FINISH_TIME (65)
#define DIAG_SCALAPACK_INIT (66)
#define DIAG_DISTMAT (67)
#define REAL_SUM_ALL_TIME (68)
#define GLOBAL_SUMS_TIME (69)
#define DIAG_BCAST_EIGS (70)
#define READ_PSEUDO_TIME 71
#define READ_CONTROL_TIME 72

#define LAST_TIME (80)


/* Occupation flags */
#define OCC_NONE 0
#define OCC_FD 1
#define OCC_GS 2
#define OCC_EF 3


/* SMP directives for the threads */
#define     SMP_EIG       1
#define     SMP_ORTHO1    2
#define     SMP_ORTHO2    3
#define     SMP_GET_RHO   4
#define     SMP_SORT_PSI  5
#define     SMP_SKIP      6
#define     SMP_DIAG1     7
#define     SMP_DIAG2     8
#define     SMP_NLFORCE   9
#define     SMP_GETNLOP  10


/* Crystal lattice types */
/** Simple cubic lattice type.
 *  @doc Set input file value = 1 */
#define CUBIC_PRIMITIVE 	1

/** Face centered cubic lattice type. 
 *  @doc Set input file value = 2 */
#define CUBIC_FC		2

/** Bodycentered cubic lattice type. 
 *  @doc Set input file value = 3 */
#define CUBIC_BC		3

/** Hexagonal lattice type. 
 *  @doc Set input file value = 4 */
#define HEXAGONAL		4

#define TRIGONAL_PRIMITIVE	5
#define TETRAGONAL_PRIMITIVE	6
#define TETRAGONAL_BC           7

/** Orthorhombic lattice type. 
 *  @doc Set input file value = 8 */
#define ORTHORHOMBIC_PRIMITIVE  8

#define ORTHORHOMBIC_BASE_CENTRED 9
#define ORTHORHOMBIC_BC         10
#define ORTHORHOMBIC_FC 11
#define MONOCLINIC_PRIMITIVE 12
#define MONOCLINIC_BASE_CENTRED 13
#define TRICLINIC_PRIMITIVE 14

/* The real or imaginary part of a wavefunction */
#define PSI_REAL     0
#define PSI_IMAG     1






/* this declares macros and external functions for memory allocation */
#include "salloc.h"


/* other general macros */
#include "macros.h"

/* routines for input parsing */
#include "input.h"

/******/
