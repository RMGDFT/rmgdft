/**@name SPECIES
 * @memo Species (pseudopotential) control structure
 * @doc Structure holds data about the pseudopotentials used to
 * represent different types of atomic species. 
*/
typedef struct
{

    /* Text header with extra information. For UPF pseudopotentials it is the PP_INFO node */
    char *INFO;

    /* pseudopotential filename */
    char pseudo_filename[MAX_PATH];

    /** Description of the species (e.g Atomic carbon generated using 
     * hamann's code with  rcs=0.80 rcp=0.85 bohr
     */
    char *description;

    /** Atomic number */
    int atomic_number;

    /** Atomic symbol */
    char *atomic_symbol;

    /** Atomic mass */
    double atomic_mass;

    /** Norm conserving pp flag */
    int is_norm_conserving;

    /** Number of valence electrons */
    double zvalence;

    // Exchange correlation string as read from UPF pseudopotential format
    char functional[24];

    /** Gaussian charge parameter used for compensating the 1/r Coulomb
     * tail of the pseudopotentials
     */

    double rc;

    /* Number of grid points in the local in each coordinate direction. 
     * These used to be L0_LDIM and L0_NLDIM.
     */
    //int ldim;
    int ldim_coar;
    int ldim;
    int ldim_fine;
    int nldim;
    int nlfdim;
    int qdim;




    /* These are input parameters in the pseudopotential file. They represent the
     * real radii that are used in generating ldim and nldim.
     */
    double lradius;
    double nlradius;
    double qradius;

    /*Radius for milliken analysis*/
    double mill_radius;

    /*Radius in number of grid points*/
    int mill_dim;

    /*Number of radial atomic wave functions - these depend on l only, not on m*/
    int num_atomic_waves;

    /*l-numbers for states for which we have atomic orbitals*/
    int *atomic_wave_l;

    /* total number of atomic wave functions including m-dependence */
    int num_atomic_waves_m;

    double *atomic_wave_oc;
    
    //char *atomic_wave_label[3];

    double *atomic_rho;
    
    /* Pseudo atomic valence density read from PP file in log grid*/
    double **atomic_wave;
    /* Pseudo atomic valence density on linear grid*/
    double **awave_lig;


    /*Sum of all atomic states (with different l or m numbers*/
    //int sum_atomic_waves;

    /*This will store name of atomic wavefunctions, such as s, px, py, pz, dxx etc*/
    //char atomic_wave_symbol[20][12];


    /** Number of radial grid points in the pseudopotential file */
    int rg_points;

    /* Log mesh parameter, where aa=exp(-aasf)/Z, bb=1.0/bbsf */
    double aa, bb;

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
    double ddd0[MAX_NL][MAX_NL];
    double ddd[MAX_NL][MAX_NL];

    /*matrix qqq(nbeta,nbeta) */
    double qqq[MAX_NL][MAX_NL];

    /*the number of L=|l1-l2|.....|l1+l2|, we limit nlc <=5 */
    int nlc;

    /*the number of component in polynomial of the pseudized Q_I(r) function we limit nqf<=10 */
    int nqf;

    /*L-independent inner coutoff radii rinner for Q_I(r) function */
    double *rinner;

    /* ultrosoft Vanderbilt Qnm_rad(r) function and */
    double *qnm;
    double *qnmlig;

    /* the coefficient for pseudosation of Qnm_L(r) */
    double *qfcoef;

    /* Logarithmic radial mesh information */
    double *r;
    double *rab;


    /* Local Pseudopotentials */
    double *vloc0;

    /* Core charge radial grids */
    double *cr;



    /* Pseudo atomic core density */
    double *rspsco;

    /*the L-value for the beta function */
    int llbeta[MAX_NB];

    /*utrosoft Vanderbilt beta_n(r) function on radial grid */
    double *beta[MAX_NB];


    /* Total number of projectors */
    int nbeta;


    /* Log interpolation storage for the compensated local potential
     * and for it's radial derivative.
     */
    double localig[MAX_LOGGRID];

    /* Log interpolation storage for the core charge density */
    double rhocorelig[MAX_LOGGRID];

    /* Utrosoft Vandbelit Projectors on log interpolation grid */
    double betalig[MAX_NB][MAX_LOGGRID];

    // projectors in G-space radial log grid
    double *beta_g[MAX_NB];

    // local pseudopotential (rhoc part not included ) in G-space radial log grid.
    double *localpp_g;
    double *arho_g;
    double *rhocore_g;

    /* Array of r and xyz vals used to setup forward fft's for init_weight and init_derweight
     * the r_index array is dimensioned (nlfdim*nlfdim*nlfdim) while the others are nlfdim
     */
    double *r_index;
    double *x_index;
    double *y_index;
    double *z_index;

    /*Grid spacing for atomic charge density on linear grid*/
    double drlig_arho;
    
    /*Grid spacing for atomic wave functions on linear grid*/
    double drlig_awave;


    /* Pseudopotential filtering parameters */
    double qcut;                  /* Real space local cutoff for qfunctions */
    double rwidth;                /* Real-space width parameter */
    double gwidth;                /* G-space width parameter */

    /*Filtering parameters for atomic wavefunctions and charge density*/
    double acut; 
    double aradius; 
    double agwidth;
    double arwidth;

    /* radius of atomic wavefunctions and charge in terms of number of grid points*/
    int adim_rho;
    int adim_wave;


    /*Total number (of what exactly ???) */
    int num_projectors;

    fftw_complex *phase;

    /*This will store results of forward fourier transform on the coarse grid */
    fftw_complex *forward_beta;
    fftw_complex *forward_vnuc;
    fftw_complex *forward_rhoc;
    fftw_complex *forward_rhocore;

    /*Backwards wisdom for fftw */
    char *backward_wisdom;

    /*Some parameters for Q function*/
    int indv[MAX_NL];
    int nhtol[MAX_NL];
    int nhtom[MAX_NL];
    int nh_l2m[MAX_NL];
    int nh;

    /*Atomic charge density on linear grid*/
    double arho_lig[MAX_LOGGRID];
    

    int localidx;

} SPECIES;

