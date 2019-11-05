#ifndef RMG_SPECIES_H
#define RMG_SPECIES_H 1

#include "Pw.h"

class SPECIES
{

public:

    /* Text header with extra information. For UPF pseudopotentials it is the PP_INFO node */
    std::string INFO;

    /* pseudopotential filename */
    std::string pseudo_filename;

    /* Pseudopotential file format. Types are UPF (0) and XML (QMC) (1) */
    int ftype;

    /* Grid type, 0=log, 1=linear */
    int gtype;

    /** Description of the species (e.g Atomic carbon generated using 
     * hamann's code with  rcs=0.80 rcp=0.85 bohr
     */
    std::string description;

    /** Atomic number */
    int atomic_number;

    /** Atomic symbol */
    char *atomic_symbol;

    /** Atomic mass */
    double atomic_mass;

    /** Norm conserving pp flag */
    int is_norm_conserving;

    /** semi-local pp flag */
    int is_semi_local;

    /** is included in ldaU projections */
    int is_ldaU;

    /** Number of ldaU orbitals */
    int num_ldaU_orbitals;

    /** ldaU l value */
    int ldaU_l;

    // Hubbard U and J
    double Hubbard_U;
    double Hubbard_J;

    /** Number of valence electrons */
    double zvalence;

    // Exchange correlation string as read from UPF or XML pseudopotential format
    std::string functional;

    // Exchange correlation type string as read from XML pseudopotential format
    std::string functional_type;

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
    
    /* Pseudo atomic wavefunctions read from PP file in log grid*/
    double **atomic_wave;

    /* Pseudo atomic wavefunctions on G-space log grid */
    double *atomic_wave_g[MAX_INITWF];

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
    bool is_ddd_diagonal;   // Per species. We also have is_ddd_non_diagonal in the main control structure

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

    /* Projectors on radial grid. KB for norm-conserving and beta_n(r) for Vanderbilt Ultrasoft. */
    double *beta[MAX_NB];

    /* Difference potentials for semi-local dVl = V_l - V_local */
    double *dVl[MAX_L];

    /* l-value associated with each difference potential */
    int dVl_l[MAX_L];

    /* Total number of radial projectors */
    int nbeta;


    /* Log interpolation storage for the compensated local potential
     * and for it's radial derivative.
     */
    double localig[MAX_LOGGRID];

    /* Log interpolation storage for the core charge density */
    double rhocorelig[MAX_LOGGRID];

    // projectors in G-space radial log grid
    double *beta_g[MAX_NB];
    double *rbeta_g[MAX_NB][MAX_L+1];

    // local pseudopotential (rhoc part not included ) in G-space radial log grid.
    double *localpp_g;
    double *der_localpp_g;
    double *arho_g;
    double *rhocore_g;

    /*Grid spacing for atomic charge density on linear grid*/
    double drlig_arho;
    
    /*Grid spacing for atomic wave functions on linear grid*/
    double drlig_awave;

    // Flag indicating whether a particular orbital is included in LDA+U projectors
    bool *awave_is_ldaU;

    /* Pseudopotential filtering parameters */
    double qcut;                  /* Real space local cutoff for qfunctions */
    double rwidth;                /* Real-space width parameter */
    double gwidth;                /* G-space width parameter */

    /*Filtering parameters for atomic wavefunctions and charge density*/
    double *aradius; 
    double agwidth;
    double arwidth;

    /* radius of atomic wavefunctions and charge in terms of number of grid points*/
    int adim_rho;
    int adim_wave;


    /* Total number of 3D projectors for this atomic species */
    int num_projectors;

    /* Total number of 3D atomic orbitals for this atomic species */
    int num_orbitals;

    fftw_complex *phase;

    /*This will store results of forward fourier transform on the coarse grid */
    fftw_complex *forward_beta=NULL;
    fftw_complex *forward_beta_r[3];
    fftw_complex *forward_orbital=NULL;

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

    /* Point to plane wave object used for localized projectors */
    Pw *prj_pwave;

    /* Point to local BaseGrid object */
    BaseGrid *OG;
    
};
#endif
