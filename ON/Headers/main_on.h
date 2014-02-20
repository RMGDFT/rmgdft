/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/md.h *****
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


#define ORDER_N 1


/* Version information */
#include "main.h"


/* Compile time parameters */
#include    "params_on.h"


/* Constants and symbolic definitions */
#include    "const_on.h"

/* Fourier transformation structure definition */
#include    "fftw.h"

#include    "my_finegrid.h"



/** Size of floating point variables used in QMD */
#define     REAL    double
#include    "common_prototypes.h"


int MXLLDA, MXLCOL;
REAL *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot,
    *vtot_c;
REAL *vh_old, *vxc_old;
REAL *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *Hij_00, *Bij_00;
REAL *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
REAL *work_matrix_row, *coefficient_matrix_row, *nlarray1;
REAL *projectors, *projectors_x, *projectors_y, *projectors_z;
REAL *sg_twovpsi, *sg_res;
int *nlindex;
REAL *work_memory;
REAL *sg_orbit;
REAL *sg_orbit_res;
REAL *orbit_tem;
REAL *vtot_global;


int NPES;
int NX_GRID, NY_GRID, NZ_GRID;
int P0_BASIS;
int S0_BASIS;
int PX0_GRID;
int PY0_GRID;
int PZ0_GRID;
int state_to_ion[MAX_STATES];
int state_to_proc[MAX_STATES];
char num_loc_st[MAX_IONS];
short int overlap_or_not[MAX_IONS * MAX_IONS];
short int state_overlap_or_not[MAX_STATES * MAX_STATES];

struct b_list
{
    char *TagName;
    int Flag;
};

/* option tags list */
struct t_list
{
    char *TagName;
    char *OptName;
    int OptCount;
    char *Opt[MAX_OPTS];
};


/* Array telling for each pair of ions if the localization regions
   around touch both a same region (anywhere on any PE) */
int array_both_not_zero[MAX_IONS][MAX_IONS];


STATE states[MAX_STATES];
STATE states1[MAX_STATES];
STATE states_res[MAX_STATES];
STATE states_res1[MAX_STATES];
STATE states_tem[MAX_STATES];


/* Header file for blas routines */
#include "blas3_def.h"
#include "lapack_def.h"



/* Function prototypes */
void app_4del2 (REAL * f, REAL * work);

void app_nl (REAL * psiR, REAL * psiI, REAL * workR, REAL * workI, int state,
        int flag, int kidx, int tid);
void cholesky (REAL * a, int n);
void cross_product (REAL * a, REAL * b, REAL * c);
REAL fill (STATE * states, REAL width, REAL nel, REAL mix,
        int num_st, int occ_flag);
void genvpsi (REAL * psi, REAL * twovpsi, REAL * pvtot, REAL * pvnl,
        REAL * kd, REAL kmag, int dimx, int dimy, int dimz);
void get_index_loc (STATE *);
void get_eig (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc);
char *get_num (char *str);
void get_zdens (STATE * states, int state, REAL * zvec);
void xclda_pz81 (REAL * rho, REAL * vxc);
void exclda_pz81 (REAL * rho, REAL * exc);
void xcgga (REAL * rho, REAL * vxc, REAL * exc, int flag);
void gram (KPOINT * kpoint, REAL h, int numst, int maxst, int numpt,
        int maxpt);
REAL get_ke (STATE * sp, int tid);
void global_sums (REAL * vect, int *length, MPI_Comm comm);
void init_kbr (void);
void init_pe_on (void);
void init_pegrid (void);
void init_wf (STATE * states);
void init_wflcao (STATE * states);
void init_nuc (REAL * vnuc, REAL * rhoc, REAL * rhocore);
void init_pos ();
void init_sym (void);
void symmetrize_rho (REAL * rho);
void symforce (void);
void rmg_timings (int what, REAL time);
void mg_eig_state (STATE * states, int tid, REAL * vtot);
REAL minimage (ION * ip1, ION * ip2, REAL * xtal_r);
REAL my_crtc (void);
/* Local function prototypes */
void norm_psi (STATE * psi, REAL bv);
void ortho_full (STATE * states);
void ortho_half (STATE * states);
void ortho_bcc (STATE * states);
void output (STATE * states, int its);
void pe2xyz (int pe, int *x, int *y, int *z);
void pack_ptos (REAL * sg, REAL * pg, int dimx, int dimy, int dimz);
void pack_stop (REAL * sg, REAL * pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (REAL * sg, REAL * pg, REAL alpha, int dimx, int dimy,
        int dimz);
void pack_ptos_trade (REAL * sg, REAL * pg, int dimx, int dimy, int dimz);
void pack_vhstod (REAL * s, REAL * d, int dimx, int dimy, int dimz);
void pack_vhdtos (REAL * s, REAL * d, int dimx, int dimy, int dimz);
double radint (double *f, double *r, int n, double al);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
REAL real_sum_all (REAL x, MPI_Comm comm);
REAL real_max_all (REAL x);
void sortpsi (STATE * states);
void subdiag (STATE * states, REAL * vh, REAL * vnuc, REAL * vxc);
void trade_images (REAL * mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_mpi (REAL * mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_images_smp (REAL * mat, int dimx, int dimy, int dimz, int *nb_ids);
void set_bc (REAL * mat, int dimx, int dimy, int dimz, int images, REAL val);
void getpoi_bc (REAL * rho, REAL * vh_bc, int dimx, int dimy, int dimz);
void vol_rho (REAL * rho, int step);
void vol_wf (STATE * states, int state, int step);
void write_avgd (REAL * rho);
void write_avgv (REAL * vh, REAL * vnuc);
void write_zstates (STATE * states);
void write_header (void);
void write_pos (void);
void write_eigs (STATE * states);
void write_occ (STATE * states);
void write_timings (void);
REAL rand0 (long *idum);

void gather_psi (REAL * tmp_psiR, REAL * tmp_psiI, STATE * sp, int tid);
void scatter_psi (REAL * tmp_psiR, REAL * tmp_psiI, STATE * sp, int tid);
void get_milliken (STATE * states);

void output_wave (STATE * states, int kpt, int fhand);



/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, REAL * celldm, REAL * A0I, REAL * A1I, REAL * A2I,
        REAL * OMEGAI, int *flag);
void recips (void);
void to_cartesian (REAL crystal[], REAL cartesian[]);
void to_crystal (REAL crystal[], REAL cartesian[]);
REAL metric (REAL * crystal);

/* Md run types */
void run (STATE * states, STATE * states1);
void fastrlx( STATE *states, STATE *states1, REAL *vxc, REAL *vh, REAL *vnuc, REAL *vh_old, REAL *vxc_old,
        REAL *rho, REAL *rhocore, REAL *rhoc );
void cdfastrlx (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc,
        REAL * rho, REAL * rhocore, REAL * rhoc);
/*void moldyn (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc, REAL * rho, REAL * rhoc, REAL * rhocore);*/
void dx (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc, REAL * rho,
        REAL * rhoc);
void psidx (STATE * states, REAL * vxc, REAL * vh, REAL * vnuc, REAL * rho,
        REAL * rhoc);
void cholesky (REAL * a, int n);


REAL minimage1 (REAL aa[3], REAL bb[3]);
void init_parameter (STATE * states);
void get_mehr (void);
void make_mask_grid (REAL rcut, int level, STATE * states);
void make_mask_grid_state (int level, STATE * states);
void allocate_func (STATE * states, int inist);
void allocate_psi (STATE * states, STATE * states1);
void allocate_matrix ();
void xyz2pe (int x, int y, int z, int *pe);
void get_normKB (SPECIES * sp, double *pd);
void get_start (int cdim, double crds, double cstart, double hgrid,
        int *istart, double *nlcstart);
void get_index_array (int *pvec, int *dvec, int dim,
        int *Aix, int *Aiy, int *Aiz, int *icount,
        int index_low[3], int index_high[3]);
void get_Ai (int cdim, int *Ai, int ngrid, int istart);
void xclsd_pz81 (REAL * rho, REAL * vxc);
void global_sums_int (int *, int *);
REAL linint (REAL * y, REAL r, REAL invdr);
void ortho_norm_local (STATE * states);
void global_sums_int (int *vect, int *length);
void my_barrier (void);
void matrix_and_diag (STATE * states, STATE * states1, REAL * vxc, int flag);
void get_kbpsi (STATE *sp1, double *kbpsi);
void precond_mg (double *res, double *work1, double *work2, int istate);

void matS_cholesky_real (STATE * states);
void get_invBmat (STATE * states, double *matB);
void print_matrix (double *b, int n, int ldb);
void get_overlap_real (double *aa, int numst, int numpt,
        int lda, double *ss, int lds);
void get_cholesky_real (double *matS);
void get_all_kbpsi (STATE * states1, STATE * states);
void get_Hvnlij (double *Aij);
void genvlocpsi (REAL * psi, int st1, REAL * work1, REAL * vtot_global, STATE * states);
void genvnlpsi (double *sg_twovpsi, double *vnl,
        int dimx, int dimy, int dimz);
void get_new_rho (STATE * states, double *rho);
void get_wave(int st, STATE * states);
void add_orbit_to_wave(int st1, REAL scale, REAL * psi1, REAL * wave_global, STATE * states);
void print_wave(int wave_plot, STATE * states, int coarse_level);
void mg_eig (STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, REAL * vxc_old,
        REAL * vh_old);
void fsymforces (REAL * force, int *s, int *irg, int *irt, int *nat,
        int *ibrav, int *nsym, REAL * celldm, int *nr1, int *nr2,
        int *nr3);
void symmetry (int *ibrav, int *s, int *nsym, int *irg, int *irt, int *ftau,
        int *nat, REAL * tau, int *ityp, int *nks, REAL * xk,
        REAL * wk, REAL * celldm, int *nr1, int *nr2, int *nr3,
        int *wflag);
void symrho (REAL * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
        int *irg, int *ftau);

void sgama (int *nrot,int *nat,double *s,double *at,double *bg,double *tau,
            int *ityp,int *nsym,int *nr1,int *nr2,int *nr3,int *irg,
            int *irt,double *ftau,double *rtau,int *npk,int *nks,double *xk,
            double *wk,double *xau,double *rau,bool *invsym, int *wflag);

void line_min_three_point (STATE *, STATE *, REAL, REAL, REAL *, REAL *,
        REAL *, REAL *, REAL *);

void dot_product_orbit_orbit (STATE *orbit1, STATE *orbit2, STATE
*orbit3, double *H, double *S);

void orbit_dot_orbit (STATE * states, STATE * states1, REAL *Hij_row, REAL * Bij_row);

void app_mask (int istate, double *u, int level);
void allocate_masks (STATE * states);

void state_corner_xyz (STATE * states);

void density_orbit_X_orbit (int st1, int st2, REAL scale, REAL * psi1,
        REAL * psi2, REAL * rho_global, int mode,
        STATE * states);

char *get_symbol (int atomic_number);

void get_matB_qnm (double *Aij);
void pack_vtot_ftoc (REAL * vtot, REAL * vtot_c);
void get_qnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void qnm_beta_betapsi (STATE *sp, int ion2, REAL * prjptr, REAL * work);
void get_dnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void dnm_beta_betapsi (STATE *sp, int ion2, REAL scale, int ip2, REAL * prjptr,
        REAL * work);
void pack_rho_ctof (REAL * rho1, REAL * rho_f);
void rho_augmented (REAL * rho, REAL * global_mat_X);
void rho_Qnm_mat (double *Aij, REAL * global_mat_X);
void rho_nm_mat (double *Aij, REAL * global_mat_X);
int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz, int *ilow, int *ihi,
        int *jlow, int *jhi, int *klow, int *khi, int cdim, int pxgrid,
        int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid,
        REAL * xcstart, REAL * ycstart, REAL * zcstart);
void xcgga (REAL * rho, REAL * vxc, REAL * exc, int mode);
REAL radint1 (REAL * func, REAL * r, REAL * rab, int n);
void get_all_partial_kbpsi (STATE *states);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, REAL * global_mat_X);
void partial_QI (int ion, REAL * QI_R, ION *iptr);
void partial_nlop_s (ION * iptr, REAL * betax, REAL * betay, REAL * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void partial_nlop_p (ION * iptr, REAL * betax, REAL * betay, REAL * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void partial_nlop_d (ION * iptr, REAL * betax, REAL * betay, REAL * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void get_mat_Omega (STATE * states, double Omega[]);
void md_fastrelax(void);
void change_states_crds (STATE * states);

/*end shuchun wang */

FILE *open_xbs_movie (char *filename);


/* Occupation flags */
#define OCC_NONE 0
#define OCC_FD 1
#define OCC_GS 2
#define OCC_EF 3



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

/* The number of possible point symmetries */
#define MAX_SYMMETRY	48

/******/






REAL get_te_ion_ion ();
REAL get_sum_eig (STATE * states);
REAL get_Exc (REAL * rho, REAL * rhocore);
void correct_res (STATE *, REAL *);
void get_state_to_proc (STATE * states);


#include "overlap.h"



/* different methods to update orbitals */
void kain (int step, int N, double *xm, double *fm, int NsavedSteps);
void pulay (int step, int N, double *xm, double *fm, int NsavedSteps,
        int preconditioning);
void sd (int step, int N, double *xm, double *fm);
void pulay_mixing (int size, double *rho, int NsavedSteps);
void charge_pulay (int step, int N, double *xm, double *fm, int NsavedSteps);
void pulay_kain (int step, int N, double *xm, double *fm, int NsavedSteps);

/*  Moving center stuff  */
void get_orbit_center (STATE *state, REAL * x, REAL * y, REAL * z);
void update_orbit_centers (STATE * states);
int if_update_centers (STATE * states);

void write_states_info (char *outfile, STATE * states);
void read_states_info (char *outfile, STATE * states);

double get_gamma_precond (double *vtot, double small_eig);
void init_state_size (STATE * states);
void diag_eig_maitrix(double *, double *, int *);



void filter_potential (REAL *potential, REAL *r, int rg_points, REAL rmax, REAL offset, REAL parm, REAL* potential_lgrid, 
	REAL *rab, int l_value, REAL dr, REAL  gwidth, int lgrid_points, REAL rcut, REAL rwidth, REAL * drpotential_lgrid);


#include "macros.h"
#include    "prototypes_on.h"
