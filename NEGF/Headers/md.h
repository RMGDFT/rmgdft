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

#include <complex.h>

#define NEGF1 1

/* Compile time parameters */
#include    "params.h"

/* GPU includes */
#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif

#include "my_mpi.h"
/* Version information */
#include    "version.h"
#include    "input.h"







/* Constants and symbolic definitions */
#include    "const.h"

/* Fourier transformation structure definition */
#include    "fftw.h"

#include    "my_finegrid.h"  /*shuchun*/

#include "method.h"

#include "my_scalapack.h"




/* Some definitions needed for using compile time values for global grid sizes */
#include "fixed_dims.h"


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

#define MAX_TRADE_IMAGES 5


/* The number of possible point symmetries */
#define MAX_SYMMETRY	48

/******/



int MXLLDA, MXLCOL;
REAL *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *vext,  *vcomp, *rhocore, *vtot, *vtot_old;
REAL *vh_old, *vxc_old; 
REAL *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
REAL *work_dis2, *zz_dis, *gamma_dis, *uu_dis; 
REAL *work_matrix, *nlarray1;
REAL *projectors, *projectors_x, *projectors_y, *projectors_z;
REAL *sg_twovpsi, *sg_res;
int *nlindex;
REAL *work_memory;
REAL *sg_orbit;
REAL *sg_orbit_res;
REAL *orbit_tem;
REAL *vtot_global;
complex double *sigma_all;


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
short int state_overlap_or_not[MAX_STATES * MAX_STATES];

#include "blas.h"
#include "blas3_def.h"
#include "lapack_def.h"



/* Function prototypes */
void app_4del2(REAL *f, REAL *work);
REAL app_del2c(REAL *a, REAL *b, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);

void app_smooth(REAL *f, REAL *work, REAL sfac);
void app_cir(REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_0(REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_ortho(REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_ortho_0(REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_bcc(REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_fcc(REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_hex(REAL *a, REAL *b, int dimx, int dimy, int dimz);
REAL app_cilr(REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_bcc(REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_fcc(REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_hex(REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cilr_ortho(REAL *a, REAL *b, REAL *c, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil(REAL *a, REAL *b, int dimx, int dimy, int dimz,
        REAL gridhx, REAL gridhy, REAL gridhz);
void app_nl(REAL *psiR, REAL *psiI, REAL *workR, REAL *workI,
        int state, int flag, int kidx, int tid);
void cholesky(REAL *a, int n);
void cross_product(REAL *a, REAL *b, REAL *c);
void eval_residual( REAL *mat, REAL *f_mat, int dimx, int dimy, int dimz, 
        REAL gridhx, REAL gridhy, REAL gridhz, REAL *res );
void solv_pois( REAL *vmat, REAL *fmat, REAL *work, 
        int dimx, int dimy, int dimz, REAL gridhx, 
        REAL gridhy,REAL gridhz, double step, double k);
REAL fill (STATE *states, REAL width, REAL nel, REAL mix, int num_st, int occ_flag);
void genvpsi(REAL *psi, REAL *twovpsi, REAL *pvtot, REAL *pvnl, REAL *kd, 
        REAL kmag, int dimx, int dimy, int dimz);
void get_index_loc(STATE *);
void get_nlop(void);
void get_phase(ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_smp(int tid);
void get_nlop_s(ION *iptr, REAL *rtptr, int ip);
void get_nlop_p(ION *iptr, REAL *rtptr, int ip);
void get_nlop_d(ION *iptr, REAL *rtptr, int ip);
void get_nlop_f(ION *iptr, REAL *rtptr, int ip);
void get_eig(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc);
char *get_num(char *str);
void get_rho(STATE *states, REAL *rho);
void get_te(double *rho, double *rhocore, STATE *states); 
void get_vxc(REAL *rho, REAL *rhocore, REAL *vxc);
void get_zdens(STATE *states, int state, REAL *zvec);
void xclda_pz81(REAL *rho, REAL *vxc);
void exclda_pz81(REAL *rho, REAL *exc);
void xcgga(REAL *rho, REAL *vxc, REAL *exc, int flag);
void gram(KPOINT *kpoint, REAL h, int numst, int maxst, int numpt, int maxpt);
REAL get_ke(STATE *sp, int tid);
void global_sums (REAL *vect, int *length, MPI_Comm comm);
void iiforce(void);
void init(REAL *vh, REAL *rho, REAL *rhocore, REAL *rhoc, STATE *states,
        STATE *states1, REAL *vnuc, REAL *vxc, REAL *vh_old, REAL *vxc_old);
void init_kbr(void);
void init_pe_on(void);
void init_pegrid(void);
void init_wf(STATE *states);
void init_wflcao(STATE *states);
void init_nuc(REAL *vnuc, REAL *rhoc, REAL *rhocore);
void init_ext (double *vext, double gbias_begin, double gbias_end, double BT, double gate_bias);
void init_pos();
void init_nonlocal_comm(void);
void init_sym(void);
void symmetrize_rho(REAL *rho);
void symforce(void);
void lforce(REAL *rho, REAL *vh);
void nlccforce(REAL *rho, REAL *vxc);
void cforce(REAL *rho, REAL *vh);
void rmg_timings(int what, REAL time);
void mg_eig_state(STATE *states, int tid, REAL *vtot);
REAL minimage(ION *ip1, ION *ip2, REAL *xtal_r);
REAL my_crtc(void);
/* Local function prototypes */
void norm_psi(STATE *psi, REAL bv);
void ortho_full(STATE *states);
void ortho_half(STATE *states);
void ortho_bcc(STATE *states);
void output(STATE *states, int its);
void pe2xyz(int pe, int *x, int *y, int *z);
void pack_ptos(REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void pack_stop(REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void pack_stop_axpy(REAL *sg, REAL *pg, REAL alpha, int dimx, int dimy, int dimz);
void pack_ptos_trade(REAL *sg, REAL *pg, int dimx, int dimy, int dimz);   
void pack_vhstod(REAL *s, REAL *d, int dimx, int dimy, int dimz);
void pack_vhdtos(REAL *s, REAL *d, int dimx, int dimy, int dimz);
double radint(double *f, double *r, int n, double al);
void radiff(double *f, double *df, double *r, int n, double al);
void ra2diff(double *f, double *df, double *r, int n, double al); 
void ranv(void);
void read_control(void);
void read_data(char *name, REAL *vh, REAL *rho, REAL *vxc, REAL *vh_old, REAL *vxc_old, STATE *states);
void read_rho_and_pot(char *name, double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void read_pseudo(void);
void read_LCR(void);
void read_cond_input(double *emin, double *emax, int *E_POINTS, double *E_imag, double *KT, int *kpoint );
REAL real_sum_all(REAL x, MPI_Comm comm);
REAL real_max_all(REAL x);
void rft(REAL *f, REAL *r, REAL *ffil, REAL al, int rg_points, int lval,
        REAL dr, REAL width, int lrg_points);
void scf(complex double *sigma_all, STATE *states, STATE *states1, REAL *vxc, REAL *vh, REAL *vnuc, REAL *vext,
        REAL *rho, REAL *rhocore, REAL *rhoc,REAL *vxc_old, REAL *vh_old, REAL *vbias, int *CONVERGENCE);
void apply_potential_drop(REAL *vbias);
void sortpsi(STATE *states);
void subdiag(STATE *states, REAL *vh, REAL *vnuc, REAL *vxc);
void trade_images( REAL *mat, int dimx, int dimy, int dimz, int *nb_ids, int type );
void trade_images_mpi( REAL *mat, int dimx, int dimy, int dimz, int *nb_ids );
void trade_images_smp( REAL *mat, int dimx, int dimy, int dimz, int *nb_ids );
void set_bc( REAL *mat, int dimx, int dimy, int dimz, int images, REAL val );
void getpoi_bc(REAL *rho, REAL *vh_bc, int dimx, int dimy, int dimz);
void trade_images2(REAL *f, REAL *w, int dimx, int dimy, int dimz);
void trade_images3(REAL *f, REAL *w);
void vol_rho(REAL *rho, int step);
void vol_wf(STATE *states, int state, int step);
void write_avgd(REAL *rho);
void write_avgv(REAL *vh, REAL *vnuc);
void write_zstates(STATE *states);
void write_data(char *name, REAL *vh, REAL *rho, REAL *vxc, REAL *vh_old, REAL *vxc_old, REAL *vbias, STATE *states);
void write_header(void);
void write_pos(void);
void write_eigs(STATE *states);
void write_occ(STATE *states);
void write_force(void);
void write_timings(void);

REAL rand0(long *idum);

void mg_restrict (REAL *full, REAL *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_restrict_6 (REAL * full, REAL * half, int dimx, int dimy, int dimz, int grid_ratio);
void mg_prolong (REAL *full, REAL *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

void gather_psi(REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid);
void scatter_psi(REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid);
void get_milliken(STATE *states);

void bandstructure(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
        REAL *rho, REAL *rhocore, REAL *rhoc);
void output_wave(STATE *states, int kpt, int fhand);

/* Blas wrappers */
void QMD_daxpy(int n, REAL alpha, REAL *x, int incx, REAL *y, int incy);
void QMD_dscal(int n, REAL alpha, REAL *x, int incx);
void QMD_dcopy(int n, REAL *x, int incx, REAL *y, int incy);
REAL QMD_ddot(int n, REAL *x, int incx, REAL *y, int incy);



int get_index(int gridpe, ION *iptr, int *Aix, int *Aiy, int *Aiz,
        int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
        int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
        int nxgrid, int nygrid, int nzgrid,
        REAL *lxcstart, REAL *lycstart, REAL *lzcstart);


/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, REAL *celldm, REAL *A0I, REAL *A1I, REAL *A2I, REAL *OMEGAI, int *flag);
void recips(void);
void to_cartesian(REAL crystal[], REAL cartesian[]);
void to_crystal(REAL crystal[], REAL cartesian[]);
REAL metric(REAL *crystal);

/* Md run types */
void run(STATE *states, STATE *states1, STATE *states_distribute); 
void quench(STATE *states, STATE *states1, STATE *states_dis, REAL *vxc, REAL *vh, REAL *vnuc, REAL *vext, REAL *vh_old, REAL *vxc_old, REAL *rho, 
        REAL *rhocore, REAL *rhoc, REAL *vbias);
void fastrlx(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
        REAL *rho, REAL *rhocore, REAL *rhoc);     
void cdfastrlx(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
        REAL *rho, REAL *rhocore, REAL *rhoc);
void moldyn(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
        REAL *rho, REAL *rhoc, REAL *rhocore);
void dx(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
        REAL *rho, REAL *rhoc);
void psidx(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
        REAL *rho, REAL *rhoc);
void cholesky(REAL *a, int n);


REAL minimage1(REAL aa[3], REAL bb[3]);
void init_parameter(STATE *states);
void get_mehr(void);
void make_mask_grid(REAL rcut, int level, STATE *states); 
void make_mask_grid_state(int level, STATE *states); 
void allocate_func(STATE *states, int inist);
void allocate_psi(STATE *states, STATE *states1);
void allocate_matrix();
void xyz2pe(int x, int y, int z, int *pe);
void get_normKB(SPECIES *sp, double *pd);
void get_start(int cdim, double crds, double cstart, double hgrid, 
        int *istart, double *nlcstart);
void get_index_array(int *pvec, int *dvec, int dim, 
        int *Aix, int *Aiy, int *Aiz, int *icount,
        int index_low[3], int index_high[3]);
void get_Ai(int cdim, int *Ai, int ngrid, int istart);
void xclsd_pz81(REAL *rho, REAL *vxc);
void global_sums_int (int*, int*);
REAL linint(REAL *y, REAL r, REAL invdr);
void ortho_norm_local(STATE *states);
void app_mask(int istate, double *u, int level);
void global_sums_int (int *vect, int *length);
void my_barrier(void);
void matrix_and_diag(STATE *states, STATE *states1, REAL *vxc);
void get_kbpsi(STATE *sp1, double *kbpsi);
void precond_mg(double *res, double *work1, double *work2, int istate);

void matS_cholesky_real(STATE *states);
void get_invBmat(STATE *states, double* matB);
void get_Hij(STATE *states1, STATE *states, double *vtot, double *Aij);
void get_hb(double matS[], double theta[], double hb[],
        int dim, int ldim);
void print_matrix(double *b, int n, int ldb);
void get_cholesky_real(double *matS);
void get_matB(STATE *states, STATE *states1, double *mat);
void get_all_kbpsi(STATE *states1, STATE *states);
void get_Hvnlij( double *Aij);
void genvlocpsi(REAL *psi, int st1, REAL *work1, REAL *vtot_global, STATE *states);
void get_vnlpsi(STATE *sp1, double *kbpsi, double *work);
void genvnlpsi(double *sg_twovpsi, double *vnl,
        int dimx, int dimy, int dimz);
void get_new_rho(STATE *states, double *rho);
void get_new_rho2(STATE *states, double *rho);
void mg_eig(STATE *states, STATE *states1, double *vxc, double *vh, double *vnuc,
        double *rho, double *rhoc, REAL *vxc_old, REAL *vh_old);
void fsymforces(REAL *force, int *s, int *irg, int *irt,
        int *nat, int *ibrav, int *nsym, 
        REAL *celldm, int *nr1, int *nr2, int *nr3);
void symmetry(int *ibrav, int *s, int *nsym, int *irg, int *irt,
        int *ftau, int *nat, REAL *tau, int *ityp, int *nks,
        REAL *xk, REAL *wk, REAL *celldm, int *nr1, int *nr2, int *nr3, 
        int *wflag);
void symrho(REAL *rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
        int *irg, int *ftau);                                              
void line_min_three_point(STATE*, STATE*, REAL, REAL, REAL*, REAL*, REAL*, REAL*, REAL*); 
void charge_density_matrix();
char *get_symbol(int atomic_number);


/* ===BEGIN=== change from ION to STATE --Qingzhong */

REAL  	dot_product_orbit_orbit(STATE *orbit1, STATE *orbit2);
void 	orbit_dot_orbit(STATE *states, STATE *states1, REAL *work_matrix);
void   	allocate_masks(STATE *states);
void 	state_corner_xyz(STATE *states);
void 	density_orbit_X_orbit(int st1, int st2, REAL scale,  REAL *psi1, REAL *psi2,
        REAL *rho_global, int mode, STATE *states);

/* ===END=== change from ION to STATE --Qingzhong */




/*begin shuchun wang */
void get_nlop_soft_s(ION *iptr, REAL *rtptr, int ip, fftwnd_plan p1, fftwnd_plan p2);
void get_nlop_soft_p(ION *iptr, REAL *rtptr, int ip, fftwnd_plan p1, fftwnd_plan p2);
void get_nlop_soft_d(ION *iptr, REAL *rtptr, int ip, fftwnd_plan p1, fftwnd_plan p2);
void get_matB_qnm(double *Aij);
void pack_vtot_ftoc(REAL *vtot, REAL *vtot_c);
void get_qnmpsi(STATE *sp, double *kbpsi_one_state, double *work);
void qnm_beta_betapsi(STATE *sp, int ion2, REAL scale, int ip2, REAL *prjptr, REAL *work);
void get_dnmpsi(STATE *sp, double *kbpsi_one_state, double *work);
void dnm_beta_betapsi(STATE *state1, int ion2, REAL scale, int ip2, REAL *prjptr, REAL *work);
void pack_rho_ctof(REAL *rho1,REAL *rho_f);
void rho_augmented(REAL *rho, REAL *global_mat_X);
void rho_Qnm_mat(double *Aij, REAL *global_mat_X);
void trade_images2(REAL *f, REAL *w, int dimx, int dimy, int dimz);
void app_grad(REAL *f, REAL *wx, REAL *wy, REAL *wz, int dimx, int dimy, int dimz);
void app6_del2(REAL *f, REAL *work, int dimx, int dimy, int dimz, REAL hxgrid, REAL hygrid, REAL hzgrid);
void get_ddd(REAL *veff);
void get_qqq();
void init_qfunct(void);
void rft1(REAL cparm, REAL *f, REAL *r, REAL *ffil, REAL *rab, int rg_points,
        int lval, REAL dr, REAL width, int lrg_points);
void get_QI(void);
void aainit(int lli, int mix, int lx, int mx, int nlx, REAL ap[][9][9], int lpx[][9], int lpl[][9][9]);
void ylmr2(double *r, double *ylm);
REAL qval(int ih, int jh, REAL r, REAL invdr, REAL *ptpr, int *nhtol, int *nhtom,
        int *indv, REAL *ylm, REAL ap[][9][9], int lpx[][9], int lpl[][9][9], SPECIES *sp);
REAL get_QnmL(int idx, int ltot, REAL r, SPECIES *sp);
void assign_weight(SPECIES *sp, fftw_complex *weptr, REAL *rtptr);
void pack_gftoc(SPECIES *sp, fftw_complex *gwptr, fftw_complex *gbptr);
void xcgga(REAL *rho, REAL *vxc, REAL *exc, int mode);
REAL radint1(REAL *func, REAL *r, REAL *rab,int n);
void allocate_matrix_soft();
void get_Hvnlij_soft( double *Aij);
void get_matB_soft(STATE *states, STATE *states1, double *mat);
void get_new_rho_soft(STATE *states, double *rho);
void get_nlop_soft(void);
void get_vh_soft(REAL *rho, REAL *rhoc, REAL *vh, REAL *vh_old, int cycles, int maxlevel);
void get_vxc_soft(REAL *rho, REAL *rhocore, REAL *vxc);
void init_nuc_soft(REAL *vnuc, REAL *rhoc, REAL *rhocore);
void init_psp_soft(void);
void init_soft(REAL *vh, REAL *rho, REAL *rhocore, REAL *rhoc, STATE *states,
        STATE *states1, REAL *vnuc, REAL *vext, REAL *vxc, REAL *vh_old, REAL *vxc_old, STATE *state_distribute);
void init_comp (REAL *vh);
void read_pseudo_soft(void);
REAL get_Exc_soft(REAL *rho, REAL *rhocore);
void fill_orbit_borders4(REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void app10_del2(REAL *f, REAL *work, int dimx, int dimy, int dimz, REAL hxgrid, REAL hygrid, REAL hzgrid);
void fill_orbit_borders3(REAL *sg, REAL *pg, int dimx, int dimy, int dimz);
void app8_del2(REAL *f, REAL *work, int dimx, int dimy, int dimz, REAL hxgrid, REAL hygrid, REAL hzgrid);
void read_matrix_tri(char *mfile, double *matrix, int n1);
void mix_vh(REAL *rho_out, REAL *rho_in, REAL mix, int steps, int mode);
void mix_vxc(REAL *rho_out, REAL *rho_in, REAL mix, int steps, int mode);
void mix_rho(REAL *rho_out, REAL *rho_in, REAL mix, int steps, int mode);
void precond_rho(double *res);
/*end shuchun wang */

REAL get_te_ion_ion();
REAL get_sum_eig(STATE *states);
REAL get_Exc(REAL *rho, REAL *rhocore);
void correct_res(STATE*, REAL*);
void mix_rho(REAL *rho_out, REAL *rho_in, REAL mix, int steps, int mode);
void get_state_to_proc(STATE *states);

/* different methods to update orbitals --Qingzhong */
void    kain(int step, int N, double *xm, double *fm, int NsavedSteps);
void    pulay(int step, int N, double *xm, double *fm, int NsavedSteps, int preconditioning);
void    sd(int step, int N, double *xm, double *fm);
void    pulay_mixing(int size, double *rho, int NsavedSteps);
void    charge_pulay(int step, int N, double *xm, double *fm, int NsavedSteps);
void    pulay_kain(int step, int N, double *xm, double *fm, int NsavedSteps);

/*  Moving center stuff  --Qingzhong */
void    get_orbit_center(STATE *state, REAL *x, REAL *y, REAL *z);
void    update_orbit_centers(STATE *states);
int     if_update_centers(STATE *states);

void    write_states_info(char *outfile, STATE *states);
void    read_states_info(char *outfile, STATE *states);

double  get_gamma(double *vtot, double small_eig);
void    init_state_size(STATE *states);

void get_H_transfer(STATE *states, STATE *states1, REAL *vh, REAL *vxc, REAL *vnuc);
void get_S_transfer(STATE *states, STATE *states1);


void    spline_p(int N, double *y, double h, double h_new, double x0, double *new_y); 
void    spline_n(int N, double *y, double h, double h_new, double x0, double *new_y); 


void    print_sum_idx(int size, double *data, char *msg); 
void    print_data_to_file(int size, double *data, char *filename); 

void   	write_global_data(int file_handle, double *data, int NX, int NY, int NZ);
void   	write_global_data_lead(int file_handle, double *data, int NX, int NY, int NZ);




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

