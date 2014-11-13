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

#include <stdbool.h>

#if __cplusplus
extern "C" {
#endif

/* Function prototypes */
void app6_del2 (double *rho, double *work, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz);
void app_smooth (double *f, double *work, int dimx, int dimy, int dimz);
void app_smooth_f (rmg_float_t * f, rmg_float_t * work, int dimx, int dimy, int dimz);
void app_smooth1 (double *f, double *work, int dimx, int dimy, int dimz);
void app_smooth1_f (rmg_float_t *f, rmg_float_t *work, int dimx, int dimy, int dimz);
void app_cir_driver (double *a, double *b, int dimx, int dimy, int dimz, int order);
void app_cir_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, int order);
void app_cir_fourth (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_fourth_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz);
void app_cir_sixth (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_sixth_f (rmg_float_t *a, rmg_float_t *b, int dimx, int dimy, int dimz);
void app_cir (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_ortho (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_bcc (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_fcc (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_hex (double *a, double *b, int dimx, int dimy, int dimz);
double app_cilr (double *a, double *b, double *c, int dimx, int dimy, int dimz,
               double gridhx, double gridhy, double gridhz);
double app_cilr_bcc (double *a, double *b, double *c, int dimx, int dimy, int dimz,
                   double gridhx, double gridhy, double gridhz);
double app_cilr_fcc (double *a, double *b, double *c, int dimx, int dimy, int dimz,
                   double gridhx, double gridhy, double gridhz);
double app_cilr_hex (double *a, double *b, double *c, int dimx, int dimy, int dimz,
                   double gridhx, double gridhy, double gridhz);
double app_cilr_ortho (double *a, double *b, double *c, int dimx, int dimy,
                     int dimz, double gridhx, double gridhy, double gridhz);
double app_cil (double *a, double *b, int dimx, int dimy, int dimz, double gridhx,
              double gridhy, double gridhz);
double app_cil_driver (double * a, double * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
double app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
double app_cil_fourth (double *a, double *b, int dimx, int dimy, int dimz, double gridhx,
              double gridhy, double gridhz);
double app_cil_fourth_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, double gridhx, 
              double gridhy, double gridhz);
double app_cil_sixth (double *psi, double *b, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz);
double app_cil_sixth_f (rmg_float_t *psi, rmg_float_t *b, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz);
void app_grad (double  * rho, double *wx, double *wy, double *wz, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);
//void app10_gradf (FS0_GRID * f, FP0_GRID * wx, FP0_GRID * wy, FP0_GRID * wz);
void constrain( void );
void cross_product (double *a, double *b, double *c);
double fill (STATE *states, double width, double nel, double mix,
           int num_st, int occ_flag);

void find_phase (int nldim, double *nlcdrs, double *phase_sin,
                 double *phase_cos);
//void genvpsi (double *psi, double *twovpsi, double *pvtot, 
//              double *kd, double kmag, int dimx, int dimy, int dimz);
void genvpsi_f (rmg_float_t * psi, rmg_float_t * sg_twovpsi, double * vtot, double * kd,
              double kmag, int dimx, int dimy, int dimz);
void get_nlop (void);
void get_weight (void);
void get_phase (ION *iptr, double *rtptr, int icount, int *dvec);
void get_nlop_smp (int tid);
char *get_num (char *str);


void get_xc (double * nrho, double * nrho_oppo,  double * vxc, double * exc, int xctype);

void get_zdens (STATE *states, int state, double *zvec);

/* new lda function incorporating both  1981 and 1994 monte carlo data */
void xclda(double *rho, double *vxc, double *exc);
void xclda_libxc (double * rho, double * vxc, double * exc);
void slater(double rs, double *ex, double *vx);
void pz(double rs, int iflag, double *ec, double *vc); 

/* lsda exchange correlation functional */
void xclsda_spin(double *rho, double *rho_oppo, double *vxc, double *exc);
void xclsda_spin_libxc (double * rho, double * rho_oppo, double * vxc, double * exc);
void slater_spin(double arhox, double zeta, double *ex, double *vx, double *vx_oppo);
void pz_spin(double rs, double zeta, double *ec, double *vc, double *vc_oppo);
void pz_polarized(double rs, double *ec, double *vc);
void pw_spin (double rs, double zeta, double *ec, double *vcup, double *vcdw);
void pw (double rs, int iflag, double *ec, double *vc) ;


double mu_pz (double rho);
double e_pz (double rho);
void xcgga (double *rho, double *vxc, double *exc, int flag);
void xcgga_libxc (double * rho, double * vxc, double * exc, int mode);
void xcgga_spin (double *rho, double *rho_oppo, double *vxc, double *exc, int flag);
void xcgga_spin_libxc(double * rho_up, double * rho_dw, double * vxc_up, double * exc, int mode); 


/* exchange correlation functional for PBE */

void gcxcpbe_spin(double rho_up, double rho_dw,
  double grad_up, double grad_dw, double grad, double *enxc,
  double *vxc1_up, double *vxc1_dw, double *vxc2_upup, double *vxc2_dwdw,
  double *vxc2_updw, double *vxc2_dwup);
void pbex (double rho, double grho, int iflag, double *sx, double *v1x, double *v2x);
void pbec_spin (double rho, double zeta, double grho, int iflag, double *sc, double *v1cup, double *v1cdw, double *v2c);

void gcxcpbe (double rho, double grad, double *enxc, double *vxc1, double *vxc2);

void pbec (double rho, double grad, int iflag, double *sc, double *v1c, double *v2c);



/* becke88 exchange and perdew86 correlation */

void becke88 ( double rho, double grho, double *sx, double *v1x, double *v2x );
void perdew86 ( double rho, double grho, double *sc, double *v1c, double *v2c );
void gcxbcp (double rho, double grad, double *enxc, double *vxc1, double *vxc2);
void perdew86_spin ( double rho, double zeta, double grho, double *sc, double *v1cup, double *v1cdw, double *v2c );
void becke88_spin ( double rho, double grho, double *sx, double *v1x, double *v2x );
void gcxbcp_spin (double rho_up, double rho_dw,
  double grad_up, double grad_dw, double grad, double *enxc,
  double *vxc1_up, double *vxc1_dw, double *vxc2_upup, double *vxc2_dwdw,
  double *vxc2_updw, double *vxc2_dwup);


/* PW91 exchange correlation */

void ggax(double rho, double grho, double *sx, double *v1x, double *v2x);
void ggac (double rho, double grho, double *sc, double *v1c, double *v2c);
void ggac_spin (double rho, double zeta, double grho, double *sc, double *v1cup, double *v1cdw, double *v2c);
void gcxcpw91 (double rho, double grad, double *enxc, double *vxc1, double *vxc2);
void gcxcpw91_spin(double rho_up, double rho_dw,
  double grad_up, double grad_dw, double grad, double *enxc,
  double *vxc1_up, double *vxc1_dw, double *vxc2_upup, double *vxc2_dwdw,
  double *vxc2_updw, double *vxc2_dwup);


/* BLYP exchange correlation */

void lyp ( double rs, double *ec, double *vc );
void glyp ( double rho, double grho, double *sc, double *v1c, double *v2c );
void lsd_lyp ( double rho, double zeta, double * elyp, double * valyp, double * vblyp );
void lsd_glyp ( double rhoa, double rhob, double grhoaa, double grhoab2, double grhobb, double *sc, double *v1ca, double *v2ca, double *v1cb, double *v2cb, double *v2cab );
void gcxcblyp (double rho, double grad, double *enxc, double *vxc1, double *vxc2);
void gcxcblyp_spin (double rho_up, double rho_dw,
  double grad_up, double grad_dw, double grad_updw2, double *enxc,
  double *vxc1_up, double *vxc1_dw, double *vxc2_upup, double *vxc2_dwdw,
  double *vxc2_updw, double *vxc2_dwup);



void gram (KPOINT *kpoint, double h, int numst, int maxst, int numpt,
           int maxpt);
int get_input (FILE *fh, char *id, void *dest, unsigned int flag, char *def);
void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target, int boundaryflag);
char *get_symbol (int atomic_number);
void global_sums (double *vect, int *length, MPI_Comm comm);
void init_derweight (void);
void init_derweight_s (SPECIES *sp, fftw_complex *rtptr_x,
                       fftw_complex *rtptr_y, fftw_complex *rtptr_z, int ip,
                       fftw_plan p1);
void init_derweight_p (SPECIES *sp, fftw_complex *rtptr_x,
                       fftw_complex *rtptr_y, fftw_complex *rtptr_z, int ip,
                       fftw_plan p1);
void init_derweight_d (SPECIES *sp, fftw_complex *rtptr_x,
                       fftw_complex *rtptr_y, fftw_complex *rtptr_z, int ip,
                       fftw_plan p1);
void init_fftw_wisdom (void);
void init_kbr (void);
void init_IO ( int argc, char **argv );
void init_pe ( int image );
void init_img_topo ( int dimensionality );
void init_pegrid (void);
STATE *init_states (void);
void init_weight (void);
void weight_shift_center(SPECIES * sp, fftw_complex * weptr);
void init_weight_s (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftw_plan p1);
void init_weight_p (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftw_plan p1);
void init_weight_d (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftw_plan p1);
void init_wf (STATE *states);
void init_nuc (double *vnuc, double *rhoc, double *rhocore);
void init_pos (void);
void init_sym (void);
void symmetrize_rho (double *rho);
void symforce (void);
void rmg_timings (int what, double time);
double minimage (ION *ip1, ION *ip2, double *xtal_r);
double my_crtc (void);
FILE *open_xbs_movie (char *filename);
FILE *open_restart_file (char *filename);
void xbsmovie (FILE *movie);
void ortho_half (STATE *states);
void ortho_bcc (STATE *states);
void ortho_ncpp(STATE *states);
void output_eigenvalues( STATE *states, int ikbs, int iscf );
void pack_ptos (double *sg, double *pg, int dimx, int dimy, int dimz);
void pack_ptos_f(rmg_float_t * sg, rmg_float_t * pg, int dimx, int dimy, int dimz);
void pack_stop (double *sg, double *pg, int dimx, int dimy, int dimz);
void pack_stop_f (rmg_float_t *sg, rmg_float_t *pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (double *sg, double *pg, double alpha, int dimx, int dimy, int dimz);
void pack_stop_axpy_f (rmg_float_t * sg, rmg_float_t * pg, double alpha, int dimx, int dimy, int dimz);
void pack_ptos_trade (double *sg, double *pg, int dimx, int dimy, int dimz);
void pack_vhstod (double *s, double *d, int dimx, int dimy, int dimz, int boundaryflag);
void pack_vhdtos (double *s, double *d, int dimx, int dimy, int dimz, int boundaryflag);
double radint (double *f, double *r, int n, double al);
double radint1 (double *f, double *r, double *dr_di, int n);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
void read_control (char *file);
void write_pdb (void);
int read_atom_line(char *species, double *crds, int *movable, FILE *fhand, char *tbuf, int index);
int assign_species (CONTROL * c, char *buf);

void read_pseudo (void);
double real_sum_all (double x, MPI_Comm comm);
double double_sum_all (double x, MPI_Comm comm);
double real_min_all (double x, MPI_Comm comm);


void sortpsi (STATE *states);
/*
void trade_images (double *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_f (rmg_float_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_imagesx (double *f, double *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_async (double *f, double *w, int dimx, int dimy, int dimz, int images);
void trade_imagesx_async_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images);
void trade_images1_async (double * f, int dimx, int dimy, int dimz);
*/
void trade_images1_async_f (rmg_float_t * f, int dimx, int dimy, int dimz);
void set_bc (double *mat, int dimx, int dimy, int dimz, int images, double val);
void set_bcx (double *mat, int dimx, int dimy, int dimz, int images, double val);
void getpoi_bc (double *rho, double *vh_bc, int dimx, int dimy, int dimz);
/*void trade_images2(S0_GRID *f, SS0_GRID *w);
//void trade_images2f(FS0_GRID *f, FSS0_GRID *w);
//void trade_images3(S0_GRID *f, S30_GRID *w);
//void trade_images5(S0_GRID *f, S50_GRID *w);*/
void vol_wf (STATE *states, int state, int step);
void write_avgd (double *rho);
void write_avgv (double *vh, double *vnuc);
void write_zstates (STATE *states);

void write_header (void);
void write_occ (STATE *states);
void write_force (void);
void write_timings (void);
void wvfn_residual(STATE *states);
double rand0 (long *idum);
void cgen_prolong(double coef[], double fraction, int order);
void mg_restrict_6 (double * full, double * half, int dimx, int dimy, int dimz, int grid_ratio);
void mg_prolong_6 (double * full, double * half, int dimx, int dimy, int dimz);
void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz, int grid_ratio, int order);
void gather_psi (double *tmp_psiR, double *tmp_psiI, STATE *sp, int tid);
void scatter_psi (double *tmp_psiR, double *tmp_psiI, STATE *sp, int tid);
void get_milliken (STATE *states);

void bandstructure( STATE *states, double *vxc, double *vh, double *vnuc );
void output_wave( STATE *states, int kpt, int fhand );


int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
               int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid, double * xcstart, double * ycstart, double * zcstart);

double linint (double *y, double rv, double invdr);
void my_barrier (void);


/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (double *celldm, double *A0I, double *A1I, double *A2I,
             double *OMEGAI, int *flag);
void latgen_f_ (int *ibrav, double *celldm, double *A0I, double *A1I, double *A2I,
             double *OMEGAI, int *flag);
void recips (void);
void to_cartesian (double crystal[], double cartesian[]);
void to_crystal (double crystal[], double cartesian[]);
double metric (double *crystal);

/* Md run types */
void relax (int steps, STATE *states, double *vxc, double *vh, double *vnuc,
              double *rho, double *rho_oppo, double *rhocore, double *rhoc);
void neb_relax (STATE *states, double *vxc, double *vh, double *vnuc,
              double *rho, double *rho_oppo, double *rhocore, double *rhoc);
void cdfastrlx (STATE *states, double *vxc, double *vh, double *vnuc,
                double *rho, double *rhocore, double *rhoc);
void moldyn (STATE *states, double *vxc, double *vh, double *vnuc,
             double *rho, double *rho_oppo, double *rhoc, double *rhocore);
void cholesky (double *a, int n);


/*the function for softpseudopotential*/
void aainit (int lli, int mix, int lx, int mx, int nlx, double ap[][9][9],
             int lpx[][9], int lpl[][9][9]);
double app_cil1 (double *a, double *b, int dimx, int dimy, int dimz,
               double gridhx, double gridhy, double gridhz);
double app_cil1_bcc (double *a, double *b, int dimx, int dimy, int dimz,
                   double gridhx, double gridhy, double gridhz);
double app_cil1_fcc (double *a, double *b, int dimx, int dimy, int dimz,
                   double gridhx, double gridhy, double gridhz);
double app_cil1_hex (double *a, double *b, int dimx, int dimy, int dimz,
                   double gridhx, double gridhy, double gridhz);
double app_cil1_ortho (double *a, double *b, int dimx, int dimy, int dimz,
                     double gridhx, double gridhy, double gridhz);
void app_nls (double * psiR, double * psiI, double * workR, double * workI, double *work2R, double *work2I, 
     double *sintR, double *sintI, int state, int kidx);
void app_nls_allstates (double * psiR, double * psiI, double * workR, double * workI, double *work2R, double *work2I, 
     double *Bns, double *BnsI, double *sintR, double *sintI, int kidx);
void app_nls_batch (STATE *sp, double *nv, double *ns, double *Bns, double *sintR);
void get_ddd (double *veff);
void get_nlop_d (ION *iptr, double *rtptr, int ip, int icount, int *dvec);
void get_nlop_p (ION *iptr, double *rtptr, int ip, int icount, int *dvec);
void get_nlop_s (ION *iptr, double *rtptr, int ip, int icount, int *dvec);
void get_QI (void);
void get_qqq (void);
void get_rho (STATE * states, double * rho, double * rhocore);
void get_new_rho (STATE * states, double * rho);
void get_pdos (STATE * states, double Emin, double Emax, int E_POINTS);
void mix_rho (double * new_rho, double * rho, double *rhocore, int length, int length_x, int length_y, int length_z);
void Output_rho_xsf(double *array_3d, MPI_Comm comm);
void init_psp (void);
void init_qfunct (void);
void mg_eig_state (STATE *sp, int tid, double *vtot_psi);
void mg_eig_state_f (STATE *sp, int tid, double *vtot_psi);
void ortho (STATE *states, int kpt);
double qval (int ih, int jh, double r, double invdr, double *ptpr, int *nhtol,
           int *nhtom, int *indv, double *ylm, double ap[][9][9], int lpx[][9],
           int lpl[][9][9], SPECIES *sp);
void reinit_ionic_pp (STATE * states, double * vnuc, double * rhocore, double * rhoc);
void init_pestr(void);
void read_init(char *meta);
int filename_increment(char *filename);

#if GAMMA_PT
void subdiag_gamma (STATE *states, double *vh, double *vnuc, double *vxc);
void subdiag_app_A (STATE * states, double * a_psi, double * s_psi, double * vtot_eig);
void subdiag_app_B (STATE * states, double * b_psi);
#else
void subdiag_nongamma (STATE * states, double * vh, double * vnuc, double * vxc);
void subdiag_app_A (STATE * states, double * a_psiR, double * a_psiI, double * s_psiR, double * s_psiI, double * vtot_eig);
void subdiag_app_B (STATE * states, double * b_psiR, double * b_psiI);
#endif
void init_subdiag(void);


void ylmr2 (double *r, double *ylm);
double gcutoff (double g1, double gcut, double width);
void rft1 (double cparm, double *f, double *r, double *ffil, double *rab,
           int rg_points, int lval, double dr, double width, int lrg_points);
void norm_psi1 (STATE *sp, int istate, int kpt);
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, double *cR, double *cI);
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, double cR, double cI);
double get_QnmL (int idx, int ltot, double r, SPECIES *sp);

void force (double *rho, double *rho_oppo, double *rhoc, double *vh, double *vxc, double *vnuc);


void iiforce (void);
void lforce (double *rho, double *vh);
void nlforce (double *veff);
void get_gamma (double * gammaR, int ion, int nh);
void partial_gamma (int ion, double * par_gammaR, double * par_omegaR, int nion, int nh,
                    double * newsintR_x, double * newsintR_y, double * newsintR_z,
                    double * newsintI_x, double * newsintI_y, double * newsintI_z);
void partial_betaxpsi (int ion, fftw_plan p2, double *newsintR_x,
                       double *newsintR_y, double *newsintR_z,
                       double *newsintI_x, double *newsintI_y,
                       double *newsintI_z, ION *iptr);
void partial_QI (int ion, double *QI_R, ION *iptr);
void qval_R (int ih, int jh, double r, double *x, double *qlig, double *drqlig,
             double invdr, int *nhtol, int *nhtom, int *indv, double *ylm,
             double *ylm_x, double *ylm_y, double *ylm_z, double ap[][9][9],
             int lpx[][9], int lpl[][9][9], double *Q_x, double *Q_y,
             double *Q_z, SPECIES *sp);
void ylmr2_x (double *r, double *ylm_x);
void ylmr2_y (double *r, double *ylm_y);
void ylmr2_z (double *r, double *ylm_z);
void nlccforce (double *rho, double *vxc);
double get_ve_nl (STATE *sta, int istate);
void pack_rho_ctof (double *rhoc, double *rhof);
void bspline_interp_full (double *rho, double *rho_f);
void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
void get_vxc_exc (double * nrho, double * nrho_oppo,  double * vxc, double * exc, int xctype);
void betaxpsi (STATE *states);
void betaxpsi1 (STATE *states, int kpt);
void assign_weight (SPECIES *sp, int ion, fftw_complex *beptr,
                    double *rtptr, double *Bweight);
void assign_weight2 (int nldim, int ion, double *beptr, double *rtptr);
void pack_gftoc (SPECIES *sp, fftw_complex *gwptr, fftw_complex *gbptr);
void debug_write_rho_z (double *rhoz);
void print_density_z_direction (int grid_x, int grid_y, double *density,
                                int px0_grid, int py0_grid, int pz0_grid,
                                double zside);
void get_derweight (int ion, double *beta_x, double *beta_y, double *beta_z,
                    ION *iptr, fftw_plan p2);
void partial_beta_fdiff (fftw_complex *beptr, int nldim, double *beta_x,
                         double *beta_y, double *beta_z);

void mulliken (STATE *states);
double ylm(int l, double *r);
int listlen (FILE * fh, char *id);
void print_matrix(double *b, int n, int ldb);
void sl_init(int *ictxt, int size);
void sl_exit(int ictxt);
void set_desca(int *desca, int *ictxt, int size);
void distribute_mat(int *desca, double *bigmat, double *dismat, int *size);
void matinit(int *desca, double *dismat, double *globmat, int size);
int matsum_packbuffer(int row, int col, double *buffer, double *globmat, int size);
void reduce_and_dist_matrix(int n, double *global_matrix, double *dist_matrix, double *work);
void init_efield (double * vnuc);
void pulay_rho(int step, int N, int N_x, int N_y, int N_z, double *rho_new, double *rho_old, int NsavedSteps, double ***hist, double ***rhist, int special_metric, double weight);
void mg_restrict_4th (double * full, double * half, int dimx, int dimy, int dimz);
void mix_betaxpsi (int mix);
int claim_ion (double *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid);
int find_grid_1d_owner(int igrid, int tgrid, int pgrid);
void find_grid_owner(int igridx, int igridy, int igridz, int nxgrid, int nygrid, int nzgrid, int *xnode, int *ynode, int *znode);
int is_loop_over_states(void);
void fastrelax (double *dt, double dt_max, double dt_inc, double dt_dec, int n_min, int *n_count);
void fire (double *step, double step_max, double f_inc, double f_dec, int n_min, int *n_count );
void quick_min (void);
int int_sum_all (int x, MPI_Comm comm);
void move_ions (double dt);
void get_extrapolation_constants (double *alpha, double *beta);
void lcao_init (void);
void init_atomic_rho_wf (void);
void lcao_get_rho (double * arho_f);
void lcao_get_awave (double *psi, ION *iptr, int awave_idx, int l, int m, double coeff);
void lcao_get_psi (STATE * states);
double mask_function (double x);
void apply_mask_function (double *f, double * r, int rg_points, double rmax, double offset);
void filter_potential (double *potential, double *r, int rg_points, double rmax, double offset, double parm, double* potential_lgrid, 
	double *rab, int l_value, double dr, double  gwidth, int lgrid_points, double rcut, double rwidth, double * drpotential_lgrid);
int test_overlap (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid);
void RMG_MPI_trade(double *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
void RMG_MPI_trade_f(rmg_float_t *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
void init_trade_imagesx_async(void);
void  get_rho_oppo (double * rho, double * rho_oppo);
void get_opposite_eigvals (STATE * states);
void get_opposite_occupancies (STATE * states);

void allocate_states();
#if __cplusplus
}
#endif

