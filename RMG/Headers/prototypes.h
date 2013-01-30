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


/* Function prototypes */
REAL app_del2c (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                REAL gridhx, REAL gridhy, REAL gridhz);
void app6_del2 (REAL *rho, REAL *work, int dimx, int dimy, int dimz,
                REAL gridhx, REAL gridhy, REAL gridhz);
void app_smooth (REAL *f, REAL *work, int dimx, int dimy, int dimz);
void app_smooth1 (REAL *f, REAL *work, int dimx, int dimy, int dimz);
void app_cir_driver (REAL *a, REAL *b, int dimx, int dimy, int dimz, int order);
void app_cir_fourth (REAL *a, REAL *b, int dimx, int dimy, int dimz);
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
REAL app_cil_driver (REAL * a, REAL * b, int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy, REAL gridhz, int order);
REAL app_cil_fourth (REAL *a, REAL *b, int dimx, int dimy, int dimz, REAL gridhx,
              REAL gridhy, REAL gridhz);
REAL app_cil_sixth (REAL *psi, REAL *b, int dimx, int dimy, int dimz,
                    REAL gridhx, REAL gridhy, REAL gridhz);
void app_grad (REAL  * rho, REAL *wx, REAL *wy, REAL *wz, int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy, REAL gridhz);
//void app10_gradf (FS0_GRID * f, FP0_GRID * wx, FP0_GRID * wy, FP0_GRID * wz);
void constrain( void );
void corlyp (REAL *dp, REAL *dm, REAL *dp1, REAL *dm1, REAL *dp2, REAL *dm2, REAL *ec,
             REAL *vcp0, REAL *vcm0, int *ndm);
void cross_product (REAL *a, REAL *b, REAL *c);
void destroy_fftw_wisdom (void);
void eval_residual (REAL *mat, REAL *f_mat, int dimx, int dimy, int dimz,
                    REAL gridhx, REAL gridhy, REAL gridhz, REAL *res);
void solv_pois (REAL *vmat, REAL *fmat, REAL *work,
                int dimx, int dimy, int dimz, REAL gridhx,
                REAL gridhy, REAL gridhz, REAL step, REAL k);
REAL fill (STATE *states, REAL width, REAL nel, REAL mix,
           int num_st, int occ_flag);

void find_phase (int nldim, REAL *nlcdrs, REAL *phase_sin,
                 REAL *phase_cos);
void finish_release_mem(STATE *states);
void genvpsi (REAL *psi, REAL *twovpsi, REAL *pvtot, REAL *pvnl,
              REAL *kd, REAL kmag, int dimx, int dimy, int dimz);
void get_nlop (void);
void get_weight (void);
void get_phase (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_smp (int tid);
char *get_num (char *str);
void get_te (REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc, REAL *vh, REAL *vxc,
             STATE *states);

void get_vxc (REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *vxc);

void get_xc (REAL * nrho, REAL * nrho_oppo,  REAL * vxc, REAL * exc, int xctype);

void get_zdens (STATE *states, int state, REAL *zvec);
void xclda_pz81 (REAL *rho, REAL *vxc);
void exclda_pz81 (REAL *rho, REAL *exc);
void get_psi_overlaps(REAL *psi_array, REAL *overlap, int numst, int maxst, int numpt, int maxpt);

/* new lda function incorporating both  1981 and 1994 monte carlo data */
void xclda(REAL *rho, REAL *vxc, REAL *exc);
void xclda_libxc (REAL * rho, REAL * vxc, REAL * exc);
void slater(REAL rs, REAL *ex, REAL *vx);
void pz(REAL rs, int iflag, REAL *ec, REAL *vc); 

/* lsda exchange correlation functional */
void xclsda_spin(REAL *rho, REAL *rho_oppo, REAL *vxc, REAL *exc);
void xclsda_spin_libxc (REAL * rho, REAL * rho_oppo, REAL * vxc, REAL * exc);
void slater_spin(REAL arhox, REAL zeta, REAL *ex, REAL *vx, REAL *vx_oppo);
void pz_spin(REAL rs, REAL zeta, REAL *ec, REAL *vc, REAL *vc_oppo);
void pz_polarized(REAL rs, REAL *ec, REAL *vc);
void pw_spin (REAL rs, REAL zeta, REAL *ec, REAL *vcup, REAL *vcdw);
void pw (REAL rs, int iflag, REAL *ec, REAL *vc) ;


double mu_pz (double rho);
double e_pz (double rho);
void xcgga (REAL *rho, REAL *vxc, REAL *exc, int flag);
void xcgga_libxc (REAL * rho, REAL * vxc, REAL * exc, int mode);
void xcgga_spin (REAL *rho, REAL *rho_oppo, REAL *vxc, REAL *exc, int flag);
void xcgga_spin_libxc(REAL * rho_up, REAL * rho_dw, REAL * vxc_up, REAL * exc, int mode); 


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
void get_vh (REAL * rho, REAL * rhoc, REAL * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, REAL rms_target);
char *get_symbol (int atomic_number);
void global_sums (REAL *vect, int *length, MPI_Comm comm);
void init (REAL *vh, REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc, STATE *states,
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
void init_pe ( int image );
void init_img_topo ( int dimensionality );
void init_pegrid (void);
STATE *init_states (void);
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
#if GAMMA_POINT
void symmetrize_rho (FP0_GRID *rho);
#endif
void symforce (void);
int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype);
void mgrid_solv (REAL *v_mat, REAL *f_mat, REAL *work,
                 int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy,
                 REAL gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, REAL step, REAL k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);
void rmg_timings (int what, REAL time);
REAL minimage (ION *ip1, ION *ip2, REAL *xtal_r);
REAL my_crtc (void);
FILE *open_xbs_movie (char *filename);
FILE *open_restart_file (char *filename);
void xbsmovie (FILE *movie);
void ortho_half (STATE *states);
void ortho_bcc (STATE *states);
void ortho_ncpp(STATE *states);
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
void read_control (char *file);
void write_pdb (void);
int read_atom_line(char *species, REAL *crds, int *movable, FILE *fhand, char *tbuf, int index);
int assign_species (CONTROL * c, char *buf);
void read_data (char *name, REAL *vh, REAL *rho, REAL *vxc,
                STATE *states);

void read_pseudo (void);
REAL real_sum_all (REAL x, MPI_Comm comm);
REAL real_min_all (REAL x, MPI_Comm comm);


void reset_timers (void);
void sortpsi (STATE *states);
void trade_images (REAL *mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_imagesx (REAL *f, REAL *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_async (REAL *f, REAL *w, int dimx, int dimy, int dimz, int images);
void trade_images1_async (REAL * f, int dimx, int dimy, int dimz);
void set_bc (REAL *mat, int dimx, int dimy, int dimz, int images, REAL val);
void set_bcx (REAL *mat, int dimx, int dimy, int dimz, int images, REAL val);
void getpoi_bc (REAL *rho, REAL *vh_bc, int dimx, int dimy, int dimz);
/*void trade_images2(S0_GRID *f, SS0_GRID *w);
//void trade_images2f(FS0_GRID *f, FSS0_GRID *w);
//void trade_images3(S0_GRID *f, S30_GRID *w);
//void trade_images5(S0_GRID *f, S50_GRID *w);*/
void vol_wf (STATE *states, int state, int step);
void write_avgd (REAL *rho);
void write_avgv (REAL *vh, REAL *vnuc);
void write_zstates (STATE *states);
void write_data (int fhand, REAL *vh, REAL *rho, REAL *rho_oppo, REAL *vxc,
                 STATE *states);

void write_header (void);
void write_occ (STATE *states);
void write_force (void);
void write_timings (void);
void wvfn_residual(STATE *states);
REAL rand0 (long *idum);
void cgen_prolong(REAL coef[], REAL fraction, int order);
void mg_restrict (REAL *full, REAL *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_restrict_6 (REAL * full, REAL * half, int dimx, int dimy, int dimz, int grid_ratio);
void mg_prolong (REAL *full, REAL *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_prolong_6 (REAL * full, REAL * half, int dimx, int dimy, int dimz);
void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz, int grid_ratio, int order);
void gather_psi (REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid);
void scatter_psi (REAL *tmp_psiR, REAL *tmp_psiI, STATE *sp, int tid);
void get_milliken (STATE *states);

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

int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
               int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid, REAL * xcstart, REAL * ycstart, REAL * zcstart,
               int pxoffset, int pyoffset, int pzoffset);

REAL linint (REAL *y, REAL rv, REAL invdr);
void my_barrier (void);


/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, REAL *celldm, REAL *A0I, REAL *A1I, REAL *A2I,
             REAL *OMEGAI, int *flag);
void recips (void);
void to_cartesian (REAL crystal[], REAL cartesian[]);
void to_crystal (REAL crystal[], REAL cartesian[]);
REAL metric (REAL *crystal);

/* Md run types */
bool quench (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc, REAL *rho,
             REAL *rho_oppo, REAL *rhocore, REAL *rhoc);
void relax (int steps, STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
              REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc);
void neb_relax (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
              REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc);
void cdfastrlx (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
                REAL *rho, REAL *rhocore, REAL *rhoc);
void moldyn (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
             REAL *rho, REAL *rho_oppo, REAL *rhoc, REAL *rhocore);
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
void app_nls (REAL * psiR, REAL * psiI, REAL * workR, REAL * workI, REAL *work2R, REAL *work2I, REAL *sintR, REAL *sintI, int state,
                 int kidx);
void get_ddd (REAL *veff);
void get_nlop_d (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_p (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_nlop_s (ION *iptr, REAL *rtptr, int ip, int icount, int *dvec);
void get_QI (void);
void get_qqq (void);
void get_rho (STATE * states, REAL * rho, REAL * rhocore);
void get_new_rho (STATE * states, REAL * rho);
void get_pdos (STATE * states, REAL Emin, REAL Emax, int E_POINTS);
void mix_rho (REAL * new_rho, REAL * rho, REAL *rhocore, int length, int length_x, int length_y, int length_z);
void init_psp (void);
void init_qfunct (void);
void mg_eig_state (STATE *sp, int tid, REAL *vtot_psi);
void ortho (STATE *states, int kpt);
REAL qval (int ih, int jh, REAL r, REAL invdr, REAL *ptpr, int *nhtol,
           int *nhtom, int *indv, REAL *ylm, REAL ap[][9][9], int lpx[][9],
           int lpl[][9][9], SPECIES *sp);
bool scf (STATE *states, REAL *vxc, REAL *vh, REAL *vnuc,
          REAL *rho, REAL *rho_oppo, REAL *rhocore, REAL *rhoc);
void reinit_ionic_pp (STATE * states, REAL * vnuc, REAL * rhocore, REAL * rhoc);

#if GAMMA_PT
void subdiag_gamma (STATE *states, REAL *vh, REAL *vnuc, REAL *vxc);
void subdiag_app_A (STATE * states, REAL * a_psi, REAL * s_psi, REAL * vtot_eig);
void subdiag_app_B (STATE * states, REAL * b_psi);
void subdiag_app_B_one (STATE *sp, REAL * b_psi);
void subdiag_app_A_one (STATE *sp, REAL * a_psi, REAL * s_psi, REAL * vtot_eig);
#else
void subdiag_nongamma (STATE * states, REAL * vh, REAL * vnuc, REAL * vxc);
void subdiag_app_A (STATE * states, REAL * a_psiR, REAL * a_psiI, REAL * s_psiR, REAL * s_psiI, REAL * vtot_eig);
void subdiag_app_B (STATE * states, REAL * b_psiR, REAL * b_psiI);
#endif
void betaxpsi1_calculate_one(STATE *st, int ion, int nion, REAL *sintR, REAL *sintI, int kpt);
void init_subdiag(void);


void ylmr2 (double *r, double *ylm);
REAL gcutoff (REAL g1, REAL gcut, REAL width);
void rft1 (REAL cparm, REAL *f, REAL *r, REAL *ffil, REAL *rab,
           int rg_points, int lval, REAL dr, REAL width, int lrg_points);
void norm_psi1 (STATE *sp, int istate, int kpt);
#if FAST_ORTHO
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL *cR, REAL *cI, REAL *Oij);
#else
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL *cR, REAL *cI);
#endif
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL cR, REAL cI);
REAL get_QnmL (int idx, int ltot, REAL r, SPECIES *sp);

void force (REAL *rho, REAL *rho_oppo, REAL *rhoc, REAL *vh, REAL *vxc, REAL *vnuc,
            STATE *states);


void iiforce (void);
void lforce (REAL *rho, REAL *vh);
void nlforce (REAL *veff);
void get_gamma (REAL * gammaR, int ion, int nh);
void partial_gamma (int ion, REAL * par_gammaR, REAL * par_omegaR, int nion, int nh,
                    REAL * newsintR_x, REAL * newsintR_y, REAL * newsintR_z,
                    REAL * newsintI_x, REAL * newsintI_y, REAL * newsintI_z);
void partial_betaxpsi (int ion, fftwnd_plan p2, REAL *newsintR_x,
                       REAL *newsintR_y, REAL *newsintR_z,
                       REAL *newsintI_x, REAL *newsintI_y,
                       REAL *newsintI_z, ION *iptr);
void nlforce_par_Q (REAL *veff, REAL *gamma, int ion, ION *iptr, int nh,
                     REAL *forces);
void nlforce_par_gamma (REAL * par_gamma, int ion, int nh, REAL *force);
void nlforce_par_omega (REAL * par_omega, int ion, int nh, REAL *force);
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
void get_vtot_psi (REAL * vtot_psi, REAL * vtot, int grid_ratio);
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
void print_matrix(double *b, int n, int ldb);
void sl_init(int *ictxt, int size);
void sl_exit(int ictxt);
void set_desca(int *desca, int *ictxt, int size);
void distribute_mat(int *desca, double *bigmat, double *dismat, int *size);
void matinit(int *desca, double *dismat, double *globmat, int size);
int matsum_packbuffer(int row, int col, double *buffer, double *globmat, int size);
void reduce_and_dist_matrix(int n, REAL *global_matrix, REAL *dist_matrix, REAL *work);
void print_distribute_mat(double *dismat, int *desca, int size);
void init_efield (REAL * vnuc);
void pulay_rho(int step, int N, int N_x, int N_y, int N_z, double *rho_new, double *rho_old, int NsavedSteps, REAL ***hist, REAL ***rhist, int special_metric, REAL weight);
void mg_restrict_4th (REAL * full, REAL * half, int dimx, int dimy, int dimz);
void mix_betaxpsi (int mix);
void mix_betaxpsi1 (STATE *sp);
int claim_ion (REAL *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid);
int find_grid_1d_owner(int igrid, int tgrid, int pgrid);
int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset);
int is_loop_over_states(void);
int rmg_is_open_mp_safe(void);
void RMG_set_omp_parallel(void);
void RMG_set_omp_single(void);
void Papi_init_omp_threads(int ithread);
void Papi_finalize_omp_threads(int ithread);
void fastrelax (REAL *dt, REAL dt_max, REAL dt_inc, REAL dt_dec, int n_min, int *n_count);
void fire (REAL *step, REAL step_max, REAL f_inc, REAL f_dec, int n_min, int *n_count );
void quick_min (void);
int int_sum_all (int x, MPI_Comm comm);
void move_ions (REAL dt);
void get_extrapolation_constants (REAL *alpha, REAL *beta);
void lcao_init (void);
void init_atomic_rho_wf (void);
void lcao_get_rho (REAL * arho_f);
void lcao_get_awave (REAL *psi, ION *iptr, int awave_idx, int l, int m, double coeff);
void lcao_get_psi (STATE * states);
REAL mask_function (REAL x);
void apply_mask_function (REAL *f, REAL * r, int rg_points, REAL rmax, REAL offset);
void filter_potential (REAL *potential, REAL *r, int rg_points, REAL rmax, REAL offset, REAL parm, REAL* potential_lgrid, 
	REAL *rab, int l_value, REAL dr, REAL  gwidth, int lgrid_points, REAL rcut, REAL rwidth, REAL * drpotential_lgrid);
int test_overlap (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid);
void RMG_MPI_trade(REAL *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
void init_trade_imagesx_async(void);
void  get_rho_oppo (REAL * rho, REAL * rho_oppo);
void get_opposite_eigvals (STATE * states);
void get_opposite_occupancies (STATE * states);

#if GPU_ENABLED
void init_gpu (void);
void finalize_gpu (void);
void app_cil_sixth_gpu(const double *psi,
                                                double *b,
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const double gridhx,
                                                const double gridhy,
                                                const double gridhz,
                                                const double xside,
                                                const double yside,
                                                const double zside,
                                                cudaStream_t cstream);
void app_cir_sixth_gpu(const double *psi,
                                                double *b,
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                cudaStream_t cstream);

#endif
