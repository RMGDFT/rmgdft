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

void CompatRmgTimerPrint(const char *outfile, int steps);
void betaxpsi1_calculate_one(STATE *st, int ion, int nion, rmg_double_t *sintR, rmg_double_t *sintI, int kpt, rmg_double_t *weiptr_base);
void get_te (rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *rhoc, rmg_double_t *vh, rmg_double_t *vxc,
             STATE *states, int ii_flag);
void get_vxc (rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *vxc);
void init (rmg_double_t *vh, rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *rhoc, STATE *states,
           rmg_double_t *vnuc, rmg_double_t *vxc);
void mg_eig_state_driver (STATE * sp, int tid, rmg_double_t * vtot_psi);
void nlforce_par_gamma (rmg_double_t * par_gamma, int ion, int nh, rmg_double_t *force);
void nlforce_par_omega (rmg_double_t * par_omega, int ion, int nh, rmg_double_t *force);
void nlforce_par_Q (rmg_double_t *veff, rmg_double_t *gamma, int ion, ION *iptr, int nh,
                     rmg_double_t *forces);
bool quench (STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc, rmg_double_t *rho,
             rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *rhoc);
void read_data (char *name, rmg_double_t *vh, rmg_double_t *rho, rmg_double_t *vxc,
                STATE *states);
bool scf (STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc,
          rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *rhoc);
void subdiag_app_B_one (STATE *sp, rmg_double_t * b_psi);
void subdiag_app_A_one (STATE *sp, rmg_double_t * a_psi, rmg_double_t * s_psi, rmg_double_t * vtot_eig);
void subdiag_app_AB_one (STATE *sp, rmg_double_t * a_psi, rmg_double_t * b_psi, rmg_double_t * vtot_eig_s);
void subdiag_app_AB (STATE * states, rmg_double_t * a_psi, rmg_double_t * b_psi, rmg_double_t * vtot_eig);
void write_data (int fhand, rmg_double_t *vh, rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *vxc,
                 STATE *states);
#if 0
/* Function prototypes */
void app6_del2 (rmg_double_t *rho, rmg_double_t *work, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void app_smooth (rmg_double_t *f, rmg_double_t *work, int dimx, int dimy, int dimz);
void app_smooth_f (rmg_float_t * f, rmg_float_t * work, int dimx, int dimy, int dimz);
void app_smooth1 (rmg_double_t *f, rmg_double_t *work, int dimx, int dimy, int dimz);
void app_smooth1_f (rmg_float_t *f, rmg_float_t *work, int dimx, int dimy, int dimz);
void app_cir_driver (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz, int order);
void app_cir_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, int order);
void app_cir_fourth (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_fourth_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz);
void app_cir_sixth (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_sixth_f (rmg_float_t *a, rmg_float_t *b, int dimx, int dimy, int dimz);
void app_cir (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_ortho (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_bcc (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_fcc (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_hex (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
rmg_double_t app_cilr (rmg_double_t *a, rmg_double_t *b, rmg_double_t *c, int dimx, int dimy, int dimz,
               rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cilr_bcc (rmg_double_t *a, rmg_double_t *b, rmg_double_t *c, int dimx, int dimy, int dimz,
                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cilr_fcc (rmg_double_t *a, rmg_double_t *b, rmg_double_t *c, int dimx, int dimy, int dimz,
                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cilr_hex (rmg_double_t *a, rmg_double_t *b, rmg_double_t *c, int dimx, int dimy, int dimz,
                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cilr_ortho (rmg_double_t *a, rmg_double_t *b, rmg_double_t *c, int dimx, int dimy,
                     int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx,
              rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil_driver (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order);
rmg_double_t app_cil_driver_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, int order);
rmg_double_t app_cil_fourth (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx,
              rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil_fourth_f (rmg_float_t * a, rmg_float_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, 
              rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil_sixth (rmg_double_t *psi, rmg_double_t *b, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil_sixth_f (rmg_float_t *psi, rmg_float_t *b, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void app_grad (rmg_double_t  * rho, rmg_double_t *wx, rmg_double_t *wy, rmg_double_t *wz, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
//void app10_gradf (FS0_GRID * f, FP0_GRID * wx, FP0_GRID * wy, FP0_GRID * wz);
void constrain( void );
void cross_product (rmg_double_t *a, rmg_double_t *b, rmg_double_t *c);
void destroy_fftw_wisdom (void);
rmg_double_t fill (STATE *states, rmg_double_t width, rmg_double_t nel, rmg_double_t mix,
           int num_st, int occ_flag);

void find_phase (int nldim, rmg_double_t *nlcdrs, rmg_double_t *phase_sin,
                 rmg_double_t *phase_cos);
void finish_release_mem(STATE *states);
void genvpsi (rmg_double_t *psi, rmg_double_t *twovpsi, rmg_double_t *pvtot, rmg_double_t *pvnl,
              rmg_double_t *kd, rmg_double_t kmag, int dimx, int dimy, int dimz);
void genvpsi_f (rmg_float_t * psi, rmg_float_t * sg_twovpsi, rmg_double_t * vtot, rmg_double_t * vnl, rmg_double_t * kd,
              rmg_double_t kmag, int dimx, int dimy, int dimz);
void get_nlop (void);
void get_weight (void);
void get_phase (ION *iptr, rmg_double_t *rtptr, int ip, int icount, int *dvec);
void get_nlop_smp (int tid);
char *get_num (char *str);


void get_xc (rmg_double_t * nrho, rmg_double_t * nrho_oppo,  rmg_double_t * vxc, rmg_double_t * exc, int xctype);

void get_zdens (STATE *states, int state, rmg_double_t *zvec);
void get_psi_overlaps(rmg_double_t *psi_array, rmg_double_t *overlap, int numst, int maxst, int numpt, int maxpt);

/* new lda function incorporating both  1981 and 1994 monte carlo data */
void xclda(rmg_double_t *rho, rmg_double_t *vxc, rmg_double_t *exc);
void xclda_libxc (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc);
void slater(rmg_double_t rs, rmg_double_t *ex, rmg_double_t *vx);
void pz(rmg_double_t rs, int iflag, rmg_double_t *ec, rmg_double_t *vc); 

/* lsda exchange correlation functional */
void xclsda_spin(rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *vxc, rmg_double_t *exc);
void xclsda_spin_libxc (rmg_double_t * rho, rmg_double_t * rho_oppo, rmg_double_t * vxc, rmg_double_t * exc);
void slater_spin(rmg_double_t arhox, rmg_double_t zeta, rmg_double_t *ex, rmg_double_t *vx, rmg_double_t *vx_oppo);
void pz_spin(rmg_double_t rs, rmg_double_t zeta, rmg_double_t *ec, rmg_double_t *vc, rmg_double_t *vc_oppo);
void pz_polarized(rmg_double_t rs, rmg_double_t *ec, rmg_double_t *vc);
void pw_spin (rmg_double_t rs, rmg_double_t zeta, rmg_double_t *ec, rmg_double_t *vcup, rmg_double_t *vcdw);
void pw (rmg_double_t rs, int iflag, rmg_double_t *ec, rmg_double_t *vc) ;


double mu_pz (double rho);
double e_pz (double rho);
void xcgga (rmg_double_t *rho, rmg_double_t *vxc, rmg_double_t *exc, int flag);
void xcgga_libxc (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc, int mode);
void xcgga_spin (rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *vxc, rmg_double_t *exc, int flag);
void xcgga_spin_libxc(rmg_double_t * rho_up, rmg_double_t * rho_dw, rmg_double_t * vxc_up, rmg_double_t * exc, int mode); 


/* exchange correlation functional for PBE */

void gcxcpbe_spin(rmg_double_t rho_up, rmg_double_t rho_dw,
  rmg_double_t grad_up, rmg_double_t grad_dw, rmg_double_t grad, rmg_double_t *enxc,
  rmg_double_t *vxc1_up, rmg_double_t *vxc1_dw, rmg_double_t *vxc2_upup, rmg_double_t *vxc2_dwdw,
  rmg_double_t *vxc2_updw, rmg_double_t *vxc2_dwup);
void pbex (rmg_double_t rho, rmg_double_t grho, int iflag, rmg_double_t *sx, rmg_double_t *v1x, rmg_double_t *v2x);
void pbec_spin (rmg_double_t rho, rmg_double_t zeta, rmg_double_t grho, int iflag, rmg_double_t *sc, rmg_double_t *v1cup, rmg_double_t *v1cdw, rmg_double_t *v2c);

void gcxcpbe (rmg_double_t rho, rmg_double_t grad, rmg_double_t *enxc, rmg_double_t *vxc1, rmg_double_t *vxc2);

void pbec (rmg_double_t rho, rmg_double_t grad, int iflag, rmg_double_t *sc, rmg_double_t *v1c, rmg_double_t *v2c);



/* becke88 exchange and perdew86 correlation */

void becke88 ( rmg_double_t rho, rmg_double_t grho, rmg_double_t *sx, rmg_double_t *v1x, rmg_double_t *v2x );
void perdew86 ( rmg_double_t rho, rmg_double_t grho, rmg_double_t *sc, rmg_double_t *v1c, rmg_double_t *v2c );
void gcxbcp (rmg_double_t rho, rmg_double_t grad, rmg_double_t *enxc, rmg_double_t *vxc1, rmg_double_t *vxc2);
void perdew86_spin ( rmg_double_t rho, rmg_double_t zeta, rmg_double_t grho, rmg_double_t *sc, rmg_double_t *v1cup, rmg_double_t *v1cdw, rmg_double_t *v2c );
void becke88_spin ( rmg_double_t rho, rmg_double_t grho, rmg_double_t *sx, rmg_double_t *v1x, rmg_double_t *v2x );
void gcxbcp_spin (rmg_double_t rho_up, rmg_double_t rho_dw,
  rmg_double_t grad_up, rmg_double_t grad_dw, rmg_double_t grad, rmg_double_t *enxc,
  rmg_double_t *vxc1_up, rmg_double_t *vxc1_dw, rmg_double_t *vxc2_upup, rmg_double_t *vxc2_dwdw,
  rmg_double_t *vxc2_updw, rmg_double_t *vxc2_dwup);


/* PW91 exchange correlation */

void ggax(rmg_double_t rho, rmg_double_t grho, rmg_double_t *sx, rmg_double_t *v1x, rmg_double_t *v2x);
void ggac (rmg_double_t rho, rmg_double_t grho, rmg_double_t *sc, rmg_double_t *v1c, rmg_double_t *v2c);
void ggac_spin (rmg_double_t rho, rmg_double_t zeta, rmg_double_t grho, rmg_double_t *sc, rmg_double_t *v1cup, rmg_double_t *v1cdw, rmg_double_t *v2c);
void gcxcpw91 (rmg_double_t rho, rmg_double_t grad, rmg_double_t *enxc, rmg_double_t *vxc1, rmg_double_t *vxc2);
void gcxcpw91_spin(rmg_double_t rho_up, rmg_double_t rho_dw,
  rmg_double_t grad_up, rmg_double_t grad_dw, rmg_double_t grad, rmg_double_t *enxc,
  rmg_double_t *vxc1_up, rmg_double_t *vxc1_dw, rmg_double_t *vxc2_upup, rmg_double_t *vxc2_dwdw,
  rmg_double_t *vxc2_updw, rmg_double_t *vxc2_dwup);


/* BLYP exchange correlation */

void lyp ( rmg_double_t rs, rmg_double_t *ec, rmg_double_t *vc );
void glyp ( rmg_double_t rho, rmg_double_t grho, rmg_double_t *sc, rmg_double_t *v1c, rmg_double_t *v2c );
void lsd_lyp ( rmg_double_t rho, rmg_double_t zeta, rmg_double_t * elyp, rmg_double_t * valyp, rmg_double_t * vblyp );
void lsd_glyp ( rmg_double_t rhoa, rmg_double_t rhob, rmg_double_t grhoaa, rmg_double_t grhoab2, rmg_double_t grhobb, rmg_double_t *sc, rmg_double_t *v1ca, rmg_double_t *v2ca, rmg_double_t *v1cb, rmg_double_t *v2cb, rmg_double_t *v2cab );
void gcxcblyp (rmg_double_t rho, rmg_double_t grad, rmg_double_t *enxc, rmg_double_t *vxc1, rmg_double_t *vxc2);
void gcxcblyp_spin (rmg_double_t rho_up, rmg_double_t rho_dw,
  rmg_double_t grad_up, rmg_double_t grad_dw, rmg_double_t grad_updw2, rmg_double_t *enxc,
  rmg_double_t *vxc1_up, rmg_double_t *vxc1_dw, rmg_double_t *vxc2_upup, rmg_double_t *vxc2_dwdw,
  rmg_double_t *vxc2_updw, rmg_double_t *vxc2_dwup);



void gram (KPOINT *kpoint, rmg_double_t h, int numst, int maxst, int numpt,
           int maxpt);
int get_input (FILE *fh, char *id, void *dest, unsigned int flag, char *def);
rmg_double_t get_ke (STATE *sp, int tid);
void get_vh (rmg_double_t * rho, rmg_double_t * rhoc, rmg_double_t * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, rmg_double_t rms_target);
char *get_symbol (int atomic_number);
void global_sums (rmg_double_t *vect, int *length, MPI_Comm comm);
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
void init_weight_s (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftw_plan p1);
void init_weight_p (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftw_plan p1);
void init_weight_d (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftw_plan p1);
void init_wf (STATE *states);
void init_nuc (rmg_double_t *vnuc, rmg_double_t *rhoc, rmg_double_t *rhocore);
void init_pos (void);
void init_sym (void);
#if GAMMA_POINT
void symmetrize_rho (FP0_GRID *rho);
#endif
void symforce (void);
rmg_double_t minimage (ION *ip1, ION *ip2, rmg_double_t *xtal_r);
rmg_double_t my_crtc (void);
FILE *open_xbs_movie (char *filename);
FILE *open_restart_file (char *filename);
void xbsmovie (FILE *movie);
void ortho_half (STATE *states);
void ortho_bcc (STATE *states);
void ortho_ncpp(STATE *states);
void output_eigenvalues( STATE *states, int ikbs, int iscf );
void pack_ptos (rmg_double_t *sg, rmg_double_t *pg, int dimx, int dimy, int dimz);
void pack_ptos_f(rmg_float_t * sg, rmg_float_t * pg, int dimx, int dimy, int dimz);
void pack_stop (rmg_double_t *sg, rmg_double_t *pg, int dimx, int dimy, int dimz);
void pack_stop_f (rmg_float_t *sg, rmg_float_t *pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (rmg_double_t *sg, rmg_double_t *pg, rmg_double_t alpha, int dimx, int dimy, int dimz);
void pack_stop_axpy_f (rmg_float_t * sg, rmg_float_t * pg, rmg_double_t alpha, int dimx, int dimy, int dimz);
void pack_ptos_trade (rmg_double_t *sg, rmg_double_t *pg, int dimx, int dimy, int dimz);
void pack_vhstod (rmg_double_t *s, rmg_double_t *d, int dimx, int dimy, int dimz);
void pack_vhdtos (rmg_double_t *s, rmg_double_t *d, int dimx, int dimy, int dimz);
double radint (double *f, double *r, int n, double al);
rmg_double_t radint1 (rmg_double_t *f, rmg_double_t *r, rmg_double_t *dr_di, int n);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
void read_control (char *file);
void write_pdb (void);
int read_atom_line(char *species, rmg_double_t *crds, int *movable, FILE *fhand, char *tbuf, int index);
int assign_species (CONTROL * c, char *buf);

void read_pseudo (void);
rmg_double_t real_sum_all (rmg_double_t x, MPI_Comm comm);
rmg_double_t real_min_all (rmg_double_t x, MPI_Comm comm);


void reset_timers (void);
void sortpsi (STATE *states);
/*
void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_f (rmg_float_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_imagesx (rmg_double_t *f, rmg_double_t *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_async (rmg_double_t *f, rmg_double_t *w, int dimx, int dimy, int dimz, int images);
void trade_imagesx_async_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images);
void trade_images1_async (rmg_double_t * f, int dimx, int dimy, int dimz);
*/
void trade_images1_async_f (rmg_float_t * f, int dimx, int dimy, int dimz);
void set_bc (rmg_double_t *mat, int dimx, int dimy, int dimz, int images, rmg_double_t val);
void set_bcx (rmg_double_t *mat, int dimx, int dimy, int dimz, int images, rmg_double_t val);
void getpoi_bc (rmg_double_t *rho, rmg_double_t *vh_bc, int dimx, int dimy, int dimz);
/*void trade_images2(S0_GRID *f, SS0_GRID *w);
//void trade_images2f(FS0_GRID *f, FSS0_GRID *w);
//void trade_images3(S0_GRID *f, S30_GRID *w);
//void trade_images5(S0_GRID *f, S50_GRID *w);*/
void vol_wf (STATE *states, int state, int step);
void write_avgd (rmg_double_t *rho);
void write_avgv (rmg_double_t *vh, rmg_double_t *vnuc);
void write_zstates (STATE *states);

void write_header (void);
void write_occ (STATE *states);
void write_force (void);
void write_timings (void);
void wvfn_residual(STATE *states);
rmg_double_t rand0 (long *idum);
void cgen_prolong(rmg_double_t coef[], rmg_double_t fraction, int order);
void mg_restrict_6 (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz, int grid_ratio);
void mg_prolong_6 (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz);
void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz, int grid_ratio, int order);
void gather_psi (rmg_double_t *tmp_psiR, rmg_double_t *tmp_psiI, STATE *sp, int tid);
void scatter_psi (rmg_double_t *tmp_psiR, rmg_double_t *tmp_psiI, STATE *sp, int tid);
void get_milliken (STATE *states);

void bandstructure( STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc );
void output_wave( STATE *states, int kpt, int fhand );

int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
               int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid, rmg_double_t * xcstart, rmg_double_t * ycstart, rmg_double_t * zcstart);

rmg_double_t linint (rmg_double_t *y, rmg_double_t rv, rmg_double_t invdr);
void my_barrier (void);


/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, rmg_double_t *celldm, rmg_double_t *A0I, rmg_double_t *A1I, rmg_double_t *A2I,
             rmg_double_t *OMEGAI, int *flag);
void recips (void);
void to_cartesian (rmg_double_t crystal[], rmg_double_t cartesian[]);
void to_crystal (rmg_double_t crystal[], rmg_double_t cartesian[]);
rmg_double_t metric (rmg_double_t *crystal);

/* Md run types */
void relax (int steps, STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc,
              rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *rhoc);
void neb_relax (STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc,
              rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhocore, rmg_double_t *rhoc);
void cdfastrlx (STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc,
                rmg_double_t *rho, rmg_double_t *rhocore, rmg_double_t *rhoc);
void moldyn (STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc,
             rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhoc, rmg_double_t *rhocore);
void cholesky (rmg_double_t *a, int n);


/*the function for softpseudopotential*/
void aainit (int lli, int mix, int lx, int mx, int nlx, double ap[][9][9],
             int lpx[][9], int lpl[][9][9]);
rmg_double_t app_cil1 (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
               rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil1_bcc (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil1_fcc (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil1_hex (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
                   rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil1_ortho (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
                     rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void app_nls (rmg_double_t * psiR, rmg_double_t * psiI, rmg_double_t * workR, rmg_double_t * workI, rmg_double_t *work2R, rmg_double_t *work2I, 
     rmg_double_t *sintR, rmg_double_t *sintI, int state, int kidx);
void app_nls_allstates (rmg_double_t * psiR, rmg_double_t * psiI, rmg_double_t * workR, rmg_double_t * workI, rmg_double_t *work2R, rmg_double_t *work2I, 
     rmg_double_t *Bns, rmg_double_t *BnsI, rmg_double_t *sintR, rmg_double_t *sintI, int kidx);
void app_nls_batch (STATE *sp, rmg_double_t *nv, rmg_double_t *ns, rmg_double_t *Bns, rmg_double_t *sintR);
void get_ddd (rmg_double_t *veff);
void get_nlop_d (ION *iptr, rmg_double_t *rtptr, int ip, int icount, int *dvec);
void get_nlop_p (ION *iptr, rmg_double_t *rtptr, int ip, int icount, int *dvec);
void get_nlop_s (ION *iptr, rmg_double_t *rtptr, int ip, int icount, int *dvec);
void get_QI (void);
void get_qqq (void);
void get_rho (STATE * states, rmg_double_t * rho, rmg_double_t * rhocore);
void get_new_rho (STATE * states, rmg_double_t * rho);
void get_pdos (STATE * states, rmg_double_t Emin, rmg_double_t Emax, int E_POINTS);
void mix_rho (rmg_double_t * new_rho, rmg_double_t * rho, rmg_double_t *rhocore, int length, int length_x, int length_y, int length_z);
void init_psp (void);
void init_qfunct (void);
void mg_eig_state (STATE *sp, int tid, rmg_double_t *vtot_psi);
void mg_eig_state_f (STATE *sp, int tid, rmg_double_t *vtot_psi);
void ortho (STATE *states, int kpt);
rmg_double_t qval (int ih, int jh, rmg_double_t r, rmg_double_t invdr, rmg_double_t *ptpr, int *nhtol,
           int *nhtom, int *indv, rmg_double_t *ylm, rmg_double_t ap[][9][9], int lpx[][9],
           int lpl[][9][9], SPECIES *sp);
void reinit_ionic_pp (STATE * states, rmg_double_t * vnuc, rmg_double_t * rhocore, rmg_double_t * rhoc);

#if GAMMA_PT
void subdiag_gamma (STATE *states, rmg_double_t *vh, rmg_double_t *vnuc, rmg_double_t *vxc);
void subdiag_app_A (STATE * states, rmg_double_t * a_psi, rmg_double_t * s_psi, rmg_double_t * vtot_eig);
void subdiag_app_B (STATE * states, rmg_double_t * b_psi);
#else
void subdiag_nongamma (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
void subdiag_app_A (STATE * states, rmg_double_t * a_psiR, rmg_double_t * a_psiI, rmg_double_t * s_psiR, rmg_double_t * s_psiI, rmg_double_t * vtot_eig);
void subdiag_app_B (STATE * states, rmg_double_t * b_psiR, rmg_double_t * b_psiI);
#endif
void init_subdiag(void);


void ylmr2 (double *r, double *ylm);
rmg_double_t gcutoff (rmg_double_t g1, rmg_double_t gcut, rmg_double_t width);
void rft1 (rmg_double_t cparm, rmg_double_t *f, rmg_double_t *r, rmg_double_t *ffil, rmg_double_t *rab,
           int rg_points, int lval, rmg_double_t dr, rmg_double_t width, int lrg_points);
void norm_psi1 (STATE *sp, int istate, int kpt);
#if FAST_ORTHO
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t *cR, rmg_double_t *cI, rmg_double_t *Oij);
#else
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t *cR, rmg_double_t *cI);
#endif
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t cR, rmg_double_t cI);
rmg_double_t get_QnmL (int idx, int ltot, rmg_double_t r, SPECIES *sp);

void force (rmg_double_t *rho, rmg_double_t *rho_oppo, rmg_double_t *rhoc, rmg_double_t *vh, rmg_double_t *vxc, rmg_double_t *vnuc,
            STATE *states);


void iiforce (void);
void lforce (rmg_double_t *rho, rmg_double_t *vh);
void nlforce (rmg_double_t *veff);
void get_gamma (rmg_double_t * gammaR, int ion, int nh);
void partial_gamma (int ion, rmg_double_t * par_gammaR, rmg_double_t * par_omegaR, int nion, int nh,
                    rmg_double_t * newsintR_x, rmg_double_t * newsintR_y, rmg_double_t * newsintR_z,
                    rmg_double_t * newsintI_x, rmg_double_t * newsintI_y, rmg_double_t * newsintI_z);
void partial_betaxpsi (int ion, fftw_plan p2, rmg_double_t *newsintR_x,
                       rmg_double_t *newsintR_y, rmg_double_t *newsintR_z,
                       rmg_double_t *newsintI_x, rmg_double_t *newsintI_y,
                       rmg_double_t *newsintI_z, ION *iptr);
void partial_QI (int ion, rmg_double_t *QI_R, ION *iptr);
void qval_R (int ih, int jh, rmg_double_t r, rmg_double_t *x, rmg_double_t *qlig, rmg_double_t *drqlig,
             rmg_double_t invdr, int *nhtol, int *nhtom, int *indv, rmg_double_t *ylm,
             rmg_double_t *ylm_x, rmg_double_t *ylm_y, rmg_double_t *ylm_z, rmg_double_t ap[][9][9],
             int lpx[][9], int lpl[][9][9], rmg_double_t *Q_x, rmg_double_t *Q_y,
             rmg_double_t *Q_z, SPECIES *sp);
void ylmr2_x (double *r, double *ylm_x);
void ylmr2_y (double *r, double *ylm_y);
void ylmr2_z (double *r, double *ylm_z);
void nlccforce (rmg_double_t *rho, rmg_double_t *vxc);
rmg_double_t get_ve_nl (STATE *sta, int istate);
void pack_rho_ctof (rmg_double_t *rhoc, rmg_double_t *rhof);
void bspline_interp_full (rmg_double_t *rho, rmg_double_t *rho_f);
void get_vtot_psi (rmg_double_t * vtot_psi, rmg_double_t * vtot, int grid_ratio);
void betaxpsi (STATE *states);
void betaxpsi1 (STATE *states, int kpt);
void assign_weight (SPECIES *sp, int ion, fftw_complex *beptr,
                    rmg_double_t *rtptr, rmg_double_t *Bweight);
void assign_weight2 (int nldim, int ion, rmg_double_t *beptr, rmg_double_t *rtptr);
void pack_gftoc (SPECIES *sp, fftw_complex *gwptr, fftw_complex *gbptr);
void debug_write_rho_z (rmg_double_t *rhoz);
void print_density_z_direction (int grid_x, int grid_y, rmg_double_t *density,
                                int px0_grid, int py0_grid, int pz0_grid,
                                rmg_double_t zside);
void get_derweight (int ion, rmg_double_t *beta_x, rmg_double_t *beta_y, rmg_double_t *beta_z,
                    ION *iptr, fftw_plan p2);
void partial_beta_fdiff (fftw_complex *beptr, int nldim, rmg_double_t *beta_x,
                         rmg_double_t *beta_y, rmg_double_t *beta_z);

void mulliken (STATE *states);
rmg_double_t ylm(int l, rmg_double_t *r);
int listlen (FILE * fh, char *id);
int del_space (FILE * fh, int *tchr, int isdata);
void print_matrix(double *b, int n, int ldb);
void sl_init(int *ictxt, int size);
void sl_exit(int ictxt);
void set_desca(int *desca, int *ictxt, int size);
void distribute_mat(int *desca, double *bigmat, double *dismat, int *size);
void matinit(int *desca, double *dismat, double *globmat, int size);
int matsum_packbuffer(int row, int col, double *buffer, double *globmat, int size);
void reduce_and_dist_matrix(int n, rmg_double_t *global_matrix, rmg_double_t *dist_matrix, rmg_double_t *work);
void print_distribute_mat(double *dismat, int *desca, int size);
void init_efield (rmg_double_t * vnuc);
void pulay_rho(int step, int N, int N_x, int N_y, int N_z, double *rho_new, double *rho_old, int NsavedSteps, rmg_double_t ***hist, rmg_double_t ***rhist, int special_metric, rmg_double_t weight);
void mg_restrict_4th (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz);
void mix_betaxpsi (int mix);
void mix_betaxpsi1 (STATE *sp);
int claim_ion (rmg_double_t *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid);
int find_grid_1d_owner(int igrid, int tgrid, int pgrid);
void find_grid_owner(int igridx, int igridy, int igridz, int nxgrid, int nygrid, int nzgrid, int *xnode, int *ynode, int *znode);
int is_loop_over_states(void);
int rmg_is_open_mp_safe(void);
void RMG_set_omp_parallel(void);
void RMG_set_omp_single(void);
void Papi_init_omp_threads(int ithread);
void Papi_finalize_omp_threads(int ithread);
void fastrelax (rmg_double_t *dt, rmg_double_t dt_max, rmg_double_t dt_inc, rmg_double_t dt_dec, int n_min, int *n_count);
void fire (rmg_double_t *step, rmg_double_t step_max, rmg_double_t f_inc, rmg_double_t f_dec, int n_min, int *n_count );
void quick_min (void);
int int_sum_all (int x, MPI_Comm comm);
void move_ions (rmg_double_t dt);
void get_extrapolation_constants (rmg_double_t *alpha, rmg_double_t *beta);
void lcao_init (void);
void init_atomic_rho_wf (void);
void lcao_get_rho (rmg_double_t * arho_f);
void lcao_get_awave (rmg_double_t *psi, ION *iptr, int awave_idx, int l, int m, double coeff);
void lcao_get_psi (STATE * states);
rmg_double_t mask_function (rmg_double_t x);
void apply_mask_function (rmg_double_t *f, rmg_double_t * r, int rg_points, rmg_double_t rmax, rmg_double_t offset);
void filter_potential (rmg_double_t *potential, rmg_double_t *r, int rg_points, rmg_double_t rmax, rmg_double_t offset, rmg_double_t parm, rmg_double_t* potential_lgrid, 
	rmg_double_t *rab, int l_value, rmg_double_t dr, rmg_double_t  gwidth, int lgrid_points, rmg_double_t rcut, rmg_double_t rwidth, rmg_double_t * drpotential_lgrid);
int test_overlap (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid);
void RMG_MPI_trade(rmg_double_t *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
void RMG_MPI_trade_f(rmg_float_t *buf, int count, int type, int pe_x_offset, int pe_y_offset, int pe_z_offset, MPI_Comm comm, int tag, MPI_Request *req);
void init_trade_imagesx_async(void);
void  get_rho_oppo (rmg_double_t * rho, rmg_double_t * rho_oppo);
void get_opposite_eigvals (STATE * states);
void get_opposite_occupancies (STATE * states);
#endif
#if __cplusplus
}
#endif

