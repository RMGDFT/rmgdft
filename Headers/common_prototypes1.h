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
void get_dipole (double * rho);
void app6_del2 (double *rho, double *work, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz);
void app_smooth (double *f, double *work, int dimx, int dimy, int dimz);
void app_smooth_f (float * f, float * work, int dimx, int dimy, int dimz);
void app_smooth1 (double *f, double *work, int dimx, int dimy, int dimz);
void app_smooth1_f (float *f, float *work, int dimx, int dimy, int dimz);
void app_cir_driver (double *a, double *b, int dimx, int dimy, int dimz, int order);
void app_cir_fourth (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_sixth (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir (double *a, double *b, int dimx, int dimy, int dimz);
double app_cil (double *a, double *b, int dimx, int dimy, int dimz, double gridhx,
              double gridhy, double gridhz);
double app_cil_driver (double * a, double * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);
double app_cil_fourth (double *a, double *b, int dimx, int dimy, int dimz, double gridhx,
              double gridhy, double gridhz);
double app_cil_sixth (double *psi, double *b, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz);
void app_grad (double  * rho, double *wx, double *wy, double *wz, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);
void constrain( void );
void cross_product (double *a, double *b, double *c);
double Fill_on (STATE *states, double width, double nel, double mix,
           int num_st, int occ_flag, int mp_order);

void find_phase (int nlxdim, int nlydim, int nlzdim, double * nlcdrs, double ** phase_sin, double ** phase_cos);
void get_phase (ION *iptr, double *rtptr, int icount, int *dvec);
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

void mgga_libxc (double * rho, double * tau, double * vxc, double * exc, int mode);
void get_mgga_exc_vxc (double * rho, double * rho_oppo, double * tau, double * vxc, double * exc);


double get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target, int boundaryflag);
void global_sums (double *vect, int *length, MPI_Comm comm);
void init_pe ( int image );
STATE *init_states (void);
void init_weight (void);
void weight_shift_center(SPECIES * sp, fftw_complex * weptr);
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
void output_eigenvalues( STATE *states, int ikbs, int iscf );
void pack_ptos (double *sg, double *pg, int dimx, int dimy, int dimz);
void pack_stop (double *sg, double *pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (double *sg, double *pg, double alpha, int dimx, int dimy, int dimz);
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

double real_sum_all (double x, MPI_Comm comm);
double double_sum_all (double x, MPI_Comm comm);
void FftFilterFine(double *x,  double factor);


void sortpsi (STATE *states);
void set_bc (double *mat, int dimx, int dimy, int dimz, int images, double val);
void set_bcx (double *mat, int dimx, int dimy, int dimz, int images, double val);
void getpoi_bc (double *rho, double *vh_bc, int dimx, int dimy, int dimz);
void vol_wf (STATE *states, int state, int step);
void write_avgd (double *rho);
void write_avgv (double *vh, double *vnuc);
void write_zstates (STATE *states);

void write_header (void);
void write_force (void);
void write_timings (void);
void wvfn_residual(STATE *states);
double rand0 (long *idum);
void cgen_prolong(double coef[], double fraction, int order);
void mg_restrict_6 (double * full, double * half, int dimx, int dimy, int dimz, int grid_ratio);
void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz, int grid_ratio, int order);
void get_milliken (STATE *states);


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
void get_ddd (double *veff);
void get_QI (void);
void get_qqq (void);
void get_rho (STATE * states, double * rho, double * rhocore);
void get_new_rho (STATE * states, double * rho);
void get_pdos (STATE * states, double Emin, double Emax, int E_POINTS);
void mix_rho (double * new_rho, double * rho, double *rhocore, int length, int length_x, int length_y, int length_z);
void Output_rho_xsf(double *array_3d, MPI_Comm comm);
void init_qfunct (void);
void mg_eig_state (STATE *sp, int tid, double *vtot_psi);
void ortho (STATE *states, int kpt);
void init_pestr(void);
void read_init(char *meta);
int filename_increment(char *filename);
void init_subdiag(void);


void ylmr2 (double *r, double *ylm);
double gcutoff (double g1, double gcut, double width);
void rft1 (double cparm, double *f, double *r, double *ffil, double *rab,
           int rg_points, int lval, double dr, double width, int lrg_points);
void norm_psi1 (STATE *sp, int istate, int kpt);
double get_QnmL (int idx, int ltot, double r, SPECIES *sp);

void force (double *rho, double *rho_oppo, double *rhoc, double *vh, double *vxc, double *vnuc);


void iiforce (double *force);
void lforce (double *rho, double *vh, double *force);
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
void nlccforce (double *rho, double *vxc, double *force);
void pack_rho_ctof (double *rhoc, double *rhof);
void bspline_interp_full (double *rho, double *rho_f);
void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
void get_vxc_exc (double * nrho, double * nrho_oppo,  double * vxc, double * exc, int xctype);
void betaxpsi (STATE *states);
void pack_gftoc (SPECIES *sp, fftw_complex *gwptr, fftw_complex *gbptr);
void debug_write_rho_z (double *rhoz);
void print_density_z_direction (int grid_x, int grid_y, double *density,
                                int px0_grid, int py0_grid, int pz0_grid,
                                double zside);

void mulliken (STATE *states);
double ylm(int l, double *r);
int listlen (FILE * fh, char *id);
void print_matrix(double *b, int n, int ldb);
void sl_init(int *ictxt, int size);
void sl_exit(int ictxt, int NotDone);
void set_desca(int *desca, int *ictxt, int size);
void distribute_mat(int *desca, double *bigmat, double *dismat, int *size);
void matinit(int *desca, double *dismat, double *globmat, int size);
int matsum_packbuffer(int row, int col, double *buffer, double *globmat, int size);
void reduce_and_dist_matrix(int n, double *global_matrix, double *dist_matrix, double *work);
void init_efield (double * vnuc);
void pulay_rho(int step, int N, int N_x, int N_y, int N_z, double *rho_new, double *rho_old, int NsavedSteps, double ***hist, double ***rhist, int special_metric, double weight);
void mix_betaxpsi (int mix);
int claim_ion (double *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid);
int is_loop_over_states(void);
void fastrelax (double *dt, double dt_max, double dt_inc, double dt_dec, int n_min, int *n_count);
void fire (double *step, double step_max, double f_inc, double f_dec, int n_min, int *n_count );
int int_sum_all (int x, MPI_Comm comm);
void move_ions (double dt);
void lcao_init (void);
void filter_potential (double *potential, double *r, int rg_points, double rmax, double offset, double parm, double* potential_lgrid, 
	double *rab, int l_value, double dr, double  gwidth, int lgrid_points, double rcut, double rwidth, double * drpotential_lgrid);
int test_overlap (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid);
void init_trade_imagesx_async(void);
void  get_rho_oppo (double * rho, double * rho_oppo);
void get_opposite_eigvals (STATE * states);
void get_opposite_occupancies (STATE * states);
void get_tf_rho (double * tf_rho);
void get_dipole (double * rho);
void set_pbc(double *position, int num_ions, int num_images);

#if __cplusplus
}
#endif

