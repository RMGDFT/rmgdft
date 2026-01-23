#pragma once
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

#ifndef RMG_COMMON_PROTOTYPES1_H
#define RMG_COMMON_PROTOTYPES1_H 1


#include <stdbool.h>

/* Function prototypes */
void get_dipole (double * rho, double *dipole);
void get_dipole (double * rho, double *, double *dipole);
void app6_del2 (double *rho, double *work, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz);
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
double Fill_on (std::vector<double> &eig_all, std::vector<double> &kweight, std::vector<double> &occ, double width, double nel, double mix,
           int occ_flag, int mp_order);

void find_phase (int nlxdim, int nlydim, int nlzdim, double * nlcdrs, double ** phase_sin, double ** phase_cos);
void get_phase (ION *iptr, double *rtptr, int icount, int *dvec);
char *get_num (char *str);


void get_xc (double * nrho, double * nrho_oppo,  double * vxc, double * exc, int xctype);

void get_zdens (STATE *states, int state, double *zvec);


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
void symmetrize_tensor (double *mat);
void symforce (void);
void rmg_timings (int what, double time);
double minimage (ION *ip1, ION *ip2, double *xtal_r);
double my_crtc (void);
FILE *open_restart_file (char *filename);
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


void sortpsi (STATE *states);
void set_bc (double *mat, int dimx, int dimy, int dimz, int images, double val);
void set_bcx (double *mat, int dimx, int dimy, int dimz, int images, double val);
void vol_wf (STATE *states, int state, int step);
void write_avgd (double *rho);
void write_avgv (double *vh, double *vnuc);

void write_header (void);
void write_force (void);
void write_timings (void);
double rand0 (long *idum);
void cgen_prolong(double coef[], double fraction, int order);
void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz, int grid_ratio, int order);
void get_milliken (STATE *states);


int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz,
               int *ilow, int *ihi, int *jlow, int *jhi, int *klow,
               int *khi, int cdim, int pxgrid, int pygrid, int pzgrid,
               int nxgrid, int nygrid, int nzgrid, double * xcstart, double * ycstart, double * zcstart);

double linint (double *y, double rv, double invdr);


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
void cdfastrlx (STATE *states, double *vxc, double *vh, double *vnuc,
                double *rho, double *rhocore, double *rhoc);
void moldyn (STATE *states, double *vxc, double *vh, double *vnuc,
             double *rho, double *rho_oppo, double *rhoc, double *rhocore);
void cholesky (double *a, int n);


/*the function for softpseudopotential*/
void get_ddd (double *veff);
void get_QI (void);
void get_qqq (void);
void get_qqq_dk (double dk[3], std::complex<double> *qqq_dk, std::complex<double> *qqq_dk_so);
void get_rho (STATE * states, double * rho, double * rhocore);
void get_new_rho (STATE * states, double * rho);
void mix_rho (double * new_rho, double * rho, double *rhocore, int length, int length_x, int length_y, int length_z);
void mg_eig_state (STATE *sp, int tid, double *vtot_psi);
//void ortho (STATE *states, int kpt);
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
void bspline_interp_full (double *rho, double *rho_f);
void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
void get_vxc_exc (double * nrho, double * nrho_oppo,  double * vxc, double * exc, int xctype);
void betaxpsi (STATE *states);
void pack_gftoc (SPECIES *sp, fftw_complex *gwptr, fftw_complex *gbptr);
void print_density_z_direction (int grid_x, int grid_y, double *density,
                                int px0_grid, int py0_grid, int pz0_grid,
                                double zside);

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
void init_efield (double * vnuc, double efeild[3]);
int claim_ion (double *xtal,  int pxgrid, int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid);
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
void inline progress_tag(void)
{
    if(pct.gridpe==0)
        fprintf(ct.logfile, "[ %3d %3d %4d %8.0f ] %s: ",
           ct.md_steps, ct.scf_steps, ct.total_scf_steps,
           my_crtc() - ct.time0, (strrchr(__FILE__, '/')+1));
}

#endif
