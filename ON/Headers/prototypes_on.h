#define GAMMA_PT 1
//#include "typedefs.h"
void KbpsiComm();
void InitNonlocalComm();
void GetHS(STATE * states, STATE * states1, double *vtot_c, double *Aij, double *Bij);
void GetHvnlij (double *Aij, double *Bij);
void OrbitalOptimize (STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, double * vxc_old,
        double * vh_old);
void GetNewRho_on (STATE * states, double *rho, double *rho_matrix);
void RhoQnmMat (double *Aij, double * global_mat_X);
void RhoAugmented (double * rho, double * global_mat_X);
void UpdatePot(double *vxc, double *vh, double * vxc_old, double * vh_old,
        double *vnuc, double *rho, double *rho_oppo, double *rhoc, double *rhocore);
void KbpsiUpdate(STATE *states);
void FillOrbitalBorders(double *, double *, int, int, int, int);
void VhDriver(double *rho, double *rhoc, double *vh);
void CheckConvergence(double *vxc, double *vh, double * vxc_old, double * vh_old, double *, double *rho_old, int *CONVERGENCE);
double gaussintegral(double, int);
void DipoleCorrection(double *rho, double *rhoc, double *vcorr, double *vhx, double *vhy, double *vhz);
void VhcorrPeriodicPart(double *vh_x, double *vh_y, double *vh_z, double alpha, double *r0);
void VhcorrDipoleInit(double *vh_x, double *vh_y, double *vh_z, double *rhoc);
void GetNlop_on(void);


#ifdef __cplusplus

#include <unordered_map>
#include "InputKey.h"
#include "FiniteDiff.h"

extern "C" void precond_mg_c(double *res, double *work1, double *work2, int istate);
void InitON(double * vh, double * rho, double *rho_oppo,  double * rhocore, double * rhoc,
          STATE * states, STATE * states1, double * vnuc, double * vxc, double * vh_old, 
          double * vxc_old, std::unordered_map<std::string, InputKey *>& ControlMap);
void PrecondMg(double *psiR, double *work1, STATE *sp);
void Precond(double *x);
void ZeroBoundary(double *a, int ixx, int iyy, int izz);
void ZeroBoundary(double *a, int ixx, int iyy, int izz, int width);
void MgridSolvLocal(double * v_mat, double * f_mat, double * work,
                int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                int *post_cyc, double step,
                int mu_cyc, int istate, int *iion, double Zfac);
void EvalResidualLocal (FiniteDiff *FD, double *mat, double *f_mat, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, double *res);
void SolvPoisLocal (FiniteDiff *FD, double *vmat, double *fmat, double *work,
                int dimx, int dimy, int dimz, double gridhx,
                double gridhy, double gridhz, double step, double Zfac, double k);


extern "C" {
#endif


void write_rho_x(double *, char*);
void write_rho_y(double *, char*);
void write_rho_z(double *, char*);
void orbital_comm(STATE *);
void init_wf_atom(STATE *);
void init_wf_lcao(STATE *);

void nlforce_par_omega(double *, double *, double *, int, int);
void nlforce_par_rho(double *, double *, double *, int, int);
void nlforce_par_Q(double *, double *, int, int, double *);
int prime_factors(int, int *);
void get_vxc(double*, double *, double *, double*);
void assign_weight_on(SPECIES *, fftw_complex *, double *);
void quench(STATE *, STATE *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
void app10_del2(double *, double *, int, int, int, double, double, double);
void sl_init_comm(int *, int, int, MPI_Comm comm);
void read_data(char *, double *, double *, double *, double *, double *, double *, STATE *);
double get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target, int boundaryflag);
double app_cil_orbital (double * a, double * b, int dimx, int dimy, int dimz,
              double gridhx, double gridhy, double gridhz);
double app_cil_orbital6 (double * a, double * b, int dimx, int dimy, int dimz,
              double gridhx, double gridhy, double gridhz);
void app_cir_fourth (double *a, double *b, int dimx, int dimy, int dimz);
void app_cir_sixth (double *a, double *b, int dimx, int dimy, int dimz);
double app_cil_fourth (double *a, double *b, int dimx, int dimy, int dimz, double gridhx,
              double gridhy, double gridhz);
double app_cil_sixth (double *psi, double *b, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz);
void trade_imagesx (double *f, double *w, int dimx, int dimy, int dimz, int images, int type);
void constrain( void );
void init_pe ( int image );
void mg_restrict_6 (double * full, double * half, int dimx, int dimy, int dimz, int grid_ratio);

void fill_orbit_borders4(double * sg, double * pg, int dimx, int dimy, int dimz);
int open_wave_file (char *filename);
void get_vtot_psi (double * vtot_psi, double * vtot, int grid_ratio);
void init_nonlocal_comm();
void normalize_orbits( STATE *states);
void init_dimension(int *, int* );
void precond(double *f);
int int_max_all(int);

void pulay_weighted (int step0, int N, double *xm, double *fm, int NsavedSteps,
                int Nrefresh, double scale, int preconditioning);
void get_invmat(double *);
void distribute_to_global(double *vtot_c, double *vtot_global);

void init_IO ( int argc, char **argv );
void get_dm_diag_p(STATE *states, double *l_s, double *mat_X, double *Hij);

void get_HS(STATE * states, STATE * states1, double *vtot_c, double *Aij, double *Bij);
void interpolate_atom_density(double *rho_tem, double *rho_out, int ixmin, int ixmax,
                              int iymin, int iymax, int izmin, int izmax,
                              double *crds, double *crds1, double hxgrid,
                              double hygrid, double hzgrid,
                              double hxgrid_new, double hygrid_new, double hzgrid_new);


void sl_init_on(int *ictxt, int nprow, int npcol);

void init_rho_atom(double *rho);

void init_wf_gaussian(STATE * states);


void global_to_distribute(double * global_array, double * distr_array);

void diag_eig_matrix(double *eig_matrix, double *eig, int *desca);

double mask_function (double x);
void apply_mask_function (double *f, double * r, int rg_points, double rmax, double offset);
void backout_mask_function (double *f, double dr, int lg_points, double rmax);
void thread_barrier_wait(void);
int is_loop_over_states(void);
void cgen_prolong(double coef[], double fraction, int order);
void xbsmovie (FILE *movie);

void app_grad6(double * f, double * wx, double * wy, double * wz, int dimx, int dimy, int dimz);


STATE *init_states();

void mgrid_solv_local(double * v_mat, double * f_mat, double * work,
                int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                int *post_cyc, double step, int mu_cyc, 
                int istate, int *iion, int flag_local, double Zfac);

/* Function prototypes */
void app_4del2 (double * f, double * work);

void app_nl (double * psiR, double * psiI, double * workR, double * workI, int state,
        int flag, int kidx, int tid);
void cholesky (double * a, int n);
void cross_product (double * a, double * b, double * c);
double fill (STATE * states, double width, double nel, double mix,
        int num_st, int occ_flag);
void genvpsi (double * psi, double * twovpsi, double * pvtot, double * pvnl,
        double * kd, double kmag, int dimx, int dimy, int dimz);
void get_index_loc (STATE *);
void get_eig (STATE * states, double * vxc, double * vh, double * vnuc);
char *get_num (char *str);
void get_zdens (STATE * states, int state, double * zvec);
void xcgga (double * rho, double * vxc, double * exc, int flag);
void gram (KPOINT * kpoint, double h, int numst, int maxst, int numpt,
        int maxpt);
double get_ke (STATE * sp, int tid);
void global_sums (double * vect, int *length, MPI_Comm comm);
void init_pe_on (void);
void init_pegrid (void);
void init_wf (STATE * states);
void init_wflcao (STATE * states);
void init_nuc (double * vnuc, double * rhoc, double * rhocore);
void init_pos ();
void init_sym (void);
void symmetrize_rho (double * rho);
void symforce (void);
void rmg_timings (int what, double time);
void mg_eig_state (STATE * states, int tid, double * vtot);
double minimage (ION * ip1, ION * ip2, double * xtal_r);
double my_crtc (void);
/* Local function prototypes */
void norm_psi (STATE * psi, double bv);
void ortho_full (STATE * states);
void ortho_half (STATE * states);
void ortho_bcc (STATE * states);
void output (STATE * states, int its);
void pe2xyz (int pe, int *x, int *y, int *z);
void pack_ptos (double * sg, double * pg, int dimx, int dimy, int dimz);
void pack_stop (double * sg, double * pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (double * sg, double * pg, double alpha, int dimx, int dimy,
        int dimz);
void pack_ptos_trade (double * sg, double * pg, int dimx, int dimy, int dimz);
void pack_vhstod (double * s, double * d, int dimx, int dimy, int dimz, int boundarflag);
void pack_vhdtos (double * s, double * d, int dimx, int dimy, int dimz, int boundarflag);
double radint (double *f, double *r, int n, double al);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
double real_sum_all (double x, MPI_Comm comm);
double real_max_all (double x);
void sortpsi (STATE * states);
void subdiag (STATE * states, double * vh, double * vnuc, double * vxc);
void trade_images (double * mat, int dimx, int dimy, int dimz, int type);
void trade_images_mpi (double * mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_images_smp (double * mat, int dimx, int dimy, int dimz, int *nb_ids);
void set_bc (double * mat, int dimx, int dimy, int dimz, int images, double val);
void getpoi_bc (double * rho, double * vh_bc, int dimx, int dimy, int dimz);
void vol_rho (double * rho, int step);
void vol_wf (STATE * states, int state, int step);
void write_avgd (double * rho);
void write_avgv (double * vh, double * vnuc);
void write_zstates (STATE * states);
void write_header (void);
void write_pos (void);
void write_eigs (STATE * states);
void write_timings (void);
double rand0 (long *idum);

void gather_psi (double * tmp_psiR, double * tmp_psiI, STATE * sp, int tid);
void scatter_psi (double * tmp_psiR, double * tmp_psiI, STATE * sp, int tid);
void get_milliken (STATE * states);

void output_wave (STATE * states, int kpt, int fhand);



/* Conversion between crystal and cartesian coordinate prototypes */
void recips (void);
void to_cartesian (double crystal[], double cartesian[]);
void to_crystal (double crystal[], double cartesian[]);
double metric (double * crystal);

/* Md run types */
void run (STATE * states, STATE * states1);
void fastrlx( STATE *states, STATE *states1, double *vxc, double *vh, double *vnuc, double *vh_old, double *vxc_old,
        double *rho, double *rhocore, double *rhoc );
void cdfastrlx (STATE * states, double * vxc, double * vh, double * vnuc,
        double * rho, double * rhocore, double * rhoc);
/*void moldyn (STATE * states, double * vxc, double * vh, double * vnuc, double * rho, double * rhoc, double * rhocore);*/
void dx (STATE * states, double * vxc, double * vh, double * vnuc, double * rho,
        double * rhoc);
void psidx (STATE * states, double * vxc, double * vh, double * vnuc, double * rho,
        double * rhoc);
void cholesky (double * a, int n);


double minimage1 (double aa[3], double bb[3]);
void init_parameter (STATE * states);
void make_mask_grid (double rcut, int level, STATE * states);
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
void xclsd_pz81 (double * rho, double * vxc);
void global_sums_int (int *, int *);
double linint (double * y, double r, double invdr);
void ortho_norm_local (STATE * states);
void global_sums_int (int *vect, int *length);
void my_barrier (void);
void matrix_and_diag (STATE * states, STATE * states1, double * vxc, int flag);
void precond_mg (double *res, double *work1, double *work2, int istate);

void matS_cholesky_real (STATE * states);
void get_invBmat (STATE * states, double *matB);
void print_matrix (double *b, int n, int ldb);
void get_overlap_real (double *aa, int numst, int numpt,
        int lda, double *ss, int lds);
void get_cholesky_real (double *matS);
void get_all_kbpsi (STATE * states1, STATE * states, ION_ORBIT_OVERLAP *, double *, double *);
void get_Hvnlij (double *Aij, double *Bij);
void genvlocpsi (double * psi, int st1, double * work1, double * vtot_global, STATE * states);
void genvnlpsi (double *sg_twovpsi, double *vnl,
        int dimx, int dimy, int dimz);
void get_wave(int st, STATE * states);
void add_orbit_to_wave(int st1, double scale, double * psi1, double * wave_global, STATE * states);
void print_wave(int wave_plot, STATE * states, int coarse_level);
void fsymforces (double * force, int *s, int *irg, int *irt, int *nat,
        int *ibrav, int *nsym, double * celldm, int *nr1, int *nr2,
        int *nr3);
void symrho (double * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
        int *irg, int *ftau);

void sgama (int *nrot,int *nat,double *s,double *at,double *bg,double *tau,
        int *ityp,int *nsym,int *nr1,int *nr2,int *nr3,int *irg,
        int *irt,double *ftau,double *rtau,int *npk,int *nks,double *xk,
        double *wk,double *xau,double *rau,bool *invsym, int *wflag);

void line_min_three_point (STATE *, STATE *, double, double, double *, double *,
        double *, double *, double *);

void dot_product_orbit_orbit (STATE *orbit1, STATE *orbit2, STATE
        *orbit3, double *H, double *S);

void orbit_dot_orbit (STATE * states, STATE * states1, double *Hij_row, double * Bij_row);

void app_mask (int istate, double *u, int level);
void allocate_masks (STATE * states);

void state_corner_xyz (STATE * states);

char *get_symbol (int atomic_number);

void get_matB_qnm (double *Aij);
void pack_vtot_ftoc (double * vtot, double * vtot_c);
void qnm_beta_betapsi (STATE *sp, int ion2, double * prjptr);
void pack_rho_ctof (double * rho1, double * rho_f);
void rho_augmented (double * rho, double * global_mat_X,
int *state_begin, int *state_end, int *num_nonlocal_ion, double *kbpsi,
int max_ion_nonlocal, double *kbpsi_comm, int *ionidx_allproc);

void rho_Qnm_mat (double *Aij, double * global_mat_X,
int *state_begin, int *state_end, int *num_nonlocal_ion, double *kbpsi,
int max_ion_nonlocal, double *kbpsi_comm, int *ionidx_allproc);
void rho_nm_mat (double *Aij, double * global_mat_X);
int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz, int *ilow, int *ihi,
        int *jlow, int *jhi, int *klow, int *khi, int cdim, int pxgrid,
        int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid,
        double * xcstart, double * ycstart, double * zcstart);
void xcgga (double * rho, double * vxc, double * exc, int mode);
double radint1 (double * func, double * r, double * rab, int n);
void get_all_partial_kbpsi (STATE *states, ION_ORBIT_OVERLAP *,
        double *, double *, double *, double *,
        double *, double *);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, double * global_mat_X);
void partial_QI (int ion, double * QI_R, ION *iptr);
void partial_nlop_s (ION * iptr, double * betax, double * betay, double * betaz,
        int ip, fftw_plan p1, fftw_plan p2);
void partial_nlop_p (ION * iptr, double * betax, double * betay, double * betaz,
        int ip, fftw_plan p1, fftw_plan p2);
void partial_nlop_d (ION * iptr, double * betax, double * betay, double * betaz,
        int ip, fftw_plan p1, fftw_plan p2);
void get_mat_Omega (STATE * states, double Omega[]);
void md_fastrelax(void);
void change_states_crds (STATE * states);

/*end shuchun wang */

FILE *open_xbs_movie (char *filename);






double get_te_ion_ion ();
double get_sum_eig (STATE * states);
double get_Exc (double * rho, double * rhocore);
void correct_res (STATE *, double *);
void get_state_to_proc (STATE * states);



/* different methods to update orbitals */
void kain (int step, int N, double *xm, double *fm, int NsavedSteps);
void pulay (int step, int N, double *xm, double *fm, int NsavedSteps,
        int preconditioning);
void sd (int step, int N, double *xm, double *fm);
void pulay_mixing (int size, double *rho, int NsavedSteps);
void charge_pulay (int step, int N, double *xm, double *fm, int NsavedSteps);
void pulay_kain (int step, int N, double *xm, double *fm, int NsavedSteps);

/*  Moving center stuff  */
void get_orbit_center (STATE *state, double * x, double * y, double * z);
void update_orbit_centers (STATE * states);
int if_update_centers (STATE * states);

void write_states_info (char *outfile, STATE * states);
void read_states_info (char *outfile, STATE * states);

void init_state_size (STATE * states);
void diag_eig_maitrix(double *, double *, int *);



void filter_potential (double *potential, double *r, int rg_points, double rmax, double offset, double parm, double* potential_lgrid, 
        double *rab, int l_value, double dr, double  gwidth, int lgrid_points, double rcut, double rwidth, double * drpotential_lgrid);


/* Function prototypes */
void app_4del2 (double * f, double * work);

void app_nl (double * psiR, double * psiI, double * workR, double * workI, int state,
        int flag, int kidx, int tid);
void cholesky (double * a, int n);
void cross_product (double * a, double * b, double * c);
double fill (STATE * states, double width, double nel, double mix,
        int num_st, int occ_flag);
void genvpsi (double * psi, double * twovpsi, double * pvtot, double * pvnl,
        double * kd, double kmag, int dimx, int dimy, int dimz);
void get_index_loc (STATE *);
void get_eig (STATE * states, double * vxc, double * vh, double * vnuc);
char *get_num (char *str);
void get_zdens (STATE * states, int state, double * zvec);
void xcgga (double * rho, double * vxc, double * exc, int flag);
void gram (KPOINT * kpoint, double h, int numst, int maxst, int numpt,
        int maxpt);
double get_ke (STATE * sp, int tid);
void global_sums (double * vect, int *length, MPI_Comm comm);
void init_pe_on (void);
void init_pegrid (void);
void init_wf (STATE * states);
void init_wflcao (STATE * states);
void init_nuc (double * vnuc, double * rhoc, double * rhocore);
void init_pos ();
void init_sym (void);
void symmetrize_rho (double * rho);
void symforce (void);
void rmg_timings (int what, double time);
void mg_eig_state (STATE * states, int tid, double * vtot);
double minimage (ION * ip1, ION * ip2, double * xtal_r);
double my_crtc (void);
/* Local function prototypes */
void norm_psi (STATE * psi, double bv);
void ortho_full (STATE * states);
void ortho_half (STATE * states);
void ortho_bcc (STATE * states);
void output (STATE * states, int its);
void pe2xyz (int pe, int *x, int *y, int *z);
void pack_ptos (double * sg, double * pg, int dimx, int dimy, int dimz);
void pack_stop (double * sg, double * pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (double * sg, double * pg, double alpha, int dimx, int dimy,
        int dimz);
void pack_ptos_trade (double * sg, double * pg, int dimx, int dimy, int dimz);
void pack_vhstod (double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag);
void pack_vhdtos (double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag);
double radint (double *f, double *r, int n, double al);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
double real_sum_all (double x, MPI_Comm comm);
double real_max_all (double x);
void sortpsi (STATE * states);
void subdiag (STATE * states, double * vh, double * vnuc, double * vxc);
void trade_images (double * mat, int dimx, int dimy, int dimz, int type);
void trade_images_mpi (double * mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_images_smp (double * mat, int dimx, int dimy, int dimz, int *nb_ids);
void set_bc (double * mat, int dimx, int dimy, int dimz, int images, double val);
void getpoi_bc (double * rho, double * vh_bc, int dimx, int dimy, int dimz);
void vol_rho (double * rho, int step);
void vol_wf (STATE * states, int state, int step);
void write_avgd (double * rho);
void write_avgv (double * vh, double * vnuc);
void write_zstates (STATE * states);
void write_header (void);
void write_pos (void);
void write_eigs (STATE * states);
void write_occ (STATE * states);
void write_timings (void);
double rand0 (long *idum);

void gather_psi (double * tmp_psiR, double * tmp_psiI, STATE * sp, int tid);
void scatter_psi (double * tmp_psiR, double * tmp_psiI, STATE * sp, int tid);
void get_milliken (STATE * states);

void output_wave (STATE * states, int kpt, int fhand);



/* Conversion between crystal and cartesian coordinate prototypes */
void recips (void);
void to_cartesian (double crystal[], double cartesian[]);
void to_crystal (double crystal[], double cartesian[]);
double metric (double * crystal);

/* Md run types */
void run (STATE * states, STATE * states1);
void fastrlx( STATE *states, STATE *states1, double *vxc, double *vh, double *vnuc, double *vh_old, double *vxc_old,
        double *rho, double *rhocore, double *rhoc );
void cdfastrlx (STATE * states, double * vxc, double * vh, double * vnuc,
        double * rho, double * rhocore, double * rhoc);
/*void moldyn (STATE * states, double * vxc, double * vh, double * vnuc, double * rho, double * rhoc, double * rhocore);*/
void dx (STATE * states, double * vxc, double * vh, double * vnuc, double * rho,
        double * rhoc);
void psidx (STATE * states, double * vxc, double * vh, double * vnuc, double * rho,
        double * rhoc);
void cholesky (double * a, int n);


double minimage1 (double aa[3], double bb[3]);
void init_parameter (STATE * states);
void get_mehr (void);
void make_mask_grid (double rcut, int level, STATE * states);
void make_mask_grid_state (int level, STATE * states);
void allocate_func (STATE * states, int inist);
void allocate_matrix ();
void xyz2pe (int x, int y, int z, int *pe);
void get_normKB (SPECIES * sp, double *pd);
void get_start (int cdim, double crds, double cstart, double hgrid,
        int *istart, double *nlcstart);
void get_index_array (int *pvec, int *dvec, int dim,
        int *Aix, int *Aiy, int *Aiz, int *icount,
        int index_low[3], int index_high[3]);
void get_Ai (int cdim, int *Ai, int ngrid, int istart);
void xclsd_pz81 (double * rho, double * vxc);
void global_sums_int (int *, int *);
double linint (double * y, double r, double invdr);
void ortho_norm_local (STATE * states);
void global_sums_int (int *vect, int *length);
void my_barrier (void);
void matrix_and_diag (STATE * states, STATE * states1, double * vxc, int flag);
void precond_mg (double *res, double *work1, double *work2, int istate);

void matS_cholesky_real (STATE * states);
void get_invBmat (STATE * states, double *matB);
void print_matrix (double *b, int n, int ldb);
void get_overlap_real (double *aa, int numst, int numpt,
        int lda, double *ss, int lds);
void get_cholesky_real (double *matS);
void genvlocpsi (double * psi, int st1, double * work1, double * vtot_global, STATE * states);
void genvnlpsi (double *sg_twovpsi, double *vnl,
        int dimx, int dimy, int dimz);
void get_wave(int st, STATE * states);
void add_orbit_to_wave(int st1, double scale, double * psi1, double * wave_global, STATE * states);
void print_wave(int wave_plot, STATE * states, int coarse_level);
void mg_eig (STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, double * vxc_old,
        double * vh_old);
void fsymforces (double * force, int *s, int *irg, int *irt, int *nat,
        int *ibrav, int *nsym, double * celldm, int *nr1, int *nr2,
        int *nr3);
void symrho (double * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
        int *irg, int *ftau);

void sgama (int *nrot,int *nat,double *s,double *at,double *bg,double *tau,
        int *ityp,int *nsym,int *nr1,int *nr2,int *nr3,int *irg,
        int *irt,double *ftau,double *rtau,int *npk,int *nks,double *xk,
        double *wk,double *xau,double *rau,bool *invsym, int *wflag);

void line_min_three_point (STATE *, STATE *, double, double, double *, double *,
        double *, double *, double *);

void dot_product_orbit_orbit (STATE *orbit1, STATE *orbit2, STATE
        *orbit3, double *H, double *S);

void orbit_dot_orbit (STATE * states, STATE * states1, double *Hij_row, double * Bij_row);

void app_mask (int istate, double *u, int level);
void allocate_masks (STATE * states);

void state_corner_xyz (STATE * states);

void density_orbit_X_orbit (int st1, int st2, double scale, double * psi1,
        double * psi2, double * rho_global, int mode,
        STATE * states, ORBIT_ORBIT_OVERLAP *);

char *get_symbol (int atomic_number);

void get_matB_qnm (double *Aij);
void pack_vtot_ftoc (double * vtot, double * vtot_c);
void get_qnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void pack_rho_ctof (double * rho1, double * rho_f);
void rho_nm_mat (double *Aij, double * global_mat_X);
int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz, int *ilow, int *ihi,
        int *jlow, int *jhi, int *klow, int *khi, int cdim, int pxgrid,
        int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid,
        double * xcstart, double * ycstart, double * zcstart);
void xcgga (double * rho, double * vxc, double * exc, int mode);
double radint1 (double * func, double * r, double * rab, int n);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, double * global_mat_X);
void partial_QI (int ion, double * QI_R, ION *iptr);
void partial_nlop_s (ION * iptr, double * betax, double * betay, double * betaz,
        int ip, fftw_plan p1, fftw_plan p2);
void partial_nlop_p (ION * iptr, double * betax, double * betay, double * betaz,
        int ip, fftw_plan p1, fftw_plan p2);
void partial_nlop_d (ION * iptr, double * betax, double * betay, double * betaz,
        int ip, fftw_plan p1, fftw_plan p2);
void get_mat_Omega (STATE * states, double Omega[]);
void md_fastrelax(void);
void change_states_crds (STATE * states);

/*end shuchun wang */

FILE *open_xbs_movie (char *filename);



void correct_res (STATE *, double *);
void get_state_to_proc (STATE * states);



/* different methods to update orbitals */
void kain (int step, int N, double *xm, double *fm, int NsavedSteps);
void pulay (int step, int N, double *xm, double *fm, int NsavedSteps,
        int preconditioning);
void sd (int step, int N, double *xm, double *fm);
void pulay_mixing (int size, double *rho, int NsavedSteps);
void charge_pulay (int step, int N, double *xm, double *fm, int NsavedSteps);
void pulay_kain (int step, int N, double *xm, double *fm, int NsavedSteps);

/*  Moving center stuff  */
void get_orbit_center (STATE *state, double * x, double * y, double * z);
void update_orbit_centers (STATE * states);
int if_update_centers (STATE * states);

void write_states_info (char *outfile, STATE * states);
void read_states_info (char *outfile, STATE * states);

double get_gamma_precond (double *vtot, double small_eig);
void init_state_size (STATE * states);
void diag_eig_maitrix(double *, double *, int *);



void filter_potential (double *potential, double *r, int rg_points, double rmax, double offset, double parm, double* potential_lgrid, 
        double *rab, int l_value, double dr, double  gwidth, int lgrid_points, double rcut, double rwidth, double * drpotential_lgrid);


double dot_product_orbit_nl (STATE *st1, int ion2, double * psi,
        double * prjptr, ION_ORBIT_OVERLAP *);

void non_zero_pairs ();
void non_zero_pairs1 ();

void init_nl_xyz ();

void theta_phi_new (int st1, int st2, double theta_ion, double * st2_psi,
        double * state1_psi, int mode, STATE * states);

void print_status (STATE *, double *, double *, double *, double *, char *);
void print_state_projections (STATE *, char);
void print_global_function (double *, char, char *);
void print_state_sum (STATE * states);
void print_state (STATE * state);
void print_sum (int size, double *data, char *msg);
void print_sum_square (int size, double *data, char *msg);


void init_comm (STATE * states);

void init_comm_res (STATE * states);


void get_orbit_overlap_region (STATE * states);

void get_ion_orbit_overlap_nl (STATE * states);

void duplicate_states_info (STATE * states, STATE * states1);


void get_ion_ion_overlap_region_orbit ();


void is_state_overlap (STATE * states, char *);
void init_wf_lcao(STATE * states);
void write_data(char *name, double *, double *vh, double *vxc, double *vh_old,
        double *vxc_old, double *rho, STATE * states);

void init(double * vh, double * rho, double *rho_oppo,  double * rhocore, double * rhoc,
          STATE * states, STATE * states1, double * vnuc, double * vxc, double * vh_old, double * vxc_old);
void get_te (double *rho, double *rho_oppo, double *rhocore, double *rhoc, double *vh, double *vxc,
             STATE *states, int ii_flag);

void Scf_on(STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rho_oppo, double *rhoc, double *rhocore,
        double * vxc_old, double * vh_old, int *CONVERGENCE);
void pulay_rho_on (int step0, int N, double *xm, double *fm, int NsavedSteps,
                int Nrefresh, double scale, int preconditioning);



#ifdef __cplusplus
}
#endif



