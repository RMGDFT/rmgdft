void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);
rmg_double_t app_cil_orbital (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz,
              rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil_orbital6 (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz,
              rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void app_cir_fourth (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
void app_cir_sixth (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz);
rmg_double_t app_cil_fourth (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx,
              rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cil_sixth (rmg_double_t *psi, rmg_double_t *b, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void trade_imagesx (rmg_double_t *f, rmg_double_t *w, int dimx, int dimy, int dimz, int images, int type);
void constrain( void );
void init_pe ( int image );
void mg_restrict_6 (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz, int grid_ratio);

void fill_orbit_borders4(rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
int open_wave_file (char *filename);
void get_vtot_psi (rmg_double_t * vtot_psi, rmg_double_t * vtot, int grid_ratio);
void init_nonlocal_comm();
void normalize_orbits( STATE *states);
void init_dimension();
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


void global_to_distribute(rmg_double_t * global_array, rmg_double_t * distr_array);

void diag_eig_matrix(double *eig_matrix, double *eig, int *desca);

rmg_double_t mask_function (rmg_double_t x);
void apply_mask_function (rmg_double_t *f, rmg_double_t * r, int rg_points, rmg_double_t rmax, rmg_double_t offset);
void backout_mask_function (rmg_double_t *f, rmg_double_t dr, int lg_points, rmg_double_t rmax);
void scf_barrier_wait(void);
int is_loop_over_states(void);
void cgen_prolong(rmg_double_t coef[], rmg_double_t fraction, int order);
void xbsmovie (FILE *movie);
void find_phase (int nldim, rmg_double_t *nlcdrs, rmg_double_t *phase_sin,
                 rmg_double_t *phase_cos);

void corlyp (rmg_double_t *dp, rmg_double_t *dm, rmg_double_t *dp1, rmg_double_t *dm1, rmg_double_t *dp2, rmg_double_t *dm2, rmg_double_t *ec,
             rmg_double_t *vcp0, rmg_double_t *vcm0, int *ndm);

void app_grad6(double * f, double * wx, double * wy, double * wz, int dimx, int dimy, int dimz);


STATE *init_states();

void mgrid_solv_local(rmg_double_t * v_mat, rmg_double_t * f_mat, rmg_double_t * work,
                int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                int *post_cyc, int mu_cyc, int istate, int *iion, int flag_local);

/* Function prototypes */
void app_4del2 (rmg_double_t * f, rmg_double_t * work);

void app_nl (rmg_double_t * psiR, rmg_double_t * psiI, rmg_double_t * workR, rmg_double_t * workI, int state,
        int flag, int kidx, int tid);
void cholesky (rmg_double_t * a, int n);
void cross_product (rmg_double_t * a, rmg_double_t * b, rmg_double_t * c);
rmg_double_t fill (STATE * states, rmg_double_t width, rmg_double_t nel, rmg_double_t mix,
        int num_st, int occ_flag);
void genvpsi (rmg_double_t * psi, rmg_double_t * twovpsi, rmg_double_t * pvtot, rmg_double_t * pvnl,
        rmg_double_t * kd, rmg_double_t kmag, int dimx, int dimy, int dimz);
void get_index_loc (STATE *);
void get_eig (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc);
char *get_num (char *str);
void get_zdens (STATE * states, int state, rmg_double_t * zvec);
void xcgga (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc, int flag);
void gram (KPOINT * kpoint, rmg_double_t h, int numst, int maxst, int numpt,
        int maxpt);
rmg_double_t get_ke (STATE * sp, int tid);
void global_sums (rmg_double_t * vect, int *length, MPI_Comm comm);
void init_kbr (void);
void init_pe_on (void);
void init_pegrid (void);
void init_wf (STATE * states);
void init_wflcao (STATE * states);
void init_nuc (rmg_double_t * vnuc, rmg_double_t * rhoc, rmg_double_t * rhocore);
void init_pos ();
void init_sym (void);
void symmetrize_rho (rmg_double_t * rho);
void symforce (void);
void rmg_timings (int what, rmg_double_t time);
void mg_eig_state (STATE * states, int tid, rmg_double_t * vtot);
rmg_double_t minimage (ION * ip1, ION * ip2, rmg_double_t * xtal_r);
rmg_double_t my_crtc (void);
/* Local function prototypes */
void norm_psi (STATE * psi, rmg_double_t bv);
void ortho_full (STATE * states);
void ortho_half (STATE * states);
void ortho_bcc (STATE * states);
void output (STATE * states, int its);
void pe2xyz (int pe, int *x, int *y, int *z);
void pack_ptos (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
void pack_stop (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (rmg_double_t * sg, rmg_double_t * pg, rmg_double_t alpha, int dimx, int dimy,
        int dimz);
void pack_ptos_trade (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
void pack_vhstod (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz);
void pack_vhdtos (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz);
double radint (double *f, double *r, int n, double al);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
rmg_double_t real_sum_all (rmg_double_t x, MPI_Comm comm);
rmg_double_t real_max_all (rmg_double_t x);
void sortpsi (STATE * states);
void subdiag (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
void trade_images (rmg_double_t * mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_mpi (rmg_double_t * mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_images_smp (rmg_double_t * mat, int dimx, int dimy, int dimz, int *nb_ids);
void set_bc (rmg_double_t * mat, int dimx, int dimy, int dimz, int images, rmg_double_t val);
void getpoi_bc (rmg_double_t * rho, rmg_double_t * vh_bc, int dimx, int dimy, int dimz);
void vol_rho (rmg_double_t * rho, int step);
void vol_wf (STATE * states, int state, int step);
void write_avgd (rmg_double_t * rho);
void write_avgv (rmg_double_t * vh, rmg_double_t * vnuc);
void write_zstates (STATE * states);
void write_header (void);
void write_pos (void);
void write_eigs (STATE * states);
void write_occ (STATE * states);
void write_timings (void);
rmg_double_t rand0 (long *idum);

void gather_psi (rmg_double_t * tmp_psiR, rmg_double_t * tmp_psiI, STATE * sp, int tid);
void scatter_psi (rmg_double_t * tmp_psiR, rmg_double_t * tmp_psiI, STATE * sp, int tid);
void get_milliken (STATE * states);

void output_wave (STATE * states, int kpt, int fhand);



/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, rmg_double_t * celldm, rmg_double_t * A0I, rmg_double_t * A1I, rmg_double_t * A2I,
        rmg_double_t * OMEGAI, int *flag);
void recips (void);
void to_cartesian (rmg_double_t crystal[], rmg_double_t cartesian[]);
void to_crystal (rmg_double_t crystal[], rmg_double_t cartesian[]);
rmg_double_t metric (rmg_double_t * crystal);

/* Md run types */
void run (STATE * states, STATE * states1);
void fastrlx( STATE *states, STATE *states1, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc, rmg_double_t *vh_old, rmg_double_t *vxc_old,
        rmg_double_t *rho, rmg_double_t *rhocore, rmg_double_t *rhoc );
void cdfastrlx (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc,
        rmg_double_t * rho, rmg_double_t * rhocore, rmg_double_t * rhoc);
/*void moldyn (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * rho, rmg_double_t * rhoc, rmg_double_t * rhocore);*/
void dx (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * rho,
        rmg_double_t * rhoc);
void psidx (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * rho,
        rmg_double_t * rhoc);
void cholesky (rmg_double_t * a, int n);


rmg_double_t minimage1 (rmg_double_t aa[3], rmg_double_t bb[3]);
void init_parameter (STATE * states);
void get_mehr (void);
void make_mask_grid (rmg_double_t rcut, int level, STATE * states);
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
void xclsd_pz81 (rmg_double_t * rho, rmg_double_t * vxc);
void global_sums_int (int *, int *);
rmg_double_t linint (rmg_double_t * y, rmg_double_t r, rmg_double_t invdr);
void ortho_norm_local (STATE * states);
void global_sums_int (int *vect, int *length);
void my_barrier (void);
void matrix_and_diag (STATE * states, STATE * states1, rmg_double_t * vxc, int flag);
void precond_mg (double *res, double *work1, double *work2, int istate);

void matS_cholesky_real (STATE * states);
void get_invBmat (STATE * states, double *matB);
void print_matrix (double *b, int n, int ldb);
void get_overlap_real (double *aa, int numst, int numpt,
        int lda, double *ss, int lds);
void get_cholesky_real (double *matS);
void get_all_kbpsi (STATE * states1, STATE * states, ION_ORBIT_OVERLAP *, rmg_double_t *, rmg_double_t *);
void get_Hvnlij (double *Aij);
void genvlocpsi (rmg_double_t * psi, int st1, rmg_double_t * work1, rmg_double_t * vtot_global, STATE * states);
void genvnlpsi (double *sg_twovpsi, double *vnl,
        int dimx, int dimy, int dimz);
void get_new_rho (STATE * states, double *rho);
void get_wave(int st, STATE * states);
void add_orbit_to_wave(int st1, rmg_double_t scale, rmg_double_t * psi1, rmg_double_t * wave_global, STATE * states);
void print_wave(int wave_plot, STATE * states, int coarse_level);
void mg_eig (STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, rmg_double_t * vxc_old,
        rmg_double_t * vh_old);
void fsymforces (rmg_double_t * force, int *s, int *irg, int *irt, int *nat,
        int *ibrav, int *nsym, rmg_double_t * celldm, int *nr1, int *nr2,
        int *nr3);
void symmetry (int *ibrav, int *s, int *nsym, int *irg, int *irt, int *ftau,
        int *nat, rmg_double_t * tau, int *ityp, int *nks, rmg_double_t * xk,
        rmg_double_t * wk, rmg_double_t * celldm, int *nr1, int *nr2, int *nr3,
        int *wflag);
void symrho (rmg_double_t * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
        int *irg, int *ftau);

void sgama (int *nrot,int *nat,double *s,double *at,double *bg,double *tau,
        int *ityp,int *nsym,int *nr1,int *nr2,int *nr3,int *irg,
        int *irt,double *ftau,double *rtau,int *npk,int *nks,double *xk,
        double *wk,double *xau,double *rau,bool *invsym, int *wflag);

void line_min_three_point (STATE *, STATE *, rmg_double_t, rmg_double_t, rmg_double_t *, rmg_double_t *,
        rmg_double_t *, rmg_double_t *, rmg_double_t *);

void dot_product_orbit_orbit (STATE *orbit1, STATE *orbit2, STATE
        *orbit3, double *H, double *S);

void orbit_dot_orbit (STATE * states, STATE * states1, rmg_double_t *Hij_row, rmg_double_t * Bij_row);

void app_mask (int istate, double *u, int level);
void allocate_masks (STATE * states);

void state_corner_xyz (STATE * states);

char *get_symbol (int atomic_number);

void get_matB_qnm (double *Aij);
void pack_vtot_ftoc (rmg_double_t * vtot, rmg_double_t * vtot_c);
void get_qnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void qnm_beta_betapsi (STATE *sp, int ion2, rmg_double_t * prjptr, rmg_double_t * work);
void get_dnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void dnm_beta_betapsi (STATE *sp, int ion2, rmg_double_t scale, int ip2, rmg_double_t * prjptr,
        rmg_double_t * work);
void pack_rho_ctof (rmg_double_t * rho1, rmg_double_t * rho_f);
void rho_augmented (rmg_double_t * rho, rmg_double_t * global_mat_X);
void rho_Qnm_mat (double *Aij, rmg_double_t * global_mat_X);
void rho_nm_mat (double *Aij, rmg_double_t * global_mat_X);
int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz, int *ilow, int *ihi,
        int *jlow, int *jhi, int *klow, int *khi, int cdim, int pxgrid,
        int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid,
        rmg_double_t * xcstart, rmg_double_t * ycstart, rmg_double_t * zcstart);
void xcgga (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc, int mode);
rmg_double_t radint1 (rmg_double_t * func, rmg_double_t * r, rmg_double_t * rab, int n);
void get_all_partial_kbpsi (STATE *states, ION_ORBIT_OVERLAP *,
        rmg_double_t *, rmg_double_t *, rmg_double_t *, rmg_double_t *,
        rmg_double_t *, rmg_double_t *);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, rmg_double_t * global_mat_X);
void partial_QI (int ion, rmg_double_t * QI_R, ION *iptr);
void partial_nlop_s (ION * iptr, rmg_double_t * betax, rmg_double_t * betay, rmg_double_t * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void partial_nlop_p (ION * iptr, rmg_double_t * betax, rmg_double_t * betay, rmg_double_t * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void partial_nlop_d (ION * iptr, rmg_double_t * betax, rmg_double_t * betay, rmg_double_t * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void get_mat_Omega (STATE * states, double Omega[]);
void md_fastrelax(void);
void change_states_crds (STATE * states);

/*end shuchun wang */

FILE *open_xbs_movie (char *filename);






rmg_double_t get_te_ion_ion ();
rmg_double_t get_sum_eig (STATE * states);
rmg_double_t get_Exc (rmg_double_t * rho, rmg_double_t * rhocore);
void correct_res (STATE *, rmg_double_t *);
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
void get_orbit_center (STATE *state, rmg_double_t * x, rmg_double_t * y, rmg_double_t * z);
void update_orbit_centers (STATE * states);
int if_update_centers (STATE * states);

void write_states_info (char *outfile, STATE * states);
void read_states_info (char *outfile, STATE * states);

void init_state_size (STATE * states);
void diag_eig_maitrix(double *, double *, int *);



void filter_potential (rmg_double_t *potential, rmg_double_t *r, int rg_points, rmg_double_t rmax, rmg_double_t offset, rmg_double_t parm, rmg_double_t* potential_lgrid, 
        rmg_double_t *rab, int l_value, rmg_double_t dr, rmg_double_t  gwidth, int lgrid_points, rmg_double_t rcut, rmg_double_t rwidth, rmg_double_t * drpotential_lgrid);


/* Function prototypes */
void app_4del2 (rmg_double_t * f, rmg_double_t * work);

void app_nl (rmg_double_t * psiR, rmg_double_t * psiI, rmg_double_t * workR, rmg_double_t * workI, int state,
        int flag, int kidx, int tid);
void cholesky (rmg_double_t * a, int n);
void cross_product (rmg_double_t * a, rmg_double_t * b, rmg_double_t * c);
rmg_double_t fill (STATE * states, rmg_double_t width, rmg_double_t nel, rmg_double_t mix,
        int num_st, int occ_flag);
void genvpsi (rmg_double_t * psi, rmg_double_t * twovpsi, rmg_double_t * pvtot, rmg_double_t * pvnl,
        rmg_double_t * kd, rmg_double_t kmag, int dimx, int dimy, int dimz);
void get_index_loc (STATE *);
void get_eig (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc);
char *get_num (char *str);
void get_zdens (STATE * states, int state, rmg_double_t * zvec);
void xcgga (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc, int flag);
void gram (KPOINT * kpoint, rmg_double_t h, int numst, int maxst, int numpt,
        int maxpt);
rmg_double_t get_ke (STATE * sp, int tid);
void global_sums (rmg_double_t * vect, int *length, MPI_Comm comm);
void init_kbr (void);
void init_pe_on (void);
void init_pegrid (void);
void init_wf (STATE * states);
void init_wflcao (STATE * states);
void init_nuc (rmg_double_t * vnuc, rmg_double_t * rhoc, rmg_double_t * rhocore);
void init_pos ();
void init_sym (void);
void symmetrize_rho (rmg_double_t * rho);
void symforce (void);
void rmg_timings (int what, rmg_double_t time);
void mg_eig_state (STATE * states, int tid, rmg_double_t * vtot);
rmg_double_t minimage (ION * ip1, ION * ip2, rmg_double_t * xtal_r);
rmg_double_t my_crtc (void);
/* Local function prototypes */
void norm_psi (STATE * psi, rmg_double_t bv);
void ortho_full (STATE * states);
void ortho_half (STATE * states);
void ortho_bcc (STATE * states);
void output (STATE * states, int its);
void pe2xyz (int pe, int *x, int *y, int *z);
void pack_ptos (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
void pack_stop (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
void pack_stop_axpy (rmg_double_t * sg, rmg_double_t * pg, rmg_double_t alpha, int dimx, int dimy,
        int dimz);
void pack_ptos_trade (rmg_double_t * sg, rmg_double_t * pg, int dimx, int dimy, int dimz);
void pack_vhstod (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz);
void pack_vhdtos (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz);
double radint (double *f, double *r, int n, double al);
void radiff (double *f, double *df, double *r, int n, double al);
void ra2diff (double *f, double *df, double *r, int n, double al);
void ranv (void);
rmg_double_t real_sum_all (rmg_double_t x, MPI_Comm comm);
rmg_double_t real_max_all (rmg_double_t x);
void sortpsi (STATE * states);
void subdiag (STATE * states, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vxc);
void trade_images (rmg_double_t * mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_mpi (rmg_double_t * mat, int dimx, int dimy, int dimz, int *nb_ids);
void trade_images_smp (rmg_double_t * mat, int dimx, int dimy, int dimz, int *nb_ids);
void set_bc (rmg_double_t * mat, int dimx, int dimy, int dimz, int images, rmg_double_t val);
void getpoi_bc (rmg_double_t * rho, rmg_double_t * vh_bc, int dimx, int dimy, int dimz);
void vol_rho (rmg_double_t * rho, int step);
void vol_wf (STATE * states, int state, int step);
void write_avgd (rmg_double_t * rho);
void write_avgv (rmg_double_t * vh, rmg_double_t * vnuc);
void write_zstates (STATE * states);
void write_header (void);
void write_pos (void);
void write_eigs (STATE * states);
void write_occ (STATE * states);
void write_timings (void);
rmg_double_t rand0 (long *idum);

void gather_psi (rmg_double_t * tmp_psiR, rmg_double_t * tmp_psiI, STATE * sp, int tid);
void scatter_psi (rmg_double_t * tmp_psiR, rmg_double_t * tmp_psiI, STATE * sp, int tid);
void get_milliken (STATE * states);

void output_wave (STATE * states, int kpt, int fhand);



/* Conversion between crystal and cartesian coordinate prototypes */
void latgen (int *ibrav, rmg_double_t * celldm, rmg_double_t * A0I, rmg_double_t * A1I, rmg_double_t * A2I,
        rmg_double_t * OMEGAI, int *flag);
void recips (void);
void to_cartesian (rmg_double_t crystal[], rmg_double_t cartesian[]);
void to_crystal (rmg_double_t crystal[], rmg_double_t cartesian[]);
rmg_double_t metric (rmg_double_t * crystal);

/* Md run types */
void run (STATE * states, STATE * states1);
void fastrlx( STATE *states, STATE *states1, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc, rmg_double_t *vh_old, rmg_double_t *vxc_old,
        rmg_double_t *rho, rmg_double_t *rhocore, rmg_double_t *rhoc );
void cdfastrlx (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc,
        rmg_double_t * rho, rmg_double_t * rhocore, rmg_double_t * rhoc);
/*void moldyn (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * rho, rmg_double_t * rhoc, rmg_double_t * rhocore);*/
void dx (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * rho,
        rmg_double_t * rhoc);
void psidx (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * rho,
        rmg_double_t * rhoc);
void cholesky (rmg_double_t * a, int n);


rmg_double_t minimage1 (rmg_double_t aa[3], rmg_double_t bb[3]);
void init_parameter (STATE * states);
void get_mehr (void);
void make_mask_grid (rmg_double_t rcut, int level, STATE * states);
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
void xclsd_pz81 (rmg_double_t * rho, rmg_double_t * vxc);
void global_sums_int (int *, int *);
rmg_double_t linint (rmg_double_t * y, rmg_double_t r, rmg_double_t invdr);
void ortho_norm_local (STATE * states);
void global_sums_int (int *vect, int *length);
void my_barrier (void);
void matrix_and_diag (STATE * states, STATE * states1, rmg_double_t * vxc, int flag);
void precond_mg (double *res, double *work1, double *work2, int istate);

void matS_cholesky_real (STATE * states);
void get_invBmat (STATE * states, double *matB);
void print_matrix (double *b, int n, int ldb);
void get_overlap_real (double *aa, int numst, int numpt,
        int lda, double *ss, int lds);
void get_cholesky_real (double *matS);
void get_Hvnlij (double *Aij);
void genvlocpsi (rmg_double_t * psi, int st1, rmg_double_t * work1, rmg_double_t * vtot_global, STATE * states);
void genvnlpsi (double *sg_twovpsi, double *vnl,
        int dimx, int dimy, int dimz);
void get_new_rho (STATE * states, double *rho);
void get_wave(int st, STATE * states);
void add_orbit_to_wave(int st1, rmg_double_t scale, rmg_double_t * psi1, rmg_double_t * wave_global, STATE * states);
void print_wave(int wave_plot, STATE * states, int coarse_level);
void mg_eig (STATE * states, STATE * states1, double *vxc, double *vh,
        double *vnuc, double *rho, double *rhoc, rmg_double_t * vxc_old,
        rmg_double_t * vh_old);
void fsymforces (rmg_double_t * force, int *s, int *irg, int *irt, int *nat,
        int *ibrav, int *nsym, rmg_double_t * celldm, int *nr1, int *nr2,
        int *nr3);
void symmetry (int *ibrav, int *s, int *nsym, int *irg, int *irt, int *ftau,
        int *nat, rmg_double_t * tau, int *ityp, int *nks, rmg_double_t * xk,
        rmg_double_t * wk, rmg_double_t * celldm, int *nr1, int *nr2, int *nr3,
        int *wflag);
void symrho (rmg_double_t * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
        int *irg, int *ftau);

void sgama (int *nrot,int *nat,double *s,double *at,double *bg,double *tau,
        int *ityp,int *nsym,int *nr1,int *nr2,int *nr3,int *irg,
        int *irt,double *ftau,double *rtau,int *npk,int *nks,double *xk,
        double *wk,double *xau,double *rau,bool *invsym, int *wflag);

void line_min_three_point (STATE *, STATE *, rmg_double_t, rmg_double_t, rmg_double_t *, rmg_double_t *,
        rmg_double_t *, rmg_double_t *, rmg_double_t *);

void dot_product_orbit_orbit (STATE *orbit1, STATE *orbit2, STATE
        *orbit3, double *H, double *S);

void orbit_dot_orbit (STATE * states, STATE * states1, rmg_double_t *Hij_row, rmg_double_t * Bij_row);

void app_mask (int istate, double *u, int level);
void allocate_masks (STATE * states);

void state_corner_xyz (STATE * states);

void density_orbit_X_orbit (int st1, int st2, rmg_double_t scale, rmg_double_t * psi1,
        rmg_double_t * psi2, rmg_double_t * rho_global, int mode,
        STATE * states, ORBIT_ORBIT_OVERLAP *);

char *get_symbol (int atomic_number);

void get_matB_qnm (double *Aij);
void pack_vtot_ftoc (rmg_double_t * vtot, rmg_double_t * vtot_c);
void get_qnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void qnm_beta_betapsi (STATE *sp, int ion2, rmg_double_t * prjptr, rmg_double_t * work);
void get_dnmpsi (STATE *sp, double *kbpsi_one_state, double *work);
void dnm_beta_betapsi (STATE *sp, int ion2, rmg_double_t scale, int ip2, rmg_double_t * prjptr,
        rmg_double_t * work);
void pack_rho_ctof (rmg_double_t * rho1, rmg_double_t * rho_f);
void rho_augmented (rmg_double_t * rho, rmg_double_t * global_mat_X);
void rho_Qnm_mat (double *Aij, rmg_double_t * global_mat_X);
void rho_nm_mat (double *Aij, rmg_double_t * global_mat_X);
int get_index (int gridpe, ION * iptr, int *Aix, int *Aiy, int *Aiz, int *ilow, int *ihi,
        int *jlow, int *jhi, int *klow, int *khi, int cdim, int pxgrid,
        int pygrid, int pzgrid, int nxgrid, int nygrid, int nzgrid,
        rmg_double_t * xcstart, rmg_double_t * ycstart, rmg_double_t * zcstart);
void xcgga (rmg_double_t * rho, rmg_double_t * vxc, rmg_double_t * exc, int mode);
rmg_double_t radint1 (rmg_double_t * func, rmg_double_t * r, rmg_double_t * rab, int n);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, rmg_double_t * global_mat_X);
void partial_QI (int ion, rmg_double_t * QI_R, ION *iptr);
void partial_nlop_s (ION * iptr, rmg_double_t * betax, rmg_double_t * betay, rmg_double_t * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void partial_nlop_p (ION * iptr, rmg_double_t * betax, rmg_double_t * betay, rmg_double_t * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void partial_nlop_d (ION * iptr, rmg_double_t * betax, rmg_double_t * betay, rmg_double_t * betaz,
        int ip, fftwnd_plan p1, fftwnd_plan p2);
void get_mat_Omega (STATE * states, double Omega[]);
void md_fastrelax(void);
void change_states_crds (STATE * states);

/*end shuchun wang */

FILE *open_xbs_movie (char *filename);



void correct_res (STATE *, rmg_double_t *);
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
void get_orbit_center (STATE *state, rmg_double_t * x, rmg_double_t * y, rmg_double_t * z);
void update_orbit_centers (STATE * states);
int if_update_centers (STATE * states);

void write_states_info (char *outfile, STATE * states);
void read_states_info (char *outfile, STATE * states);

double get_gamma_precond (double *vtot, double small_eig);
void init_state_size (STATE * states);
void diag_eig_maitrix(double *, double *, int *);



void filter_potential (rmg_double_t *potential, rmg_double_t *r, int rg_points, rmg_double_t rmax, rmg_double_t offset, rmg_double_t parm, rmg_double_t* potential_lgrid, 
        rmg_double_t *rab, int l_value, rmg_double_t dr, rmg_double_t  gwidth, int lgrid_points, rmg_double_t rcut, rmg_double_t rwidth, rmg_double_t * drpotential_lgrid);


rmg_double_t dot_product_orbit_nl (STATE *st1, int ion2, rmg_double_t * psi,
        rmg_double_t * prjptr, ION_ORBIT_OVERLAP *);

void non_zero_pairs ();
void non_zero_pairs1 ();

void init_nl_xyz ();

void theta_phi_new (int st1, int st2, rmg_double_t theta_ion, rmg_double_t * st2_psi,
        rmg_double_t * state1_psi, int mode, STATE * states);

void print_status (STATE *, rmg_double_t *, rmg_double_t *, rmg_double_t *, rmg_double_t *, char *);
void print_state_projections (STATE *, char);
void print_global_function (rmg_double_t *, char, char *);
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


void is_state_overlap (STATE * states);

