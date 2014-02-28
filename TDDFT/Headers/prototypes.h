void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);
void mgrid_solv (REAL *v_mat, REAL *f_mat, REAL *work,
                 int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy,
                 REAL gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, REAL step, REAL k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);
REAL app_cil_orbital (REAL * a, REAL * b, int dimx, int dimy, int dimz,
              REAL gridhx, REAL gridhy, REAL gridhz);
REAL app_cil_orbital6 (REAL * a, REAL * b, int dimx, int dimy, int dimz,
              REAL gridhx, REAL gridhy, REAL gridhz);
void mg_prolong_6 (double *full, double *half, int dimx, int dimy, int dimz, int half_dimx,
                       int half_dimy, int half_dimz, int grid_ratio, int order);
void app_cir_fourth (REAL *a, REAL *b, int dimx, int dimy, int dimz);
void app_cir_sixth (REAL *a, REAL *b, int dimx, int dimy, int dimz);
REAL app_cil_fourth (REAL *a, REAL *b, int dimx, int dimy, int dimz, REAL gridhx,
              REAL gridhy, REAL gridhz);
REAL app_cil_sixth (REAL *psi, REAL *b, int dimx, int dimy, int dimz,
                    REAL gridhx, REAL gridhy, REAL gridhz);
void trade_imagesx (REAL *f, REAL *w, int dimx, int dimy, int dimz, int images, int type);
void constrain( void );
void init_pe ( int image );
void mg_restrict_6 (REAL * full, REAL * half, int dimx, int dimy, int dimz, int grid_ratio);

void fill_orbit_borders4(REAL * sg, REAL * pg, int dimx, int dimy, int dimz);
int assign_species (CONTROL * c, char *buf);

int open_wave_file (char *filename);
void get_vtot_psi (REAL * vtot_psi, REAL * vtot, int grid_ratio);
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


void sl_init(int *ictxt, int nprow, int npcol);

void init_rho_atom(double *rho);

void init_wf_gaussian(STATE * states);


void mg_prolong_MAX10 (double * full, double * half, int dimx, int dimy, int dimz, int half_dimx, int half_dimy, int half_dimz, int grid_ratio, int order);


void global_to_distribute(REAL * global_array, REAL * distr_array);

void diag_eig_matrix(double *eig_matrix, double *eig, int *desca);
void init_weight (void);
void init_weight_s (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftwnd_plan p1);
void init_weight_p (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftwnd_plan p1);
void init_weight_d (SPECIES *sp, fftw_complex *rtptr, int ip,
                    fftwnd_plan p1);

REAL mask_function (REAL x);
void apply_mask_function (REAL *f, REAL * r, int rg_points, REAL rmax, REAL offset);
void backout_mask_function (REAL *f, REAL dr, int lg_points, REAL rmax);
void thread_barrier_wait(void);
int is_loop_over_states(void);
void cgen_prolong(REAL coef[], REAL fraction, int order);
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
void xbsmovie (FILE *movie);
void find_phase (int nldim, REAL *nlcdrs, REAL *phase_sin,
                 REAL *phase_cos);

void corlyp (REAL *dp, REAL *dm, REAL *dp1, REAL *dm1, REAL *dp2, REAL *dm2, REAL *ec,
             REAL *vcp0, REAL *vcm0, int *ndm);

void app_grad6(double * f, double * wx, double * wy, double * wz, int dimx, int dimy, int dimz);


STATE *init_states();

void mgrid_solv_local(REAL * v_mat, REAL * f_mat, REAL * work,
                int dimx, int dimy, int dimz,
                REAL gridhx, REAL gridhy, REAL gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                int *post_cyc, int mu_cyc, int istate, int *iion, int flag_local);

