/*  
    $Id$
*/



void force(REAL *rho, REAL *rho_opps,REAL *rhoc, REAL *vh, REAL *vxc, REAL *vnuc,
            STATE *states);
void lforce(REAL *rho, REAL *vh);
void nlccforce(REAL *rho, REAL *vxc);
void nlforce(double *);
void nlforce_par_D (STATE *states, REAL *forces);

REAL *vnuc_x, *vnuc_y, *vnuc_z;
short int vloc_state_overlap_or_not[MAX_IONS * MAX_STATES];

void multi_GHG_munu (REAL *GHG_tri, REAL *GHG_en_tri);
double trace_AB_tri(REAL * H_tri, REAL * par_D_tri, int N, int *ni);

void tri_to_whole_complex (complex double * H_tri, complex double * Hii, int N, int * ni);
void whole_to_tri_complex (complex double * H_tri, complex double * Hii, int N, int * ni);
void partial_vlocpsi (STATE st1, int ion2, REAL * psi, REAL * prjptr, REAL *vlpsi);
void ylmr2_x(double *r, double *ylm_x);
void ylmr2_y(double *r, double *ylm_x);
void ylmr2_z(double *r, double *ylm_x);
void qval_R(int ih, int jh, REAL r, REAL * x, REAL * qlig, REAL * drqlig, REAL invdr, int *nhtol,
            int *nhtom, int *indv, REAL * ylm, REAL * ylm_x, REAL * ylm_y, REAL * ylm_z,
            REAL ap[][9][9], int lpx[][9], int lpl[][9][9], REAL * Q_x, REAL * Q_y, REAL * Q_z,
            SPECIES * sp);
void ion_partial_Hij_and_Sij (int iion, int flag,  double *Hij, double *Sij);
void tri_to_whole_real (REAL * H_tri, REAL * Hii, int N, int * ni);
void nlforce_partial_H_part2 (STATE * states, STATE * states1, REAL *GHG, REAL *force);
void get_all_partial_kbpsi (STATE * states);
void rho_nm_mat (double *Aij, REAL * global_mat_X);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, REAL * global_mat_X);
void partial_QI(int ion, REAL * QI_R, ION *iptr);
void get_partial_ddd (REAL * QnmI_R, REAL * veff, int ion, int nh);
void nlforce_par_Q(REAL * veff, REAL * rho_nm, int ion, int nh, double *forces);
void nlforce_par_rho(REAL * par_gamma_x, REAL * par_gamma_y, REAL * par_gamma_z, int ion, int nh);
void is_vloc_state_overlap (STATE *states);
void get_ion_orbit_overlap_loc (STATE * states);
void read_trans (complex_energy_integral * cei);
void get_dos (STATE * states); 
void get_3Ddos (STATE * states); 

void Stransfer_p (complex double * tot, complex double * tott, complex
        double *ch0, complex double *ch01, complex double *ch10, int iprobe);
void Sigma_p (complex double *sigma, complex double *ch0, complex double
        *ch01, complex double *ch10, complex double *green, int iprobe);
void Sgreen_cond_p (complex double *H_tri, complex double *sigma_all, int *sigma_idx,
        complex double *Green_C_leftdown, int nC, int iprobe1, int iprobe2);

void   comm_sums(double *, int*, MPI_Comm);
//void Sgreen_p (complex double * tot, complex double * tott, complex double
//        *ch0, complex double *ch1, complex double *g, int iprobe);
void gauleg (double x1, double x2, double *x, double *w, int n);
void fill_orbit_borders (REAL * sg, int dimx, int dimy, int dimz);
void get_distributed_mat (double *bigmat, double *dismat);
void whole_to_tri_p (REAL * A_tri, REAL * Aii, int N, int *ni);
void pmo_unitary_matrix(complex double *a_local, int *desca);
/*void matrix_inverse_cond_p (complex double * H_tri, int N, int * ni,
  complex double * Green_C_leftdown);*/
int int_max_all (int x);
/*void read_potrho(char *name, double *potrho, int NX, int NY, int NZ, 
  int PNX0, int PNY0, int PNZ0, int PNX1, int PNY1, int PNZ1, 
  int data_indicator, double x0_old, double x_start, double hx_old, double hx_new);*/
void read_potrho (double *vh, int iflag, int data_indicator);
void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y);
void sigma_all_energy_point (complex double * sigma_all);
void diff_hx_interpolation (double *xi, double *xi_old, int NX,
        double hx, double hx_old, double x0, double x0_old);
void find_phase (int nldim, REAL * nlcdrs, REAL * phase_sin, REAL * phase_cos);
void read_orbital (STATE * states);

void read_potrho_LCR (double *vh, double *vxc, double *rho);
void write_rho_x (REAL * rho, char *ab);
void interpolation_orbit (STATE * states);
void tri_to_whole_p (REAL * A_tri, REAL * Aii, int N, int *ni);
void tri_to_whole_complex_p (complex double * A_tri, complex double * Aii, int N, int *ni);
void global_to_distribute (REAL * global_array, REAL * distr_array);
void global_to_distribute2 (double *global_array, double *distr_array);
void global_to_distribute3 (double *global_array, double *distr_array);
void Sgreen_semi_infinite_p (complex double * green, complex double
        *ch00, complex double *ch01, int jprobe);
void find_fermi (complex double * sigma_all);
void charge_density_matrix_p (complex double * sigma_all);
void pulay_rho (int step0, int N, double *xm, double *fm, int NsavedSteps,
        int Nrefresh, double scale, int preconditioning);
void get_inverse_block_p (complex double *Hii, complex double *Gii, int *ipiv, int *desca );
void init_weight_s (SPECIES *sp, fftw_complex *rtptr, int ip, fftwnd_plan p1);
void init_weight_p (SPECIES *sp, fftw_complex *rtptr, int ip, fftwnd_plan p1);
void init_weight_d (SPECIES *sp, fftw_complex *rtptr, int ip, fftwnd_plan p1);
void confine (REAL * mat, int size_x, int size_y, int size_z, COMPASS compass, int grid_type);
void distribute_to_global (REAL * distr_array, REAL * global_array);
void distribute_to_global2 (double *distr_array, double *global_array);
void distribute_to_X_soft (REAL * distr_array, REAL * global_array);
void distribute_to_Y_soft (REAL * distr_array, REAL * global_array);
void X_to_distribute_soft (REAL * global_array, REAL * distr_array);
void Y_to_distribute_soft (REAL * global_array, REAL * distr_array);
void Sgreen_c_p (REAL * Htri, REAL * Stri, complex double * sigma_L, int *sigma_idx,
        complex double  ene, complex double * Green_C);
void Sgreen_c_noneq_p (double *H00, double *S00, complex double * sigma_L,
        int *sigma_idx, complex double ene, complex double * Green_C, int nC,
        int iprobe);
void rho_munu_p (complex double * rho_mn, complex double * green_C, complex double * sigma_L, int iprobe);
void init_derweight_s (SPECIES *sp, 
        fftw_complex *rtptr_x, 
        fftw_complex *rtptr_y,
        fftw_complex *rtptr_z, 
        int ip, fftwnd_plan p1);
void init_derweight_p (SPECIES *sp, 
        fftw_complex *rtptr_x, 
        fftw_complex *rtptr_y,
        fftw_complex *rtptr_z, 
        int ip, fftwnd_plan p1);
void init_derweight_d (SPECIES *sp, 
        fftw_complex *rtptr_x, 
        fftw_complex *rtptr_y,
        fftw_complex *rtptr_z, 
        int ip, fftwnd_plan p1);
void matrix_inverse_p (complex double * H_tri, complex double * G_tri);


void matrix_inverse_anyprobe_p (complex double *H_tri, int N, int *ni, int iprobe, complex double *Green_C);
void matrix_inverse_left_p (complex double * H_tri, int N, int * ni, complex double * Green_C);
void matrix_inverse_right_p (complex double * H_tri, int N, int * ni, complex double * Green_C);


void global_sums_X (REAL * vect, int *length);
void corlyp (REAL *dp, REAL *dm, REAL *dp1, REAL *dm1, REAL *dp2,
        REAL *dm2, REAL *ec, REAL *vcp0, REAL *vcm0, int *ndm);

#if (XT3)
void init_loc_xyz ();
void write_force();
void iiforce();
void partial_vloc ();
void init_dimension();
void init_pos();
void pmo_init ();
void lead_bandstructure ();
void get_cond_frommatrix ();
void allocate_matrix_soft ();
void writeout_matrix_p ();
void read_matrix_pp ();
void init_weight ();
void init_derweight ();
void allocate_matrix_LCR ();
void read_lead_matrix ();
void get_qqq ();
#else
void init_loc_xyz (void);
void write_force(void);
void iiforce(void);
void partial_vloc (void);
void init_dimension(void);
void init_pos(void);
void pmo_init (void);
void lead_bandstructure (void);
void get_cond_frommatrix (void);
void allocate_matrix_soft (void);
void writeout_matrix_p (void);
void read_matrix_pp (void);
void init_weight (void);
void init_derweight (void);
void allocate_matrix_LCR (void);
void read_lead_matrix (void);
void get_qqq (void);
#endif


void filter_potential (REAL *potential, REAL *r, int rg_points, REAL rmax, REAL offset, REAL parm, REAL* potential_lgrid, 
	REAL *rab, int l_value, REAL dr, REAL  gwidth, int lgrid_points, REAL rcut, REAL rwidth, REAL * drpotential_lgrid);

int min_distance_index(double *distance, int n );


void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);

void mgrid_solv_negf(REAL * v_mat, REAL * f_mat, REAL * work,
                int dimx, int dimy, int dimz,
                REAL gridhx, REAL gridhy, REAL gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, REAL step, REAL k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);
void get_vh_negf (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);

double find_new_energy_point(double *cond, double *ener1, int tot_energy_point, double simpson_tol, int *EP_final, 
                           int *energy_insert_index, double *ener1_temp);

