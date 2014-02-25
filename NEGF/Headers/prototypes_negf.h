/*  
    $Id$
*/



void force(rmg_double_t *rho, rmg_double_t *rho_opps,rmg_double_t *rhoc, rmg_double_t *vh, rmg_double_t *vxc, rmg_double_t *vnuc,
            STATE *states);
void lforce(rmg_double_t *rho, rmg_double_t *vh);
void nlccforce(rmg_double_t *rho, rmg_double_t *vxc);
void nlforce(double *);
void nlforce_par_D (STATE *states, rmg_double_t *forces);

rmg_double_t *vnuc_x, *vnuc_y, *vnuc_z;
short int vloc_state_overlap_or_not[MAX_IONS * MAX_STATES];

void multi_GHG_munu (rmg_double_t *GHG_tri, rmg_double_t *GHG_en_tri);
double trace_AB_tri(rmg_double_t * H_tri, rmg_double_t * par_D_tri, int N, int *ni);

void tri_to_whole_complex (complex double * H_tri, complex double * Hii, int N, int * ni);
void whole_to_tri_complex (complex double * H_tri, complex double * Hii, int N, int * ni);
void partial_vlocpsi (STATE st1, int ion2, rmg_double_t * psi, rmg_double_t * prjptr, rmg_double_t *vlpsi);
void ylmr2_x(double *r, double *ylm_x);
void ylmr2_y(double *r, double *ylm_x);
void ylmr2_z(double *r, double *ylm_x);
void qval_R(int ih, int jh, rmg_double_t r, rmg_double_t * x, rmg_double_t * qlig, rmg_double_t * drqlig, rmg_double_t invdr, int *nhtol,
            int *nhtom, int *indv, rmg_double_t * ylm, rmg_double_t * ylm_x, rmg_double_t * ylm_y, rmg_double_t * ylm_z,
            rmg_double_t ap[][9][9], int lpx[][9], int lpl[][9][9], rmg_double_t * Q_x, rmg_double_t * Q_y, rmg_double_t * Q_z,
            SPECIES * sp);
void ion_partial_Hij_and_Sij (int iion, int flag,  double *Hij, double *Sij);
void tri_to_whole_real (rmg_double_t * H_tri, rmg_double_t * Hii, int N, int * ni);
void nlforce_partial_H_part2 (STATE * states, STATE * states1, rmg_double_t *GHG, rmg_double_t *force);
void get_all_partial_kbpsi (STATE * states);
void rho_nm_mat (double *Aij, rmg_double_t * global_mat_X);
void partial_Mat_nm_R (double *partial_x, double *partial_y, double *partial_z, rmg_double_t * global_mat_X);
void partial_QI(int ion, rmg_double_t * QI_R, ION *iptr);
void get_partial_ddd (rmg_double_t * QnmI_R, rmg_double_t * veff, int ion, int nh);
void nlforce_par_Q(rmg_double_t * veff, rmg_double_t * rho_nm, int ion, int nh, double *forces);
void nlforce_par_rho(rmg_double_t * par_gamma_x, rmg_double_t * par_gamma_y, rmg_double_t * par_gamma_z, int ion, int nh);
void is_vloc_state_overlap (STATE *states);
void get_ion_orbit_overlap_loc (STATE * states);
void read_trans (complex_energy_integral * cei);
void get_dos (STATE * states); 
void get_3Ddos (STATE * states, double EMIN, double EMAX, int EPoints, int number); 

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
void fill_orbit_borders (rmg_double_t * sg, int dimx, int dimy, int dimz);
void get_distributed_mat (double *bigmat, double *dismat);
void whole_to_tri_p (rmg_double_t * A_tri, rmg_double_t * Aii, int N, int *ni);
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
void diff_hx_interpolation (int st, double *xi, double *xi_old, int NX,
        double hx, double hx_old, double x0, double x0_old);
void find_phase (int nldim, rmg_double_t * nlcdrs, rmg_double_t * phase_sin, rmg_double_t * phase_cos);
void read_orbital (STATE * states);

void read_potrho_LCR (double *vh, double *vxc, double *rho);
void write_rho_x (rmg_double_t * rho, char *ab);
void interpolation_orbit (STATE * states);
void tri_to_whole_p (rmg_double_t * A_tri, rmg_double_t * Aii, int N, int *ni);
void tri_to_whole_complex_p (complex double * A_tri, complex double * Aii, int N, int *ni);
void global_to_distribute (rmg_double_t * global_array, rmg_double_t * distr_array);
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
void confine (rmg_double_t * mat, int size_x, int size_y, int size_z, COMPASS compass, int grid_type);
void distribute_to_global (rmg_double_t * distr_array, rmg_double_t * global_array);
void distribute_to_global2 (double *distr_array, double *global_array);
void distribute_to_X_soft (rmg_double_t * distr_array, rmg_double_t * global_array);
void distribute_to_Y_soft (rmg_double_t * distr_array, rmg_double_t * global_array);
void X_to_distribute_soft (rmg_double_t * global_array, rmg_double_t * distr_array);
void Y_to_distribute_soft (rmg_double_t * global_array, rmg_double_t * distr_array);
void Sgreen_c_p (rmg_double_t * Htri, rmg_double_t * Stri, complex double * sigma_L, int *sigma_idx,
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


void global_sums_X (rmg_double_t * vect, int *length);
void corlyp (rmg_double_t *dp, rmg_double_t *dm, rmg_double_t *dp1, rmg_double_t *dm1, rmg_double_t *dp2,
        rmg_double_t *dm2, rmg_double_t *ec, rmg_double_t *vcp0, rmg_double_t *vcm0, int *ndm);

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


void filter_potential (rmg_double_t *potential, rmg_double_t *r, int rg_points, rmg_double_t rmax, rmg_double_t offset, rmg_double_t parm, rmg_double_t* potential_lgrid, 
	rmg_double_t *rab, int l_value, rmg_double_t dr, rmg_double_t  gwidth, int lgrid_points, rmg_double_t rcut, rmg_double_t rwidth, rmg_double_t * drpotential_lgrid);

int min_distance_index(double *distance, int n );


void get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);

void mgrid_solv_negf(rmg_double_t * v_mat, rmg_double_t * f_mat, rmg_double_t * work,
                int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, rmg_double_t step, rmg_double_t k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);
void get_vh_negf (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);

double find_new_energy_point(double *cond, double *ener1, int tot_energy_point, double simpson_tol, int *EP_final, 
                           int *energy_insert_index, double *ener1_temp);

void zgetrf( int *n1, int *, complex double *Hii, int *, int *piv, int *info);
void zgetrs(char *trans, int *n1, int *, complex double *Hii, int *,int *ipiv, complex double *Gii, int *, int *info);
void zgetrf_acc( int *n1, int *, complex double *Hii, int *, int *piv, int *info);
void zgetrs_acc(char *trans, int *n1, int *, complex double *Hii, int *,int *ipiv, complex double *Gii, int *, int *info);
void zcopy_acc(int *n1, complex double *Hii, int *, complex double *Gii, int *);



rmg_double_t dot_product_orbit_nl(STATE *st1, int ion2, rmg_double_t *psi, rmg_double_t *prjptr);

void non_zero_pairs();
void non_zero_pairs1();

void init_nl_xyz(void);

void theta_phi_new(int st1, int st2, rmg_double_t theta_ion, rmg_double_t *st2_psi,
		                rmg_double_t *state1_psi, int mode, STATE *states); 

void    print_status( STATE *, rmg_double_t *, rmg_double_t *, rmg_double_t *, rmg_double_t *, char *);
void    print_global_function(rmg_double_t *, char, char *);
void    print_state_sum(STATE *states);
void    print_state(STATE *state);
void    print_sum(int size, double *data, char *msg); 
void    print_sum_square(int size, double *data, char *msg); 


void init_comm(STATE *states); 

void init_comm_res(STATE *states); 


void get_orbit_overlap_region(STATE *states); 

void get_ion_orbit_overlap_nl(STATE *states); 

void    duplicate_states_info(STATE *states, STATE *states1); 


void is_state_overlap(STATE *states); 


void set_energy_weight(complex double *ene, complex double *weight, rmg_double_t EF, int *nenergy);
void set_energy_weight_ne(complex double *eneR, complex double *weight, rmg_double_t EF1, rmg_double_t EF2, int *nenergy);


/* do we need these declarations? I think not... 
 * */
void zgemm(char*, char*, int*, int*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, complex double*, int*); 
/* void zgemm(char*, char*, int*, int*, int*, rmg_double_t*, rmg_double_t*, int*, rmg_double_t*, int*, rmg_double_t*, rmg_double_t*, int*);
 */
void ZCOPY(int*, complex double*, int*, complex double*, int*);
void ZAXPY(int*, complex double*, complex double*, int*, complex double*, int*);
void pzgemm(_fcd, _fcd, int*, int*, int*, complex double*,
        complex double*, int*, int*, int*,
        complex double*, int*, int*, int*,
        complex double*, complex double*, int*, int*, int*);



void Stransfer_f(rmg_double_t*, rmg_double_t*, rmg_double_t*, rmg_double_t*, int*);
void mat_part(rmg_double_t*, rmg_double_t*, rmg_double_t*, int);
void part_to_distribute(rmg_double_t*, rmg_double_t*, int);
void read_LCR();
void read_data_LCR(char *nameL, char *nameC, char *nameR, double *vh,
        double *vxc, double *rho, STATE *states);
void read_data_part(char *name, double *vh, double *vxc, double *rho, int which_part);

void Sigma(rmg_double_t *sigma, rmg_double_t *HLC, rmg_double_t *SLC, rmg_double_t eneR, rmg_double_t eneI, rmg_double_t *Green, int iprobe);

void Sgreen_c_wang(rmg_double_t *Htri, rmg_double_t *Stri, complex double *sigma_all, int *sigma_idx, rmg_double_t eneR, rmg_double_t eneI, complex double *Green_C, int nC);

void Sgreen_c(rmg_double_t *H00, rmg_double_t *S00, complex double *sigma, complex
        double *, complex double ene, complex double *Green_C, int nC);
void distri_fermi(complex double ene, rmg_double_t EF, complex double *distri);
void rho_munu(complex double *rho_mn, complex double *green_C, complex double *sigma_L, int nL, int N, int *ni, int iprobe);
void modify_rho(rmg_double_t *rho, rmg_double_t *rho_old);
void modify_rho_y(rmg_double_t *rho, rmg_double_t *rho_old);
void read_data_lead(double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void write_data_lead(char *name, double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void get_cond_lead();
void whole_to_tri_real(rmg_double_t *A_tri, rmg_double_t *Aii, int N, int *ni);

void    get_inverse_block(complex double *Hii, complex double *Gii, int *ipiv, int nn);

void read_data_conductor(char *name, double *vh, double *vxc, double *rho);


//void PZGESV(int*, int*, complex double *, int*, int*, int*, int*,
//                complex double*, int*, int*, int*, int*);


void genvpsi_gpu(const double *psi, const double *veff, double *psiv, 
        const int num_local_orbital, const int p0basis);
