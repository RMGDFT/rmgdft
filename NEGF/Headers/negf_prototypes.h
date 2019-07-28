#ifndef NEGF_PROTOTYPES_H_INCLUDED
#define NEGF_PROTOTYPES_H_INCLUDED

#include <complex>
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "twoParts.h"
#include "LCR.h"
#include "blacs.h"

void H01_to_HCL (double *H01_global, double *HCL_local, int iprobe);
void Sgreen_c (double * Htri, double * Stri, std::complex<double> * sigma1, std::complex<double> * sigma2, std::complex<double> ene, std::complex<double> * Green_C, int nC);
void Sgreen_c_noneq_p (double *Htri, double *Stri, std::complex<double> * sigma, int *sigma_idx, std::complex<double> ene, std::complex<double> *Green_C, std::complex<double> *Green_C_row, std::complex<double> *Green_C_col, int nC, int iprobe);
void Sgreen_c_p (double * Htri, double * Stri, std::complex<double> * sigma, int * sigma_idx, std::complex<double> ene, std::complex<double> * Green_C);
void Sgreen_cond_p (std::complex<double> *H_tri, std::complex<double> *G_tri, std::complex<double> *sigma_all, int *sigma_idx, std::complex<double> *green_C, int nC, int iprobe1, int iprobe2);
void X_to_distribute_soft (double * global_array, double * distr_array);
void allocate_matrix_LCR (void);
void allocate_matrix_soft (void);
void apply_potential_drop (double *vbias);
void charge_density_matrix_p (std::complex<double> * sigma_all);
void diaginit (double *aa, int *desca, double *a, int lda);
void diff_hx_interpolation (int st, double *xi, double *xi_old, int NX, double hx, double hx_old, double x0, double x0_old);
void distri_fermi (std::complex<double> ene, double EF, std::complex<double> *distri);
void distribute_to_X_soft (double * distr_array, double * global_array);
void dzasum_driver(int n, std::complex<double> *A, int ia, double *sum);
void finalize_gpu (void);
void find_fermi (std::complex<double> * sigma_all);
double find_new_energy_point(double *cond, double *ener1, int tot_energy_point, double simpson_tol, int *EP_final, int *energy_insert_index, double *ener1_temp);
void gauleg (double x1, double x2, double *x, double *w, int n);
void get_3Ddos (STATE * states, double EMIN, double EMAX, int EPoints, int number);
void get_cond_frommatrix (void);
void get_cond_frommatrix_kyz (void);
void get_ddd_update (double * veff);
void get_distributed_mat (double *bigmat, double *dismat);
void get_dos (STATE * states);
void get_inverse_block (std::complex<double> *Hii, std::complex<double> *Gii, int *ipiv, int nn);
void get_inverse_block_p (std::complex<double> *Hii, std::complex<double> *Gii, int *ipiv, int *desca );
void get_ion_orbit_overlap_loc (STATE * states);
char *get_line (char *buf, FILE * fh);
char *get_num (char *str);
void get_state_to_proc (STATE * states);
void getvector_device_host (int n, int elemsize, void *x_device, int ia, void *x_host, int ib);
void global_to_distribute3 (double *global_array, double *distr_array);
void green_kpoint_phase (std::complex<double> *green, double kvecy, double kvecz, int up_and_low);
void green_lead (std::complex<double> *ch0_host, std::complex<double> *ch01_host, std::complex<double> *ch10_host, std::complex<double> *green_host, int iprobe);
void init_comp (double *vh);
void init_ext (double *vext, double gbias_begin, double gbias_end, double BT, double gate_bias);
void init_gpu (void);
void init_loc_xyz (void);
void init_pe_on (void);
void interpolation_orbit (STATE * states);
void ion_partial_Hij_and_Sij (int iion, int flag, double *Hij, double *Sij);
void is_vloc_state_overlap (STATE *states);
void kpoints(int *nkp, double *kvecx, double *kvecy, double *kvecz, int *nkp_tot1, double *kweight);
void lead_bandstructure (void);
void lead_mat_distribute (double *a_local, int *desca, double *a_global,int iprobe);
void matgather (double *aa, int *desca, double *a, int lda);
void matrix_inverse_Gauss (std::complex<double> * H_tri_host, std::complex<double> * G_tri_host);
void matrix_inverse_anyprobe (std::complex<double> * H_tri_host, int N, int * ni, int iprobe,
        std::complex<double> * G_row_host, std::complex<double> *G_col_host);
void matrix_inverse_blocknm (std::complex<double> * H_tri_host, int N_blocks, int * ni, int m, int n, std::complex<double> * Green_C_host);
void matrix_inverse_blocknm_Gauss (std::complex<double> * H_tri_host, std::complex<double> *G_tri_host, int m, int n, std::complex<double> * Green_C_host);
void matrix_inverse_driver (std::complex<double> *Hii, int *desca );
void matrix_inverse_p (std::complex<double> * H_tri_host, std::complex<double> * G_tri_host);
void matrix_inverse_rowcol (std::complex<double> * H_tri_host, int iprobe, std::complex<double> *G_tri_host, std::complex<double> *G_row_host, std::complex<double> *G_col_host);
void matrix_kpoint (int size, std::complex<double> *matrix_k, double *matrix_yz, double kvecy, double kvecz);
double *memory_ptr_host_device(double *ptr_host, double *ptr_device);
std::complex<double> *memory_ptr_host_device(std::complex<double> *ptr_host, std::complex<double> *ptr_device);
int min_distance_index(double *distance, int n );
void modify_rho (double * rho, double * rho_old);
void multi_GHG_munu (double *GHG_tri, double *GHG_en_tri);
double my_crtc (void);
void nlforce (double * veff);
void nlforce_par_D (STATE *states, double *forces);
void nlforce_par_Q(double * veff, double * rho_nm, int ion, int nh, double *forces);
void nlforce_partial_H_part2 (STATE * states, STATE * states1, double *GHG, double *force);
void ortho_norm_local (STATE *states);
void partial_vloc (void);
void partial_vlocpsi (STATE st1, int ion2, double * psi, double * prjptr, double *vlpsi);
void plane_average_rho (double *rho);
void pmo_init (void);
double pmo_trace(std::complex<double> *matrix, int *desca);
void pmo_unitary_matrix(std::complex<double> *a_local, int *desca);
void pmo_unitary_matrix_double(double *a_local, int *desca);
void print_data (int size, double *data);
void print_data_to_file (int size, double *data, char *filename);
void print_state (STATE * state);
void print_state_sum (STATE * states);
void print_states_dot_product (STATE * states);
void print_status (STATE * states, double * vh, double * vxc, double * vnuc, double * vh_old, char *msg);
void print_sum (int size, double *data, char *msg);
void print_sum_idx (int size, double *data, char *msg);
void print_sum_square (int size, double *data, char *msg);
void read_LCR (void);
void read_cond_input (double *emin, double *emax, int *E_POINTS, double *E_imag, double *KT, int *kpoint);
void read_lead_matrix (void);
void read_matrix_pp (void);
void read_orbital (STATE * states);
void read_potrho (double *vh, int iflag, char *file_ex);
void read_potrho_LCR (double *vh, double *vxc, double *rho);
void read_rho_and_pot (char *name, double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void read_trans (complex_energy_integral * cei);
void rho_munu (std::complex<double> * rho_mn_host, std::complex<double> * G_row_host, std::complex<double> *G_col_host, std::complex<double> * gamma_host, int iprobe);
void rho_nm_mat (double *Aij, double * global_mat_X);
void rmg_printout_devices( );
void row_to_tri_p (double * A_tri, double * Aii, int N, int *ni);
void scale_orbital (STATE * states);
void set_energy_weight (std::complex<double> * ene, std::complex<double> * weight, double EF, int *nenergy);
void set_energy_weight_ne (std::complex<double> * ene, std::complex<double> * weight, double EF1, double EF2, int *nenergy);
void setback_corner_matrix_H(void);
void setback_corner_matrix_S(void);
void setvector_host_device (int n, int elemsize, void *x_host, int ia, void *x_device, int ib);
void sigma_all_energy_point (std::complex<double> * sigma_all, double kvecy, double kvecz);
void sigma_one_energy_point (std::complex<double> *sigma, int jprobe, std::complex<double> ene, double kvecy, double kvecz, std::complex<double> *work);
void sl_init_on (int *ictxt, int nprow, int npcol);
void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y);
void split_matrix_center (void);
void split_matrix_lead (int iprobe);
void state_minus_state (STATE *orbit1, STATE *orbit2, double factor);
double trace_AB_tri(double * H_tri, double * par_D_tri, int N, int *ni);
void tri_to_row (double * A_tri, double * Aii_row, int N, int *ni);
void tri_to_whole_complex (std::complex<double> * H_tri, std::complex<double> * Hii, int N, int * ni);
void tri_to_whole_complex_p (std::complex<double> * A_tri, std::complex<double> * Aii, int N, int *ni);
void whole_to_tri_complex (std::complex<double> * H_tri, std::complex<double> * Hii, int N, int * ni);
void whole_to_tri_real (double * A_tri, double * Aii, int N, int *ni);
void whole_to_tri_update (double * A_tri, double * Aii, int N, int *ni);
void write_data_lead (char *name, double *vh, double *vxc, double *vh_old, double *vxc_old, double *rho);
void write_global_data (int file_handle, double *data, int fnx, int fny, int fnz);
void write_global_data_lead (int file_handle, double *data, int fnx, int fny, int fnz);
void writeout_matrix_p (void);
void zaxpy_driver (int n, std::complex<double> alpha, std::complex<double> *A, int ia, std::complex<double> *B, int ib);
void zcopy_driver (int n, std::complex<double> *A, int ia, std::complex<double> *B, int ib);
void zero_lead_image(double *tri);
void zgemm_driver (char *transa, char *transb, int m, int n, int k, std::complex<double> alpha, std::complex<double> *A, int ia, int ja, int *desca,std::complex<double> *B, int ib, int jb, int *descb, std::complex<double> beta, std::complex<double> *C, int ic, int jc, int *descc);
void zgesv_driver (std::complex<double> *A, int *desca, std::complex<double> *B, int *descb );

void Sgreen_semi_infinite_p (std::complex<double> * green_host, std::complex<double> *ch00_host,
     std::complex<double> *ch01_host, std::complex<double> *ch10_host, int jprobe);

void Sigma_p (std::complex<double> *sigma, std::complex<double> *ch, std::complex<double> *ch01,
     std::complex<double> *ch10, std::complex<double> *green, int iprobe);

void find_new_energy_point_sharp_peak(double *cond, double *ener1, int tot_energy_point,
     double critical_val, int *EP_final, int *energy_insert_index, double *ener1_temp);

void matrix_kpoint_center (std::complex<double> *H_tri, double *Stri, double *Htri,
     std::complex<double> ene, double kvecy, double kvecz);

void matrix_kpoint_lead (std::complex<double> *S00, std::complex<double> *H00, std::complex<double> *S01,
     std::complex<double> *H01, std::complex<double> *SCL, std::complex<double> *HCL, double kvecy, double kvecz, int iprobe);

void confine (double * mat, int size_x, int size_y, int size_z, COMPASS compass, int level);

void get_vh_negf (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target);

void mgrid_solv_negf(double * v_mat, double * f_mat, double * work,
                int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz,
                int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);

void local_current (void);

void Sgreen_onerow (std::complex<double> *Htri, std::complex<double> * sigma,
                     int *sigma_idx, std::complex<double> * Green_C, int nC,
                     int iprobe);

void comm_sums (double * vect, int *length, MPI_Comm COMM_TEM);
int int_max_all(int);
void global_sums_X (double * vect, int *length);
void global_sums_int (int *, int *);

#endif
