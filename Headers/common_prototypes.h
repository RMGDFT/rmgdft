/*
  Function prototypes shared by all codes.


*/

#ifndef RMG_COMMON_PROTOTYPES_H
#define RMG_COMMON_PROTOTYPES_H 1

#include "make_conf.h"

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#endif

#include "rmgtypedefs.h"


/* Blas wrappers */
void QMD_daxpy (int n, rmg_double_t alpha, rmg_double_t *x, int incx, rmg_double_t *y, int incy);
void QMD_dscal (int n, rmg_double_t alpha, rmg_double_t *x, int incx);
void QMD_dcopy (int n, rmg_double_t *x, int incx, rmg_double_t *y, int incy);
rmg_double_t QMD_ddot (int n, rmg_double_t *x, int incx, rmg_double_t *y, int incy);
void QMD_saxpy (int n, rmg_float_t alpha, rmg_float_t *x, int incx, rmg_float_t *y, int incy);
void QMD_sscal (int n, rmg_float_t alpha, rmg_float_t *x, int incx);
void QMD_scopy (int n, rmg_float_t *x, int incx, rmg_float_t *y, int incy);
rmg_float_t QMD_sdot (int n, rmg_float_t *x, int incx, rmg_float_t *y, int incy);


int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype);
void mgrid_solv (rmg_double_t *v_mat, rmg_double_t *f_mat, rmg_double_t *work,
                 int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy,
                 rmg_double_t gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, rmg_double_t step, rmg_double_t k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);
void mgrid_solv_f (rmg_float_t *v_mat, rmg_float_t *f_mat, rmg_float_t *work,
                 int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy,
                 rmg_double_t gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, rmg_double_t step, rmg_double_t k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim);
void eval_residual (rmg_double_t *mat, rmg_double_t *f_mat, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t *res);
void eval_residual_f (rmg_float_t *mat, rmg_float_t *f_mat, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_float_t *res);
void solv_pois (rmg_double_t *vmat, rmg_double_t *fmat, rmg_double_t *work,
                int dimx, int dimy, int dimz, rmg_double_t gridhx,
                rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t step, rmg_double_t k);
void solv_pois_f (rmg_float_t *vmat, rmg_float_t *fmat, rmg_float_t *work,
                int dimx, int dimy, int dimz, rmg_double_t gridhx,
                rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t step, rmg_double_t k);


void mg_restrict (rmg_double_t *full, rmg_double_t *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_restrict_f (rmg_float_t *full, rmg_float_t *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_prolong (rmg_double_t *full, rmg_double_t *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_prolong_f (rmg_float_t *full, rmg_float_t *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

rmg_double_t app_del2c (rmg_double_t *a, rmg_double_t *b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_del2c_f (rmg_float_t *a, rmg_float_t *b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void init_global_sums(void);

rmg_double_t app_cilr_sixth (rmg_double_t * psi, rmg_double_t *a_psi, rmg_double_t *b_psi, rmg_double_t *vtot_eig_s, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t app_cilr_fourth (rmg_double_t * psi, rmg_double_t *a_psi, rmg_double_t *b_psi, rmg_double_t *vtot_eig_s, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);

void app6_del2 (rmg_double_t *rho, rmg_double_t *work, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
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

void trade_imagesx_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images, int type);


void app_cir_beta_fourth (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz);
void app_cir_beta_sixth (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz);

int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset);
int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxsize, int *pysize, int *pzsize);
void read_common();
void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_f (rmg_float_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images1_central_async_f (rmg_float_t * f, int dimx, int dimy, int dimz);
void trade_images1_central_async (rmg_double_t * f, int dimx, int dimy, int dimz);
void trade_imagesx_central_async_f (rmg_float_t * f, rmg_float_t * w, int dimx, int dimy, int dimz, int images);
void trade_imagesx_central_async (rmg_double_t * f, rmg_double_t * w, int dimx, int dimy, int dimz, int images);
void pe2xyz(int pe, int *x, int *y, int *z);
int radius2grid (rmg_double_t radius, rmg_double_t mingrid_spacing);
int find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset);
int find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxsize, int *pysize, int *pzsize);


int get_PX0_GRID(void);
int get_PY0_GRID(void);
int get_PZ0_GRID(void);
int get_PX_OFFSET(void);
int get_PY_OFFSET(void);
int get_PZ_OFFSET(void);
int get_FPX_OFFSET(void);
int get_FPY_OFFSET(void);
int get_FPZ_OFFSET(void);
int get_P0_BASIS(void);
int get_FP0_BASIS(void);
int get_FPX0_GRID(void);
int get_FPY0_GRID(void);
int get_FPZ0_GRID(void);
void set_anisotropy(rmg_double_t a);
rmg_double_t get_anisotropy(void);
void set_neighbors(int *list);
int *get_neighbors(void);
void set_grids(int gridpe, int ii, int jj, int kk);
int get_ibrav_type(void);
void set_ibrav_type(int ibrav);
ION *get_ion(int ion);
rmg_double_t get_xside(void);
rmg_double_t get_yside(void);
rmg_double_t get_zside(void);
rmg_double_t get_hxgrid(void);
rmg_double_t get_hygrid(void);
rmg_double_t get_hzgrid(void);

void trade_images (rmg_double_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_images_f (rmg_float_t *mat, int dimx, int dimy, int dimz, int *nb_ids, int type);
void trade_imagesx (rmg_double_t *f, rmg_double_t *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images, int type);
void trade_imagesx_async (rmg_double_t *f, rmg_double_t *w, int dimx, int dimy, int dimz, int images);
void trade_imagesx_async_f (rmg_float_t *f, rmg_float_t *w, int dimx, int dimy, int dimz, int images);
void trade_images1_async (rmg_double_t * f, int dimx, int dimy, int dimz);
void trade_images1_async_f (rmg_float_t * f, int dimx, int dimy, int dimz);


// This set is used for the C++ interface
void FD_app_cir_sixth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz);
void FD_app_cir_sixth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz);
void FD_app_cir_sixth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b);
void FD_app_cir_sixth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b);
rmg_double_t FD_app_cil_sixth_standard_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t FD_app_cil_sixth_standard_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t FD_app_cil_sixth_global_rmg_double(rmg_double_t *rptr, rmg_double_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t FD_app_cil_sixth_global_rmg_float(rmg_float_t *rptr, rmg_float_t *b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t FD_app_del2c_rmg_double(rmg_double_t *rptr, rmg_double_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
rmg_double_t FD_app_del2c_rmg_float(rmg_float_t *rptr, rmg_float_t *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);




#if GPU_ENABLED
void init_gpu (void);
void finalize_gpu (void);
double app_cil_sixth_f_gpu(const float *psi,
                        float *b,
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
void app_cir_sixth_f_gpu(const float *psi,
                        float *b,
                        const int dimx,
                        const int dimy,
                        const int dimz,
                        cudaStream_t cstream);

double app_cil_sixth_gpu(const double *psi,
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
void app_cir_fourth_gpu(const double *psi,
                        double *b, 
                        const int dimx,
                        const int dimy,
                        const int dimz,
                        cudaStream_t cstream);
double app_cil_fourth_gpu(const double *psi,
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
double app_cil_fourth_f_gpu(const float *psi,
                       float *b,
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

#define my_fopen(_fhandle_, _filename_, _mode_) do {\
    _fhandle_ = fopen(_filename_, _mode_);\
    if (_fhandle_ == NULL)\
        error_handler("macro my_fopen: can't fopen file %s", _filename_);\
    printf("\nfopening file '%s' mode '%s'\n", _filename_, _mode_);\
} while (0)



#define my_open(_fhandle_, _filename_, _flags_, _mode_) do {\
    _fhandle_ = open(_filename_, _flags_, _mode_);\
    if (_fhandle_ < 0)\
        error_handler("macro my_open: can't open file %s", _filename_);\
    printf("\nopening file '%s' flags '%s' mode '%s'\n", _filename_, #_flags_, #_mode_);\
} while (0)

#endif

#endif
