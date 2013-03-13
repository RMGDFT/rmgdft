/*
  Function prototypes shared by all codes.


*/

/* Blas wrappers */
void QMD_daxpy (int n, REAL alpha, REAL *x, int incx, REAL *y, int incy);
void QMD_dscal (int n, REAL alpha, REAL *x, int incx);
void QMD_dcopy (int n, REAL *x, int incx, REAL *y, int incy);
REAL QMD_ddot (int n, REAL *x, int incx, REAL *y, int incy);
void QMD_saxpy (int n, rmg_float_t alpha, rmg_float_t *x, int incx, rmg_float_t *y, int incy);
void QMD_sscal (int n, rmg_float_t alpha, rmg_float_t *x, int incx);
void QMD_scopy (int n, rmg_float_t *x, int incx, rmg_float_t *y, int incy);
rmg_float_t QMD_sdot (int n, rmg_float_t *x, int incx, rmg_float_t *y, int incy);


int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype);
void mgrid_solv (REAL *v_mat, REAL *f_mat, REAL *work,
                 int dimx, int dimy, int dimz, REAL gridhx, REAL gridhy,
                 REAL gridhz, int level, int *nb_ids, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, REAL step, REAL k,
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
void eval_residual (REAL *mat, REAL *f_mat, int dimx, int dimy, int dimz,
                    REAL gridhx, REAL gridhy, REAL gridhz, REAL *res);
void eval_residual_f (rmg_float_t *mat, rmg_float_t *f_mat, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_float_t *res);
void solv_pois (REAL *vmat, REAL *fmat, REAL *work,
                int dimx, int dimy, int dimz, REAL gridhx,
                REAL gridhy, REAL gridhz, REAL step, REAL k);
void solv_pois_f (rmg_float_t *vmat, rmg_float_t *fmat, rmg_float_t *work,
                int dimx, int dimy, int dimz, rmg_double_t gridhx,
                rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t step, rmg_double_t k);


void mg_restrict (REAL *full, REAL *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_restrict_f (rmg_float_t *full, rmg_float_t *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_prolong (REAL *full, REAL *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_prolong_f (rmg_float_t *full, rmg_float_t *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

REAL app_del2c (REAL *a, REAL *b, int dimx, int dimy, int dimz,
                REAL gridhx, REAL gridhy, REAL gridhz);
rmg_double_t app_del2c_f (rmg_float_t *a, rmg_float_t *b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
void init_global_sums(void);

#if GPU_ENABLED
void init_gpu (void);
void finalize_gpu (void);
void app_cil_sixth_f_gpu(const float *psi,
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

void app_cil_sixth_gpu(const double *psi,
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

#endif

