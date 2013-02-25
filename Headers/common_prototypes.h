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


