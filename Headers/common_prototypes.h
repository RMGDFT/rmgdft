/*
  Function prototypes shared by all codes.


*/

#ifndef RMG_COMMON_PROTOTYPES_H
#define RMG_COMMON_PROTOTYPES_H 1

#include "mpi.h"

#if CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#endif

#include "rmgtypedefs.h"

template <typename T> T MyConj(T val);

void WritePoscar(FILE *FH, int md_step);

/* Blas wrappers */
void QMD_daxpy (int n, double alpha, double *x, int incx, double *y, int incy);
void QMD_dscal (int n, double alpha, double *x, int incx);
void QMD_dcopy (int n, double *x, int incx, double *y, int incy);
double QMD_ddot (int n, double *x, int incx, double *y, int incy);


int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype);
void mgrid_solv (double *v_mat, double *f_mat, double *work,
                 int dimx, int dimy, int dimz, double gridhx, double gridhy,
                 double gridhz, int level, int max_levels,
                 int *pre_cyc, int *post_cyc, int mu_cyc, double step, double k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundaryflag);
void eval_residual (double *mat, double *f_mat, double *work, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, double *res);
void solv_pois (double *vmat, double *fmat, double *work,
                int dimx, int dimy, int dimz, double gridhx,
                double gridhy, double gridhz, double step, double k);


void mg_restrict (double *full, double *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
void mg_prolong (double *full, double *half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);
double app2_del2 (double *a, double *b, int dimx, int dimy, int dimz,
                double gridhx, double gridhy, double gridhz);
void init_global_sums(void);

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

void find_node_offsets(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxoffset, int *pyoffset, int *pzoffset);
void find_node_sizes(int gridpe, int nxgrid, int nygrid, int nzgrid,
                      int *pxsize, int *pysize, int *pzsize);
void pe2xyz(int pe, int *x, int *y, int *z);
int radius2grid (double radius, double mingrid_spacing);


int get_PE_X(void);
int get_PE_Y(void);
int get_PE_Z(void);
int get_NX_GRID(void);
int get_NY_GRID(void);
int get_NZ_GRID(void);
int get_FNX_GRID(void);
int get_FNY_GRID(void);
int get_FNZ_GRID(void);
int get_FG_RATIO(void);
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
void set_anisotropy(double a);
double get_anisotropy(void);
void set_neighbors(int *list);
int *get_neighbors(void);
void set_grids(int NX_GRID, int NY_GRID, int NZ_GRID, int PE_X, int PE_Y, int PE_Z, int FG_RATIO);
void set_rank(int newgridpe, MPI_Comm comm);
int get_ibrav_type(void);
void set_ibrav_type(int ibrav);
double get_xside(void);
double get_yside(void);
double get_zside(void);
double get_hxgrid(void);
double get_hygrid(void);
double get_hzgrid(void);
double get_hxxgrid(void);
double get_hyygrid(void);
double get_hzzgrid(void);
double get_vel(void);
double get_vel_f(void);
double get_celldm(int which);
double get_a0(int which);
double get_a1(int which);
double get_a2(int which);
double get_b0(int which);
double get_b1(int which);
double get_b2(int which);
double get_omega(void);


void trade_images (double *mat, int dimx, int dimy, int dimz, int type);
void trade_imagesx (double *f, double *w, int dimx, int dimy, int dimz, int images, int type);
double AtomicInterpolate(double *f, double r);
void lbfgs (double *posion, double *force, int num_coeff);
void Vdd(double * rho);
void WriteChargeAnalysis(void);
ssize_t rmg_write(int fd, const void *buf, size_t count);


#if CUDA_ENABLED
void init_gpu (void);
void finalize_gpu (void);
#endif



#endif
