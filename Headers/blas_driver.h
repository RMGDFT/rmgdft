#ifndef RMG_blas_driver_h
#define RMG_blas_driver_h 1


void zcopy_driver (int n, std::complex<double> *A, int ia, std::complex<double> *B, int ib) ;
void zaxpy_driver (int n, std::complex<double> alpha, std::complex<double> *A, int ia, std::complex<double> *B, int ib) ;
void dzasum_driver(int n, std::complex<double> *A, int ia, double *sum);
void dcopy_driver (int n, double *A, int ia, double *B, int ib) ;
void daxpy_driver (int n, double alpha, double *A, int ia, double *B, int ib) ;
void dscal_driver(int n, double beta, double *A, int ione);
void zgemm_driver (char *transa, char *transb, int m, int n, int k,
std::complex<double> alpha, std::complex<double> *A, int ia, int ja, int *desca,
std::complex<double> *B, int ib, int jb, int *descb, std::complex<double> beta,
std::complex<double> *C, int ic, int jc, int *descc);
void dgemm_driver (char *transa, char *transb, int m, int n, int k,
double alpha, double *A, int ia, int ja, int *desca,
double *B, int ib, int jb, int *descb, double beta,
double *C, int ic, int jc, int *descc);
void my_sync_device();

#endif
