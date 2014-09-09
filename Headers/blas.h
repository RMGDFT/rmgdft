/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/blas.h *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   
 * INPUTS
 *
 * OUTPUT
 *  
 * PARENTS
 *
 * CHILDREN
 * 
 * SEE ALSO
 *  
 * SOURCE
 */
/*

			blas.h


    Lets us call blas routines from c. The Crays require that the function
    call be in Capital letters while the RS6000's require lowercase. Also
    we use the double precision form on the IBM's and the single precision
    form on the crays. (Even though a double in C is actually a fortran 
    single.)

*/

#ifdef LINUX 

#define                 sger            dger_
//#  define		global_sums     global_sums_
#  define		saxpy		daxpy_
#  define		daxpy		daxpy_
#  define		zaxpy		zaxpy_
#  define		ZAXPY		zaxpy_
#  define		scopy   	dcopy_
#  define		dcopy   	dcopy_
#  define		zcopy   	zcopy_
#  define		ZCOPY   	zcopy_
#  define		sdot    	ddot_
#  define		ddot    	ddot_
#  define		sscal   	dscal_
#  define		dscal   	dscal_
#  define		ssyrk   	dsyrk_
#  define		spotrf  	dpotrf_
#  define		dpotrf  	dpotrf_
#  define		spotri  	dpotri_
#  define		cpotrf  	zpotrf_
#  define		cpotri  	zpotri_
#  define		snrm2   	dnrm2_
#  define		sswap   	dswap_
#  define		dswap   	dswap_
#  define		ssyev		dsyev_
#  define		cheev		zheev_
#  define		xbecke		xbecke_
#  define		corlyp_f90      corlyp_f90_
#  define		corpbe  	corpbe_
#  define		exchpbe 	exchpbe_
#  define		exch		exch_
#  define		corlsd		corlsd_
#  define		corgga		corgga_
#  define		symrho  	symrho_
#  define		fsymforces  	fsymforces_
#  define		symmetry 	symmetry_
#  define		fgram   	dfgram_
#  define		dlamch  	dlamch_
#  define		dgemm  	        dgemm_
#  define		zgemm  	        zgemm_
#  define		ZGEMM  	        zgemm_
#  define		dgetrf 	        dgetrf_
#  define		dgetri 	        dgetri_
#  define		zgetrf 	        zgetrf_
#  define		zgetri 	        zgetri_
#  define		dgesv 	        dgesv_
#  define               sgesv           dgesv_
#  define		zgesv 	        zgesv_
#  define		ZGESV 	        zgesv_
#  define               dgemv           dgemv_
#  define               dsygvx          dsygvx_
#  define               zhegvx          zhegvx_
#  define               zhegst          zhegst_
#  define               zheevd          zheevd_
#  define               dsyev           dsyev_
#  define               dsyevd          dsyevd_
#  define               dsyevx          dsyevx_
#  define               dsyevr          dsyevr_
#  define               zheev           zheev_


#else
#ifdef AIX

#  define		saxpy		daxpy
#  define		scopy   	dcopy
#  define		sdot    	ddot
#  define		sscal   	dscal
#  define		ssyrk   	dsyrk
#  define		spotrf  	dpotrf
#  define		spotri  	dpotri
#  define		snrm2   	dnrm2
#  define		sswap   	dswap
#  define		ssyev		dsyev
#  define		cpotrf  	zpotrf
#  define		cpotri  	zpotri
#  define		cheev   	zheev
#  define		fgram   	dfgram

#endif
#endif


#if __cplusplus
extern "C" {
#endif

void my_copy(double *in, double *out, int length);
void my_scal(double alpha, double *vect, int length);
void my_axpy(double alpha, double *in, double *out, int length);
void my_swap(double *vec1, double *vec2, int length);


void fsymforces (double * force, int *s, int *irg, int *irt,
                 int *nat, int *ibrav, int *nsym,
                 double * celldm, int *nr1, int *nr2, int *nr3);
int ilaenv (int *ispec, char *name, char *opts, int *n1, int *n2, int *n3,
            int *n4);
void saxpy (int *n, double * alpha, double * x, int *incx, double * y, int *incy);
void sscal (int *n, double * alpha, double * x, int *incx);
void scopy (int *n, double * x, int *incx, double * y, int *incy);
double sdot (int *n, double * x, int *incx, double * y, int *incy);
double snrm2 (int *n, double * x, int *incx);
void ssyrk (const char *uplo, const char *trans, int *n, int *k, double * alpha, double * A,
            int *lda, double * beta, double * c, int *ldc);
void sswap (int *n, double * x, int *incx, double * y, int *incy);
void xbecke (double * d, double * s, double * u, double * v, double * ex, double * vx);
void corlyp_f90 (double * dp, double * dm, double * dp1, double * dm1, double * dp2,
                 double * dm2, double * ec, double * vcp0, double * vcm0, int *ndm);

void corlsd (double * rs, double * zet, double * ec, double * vcup, double * vcdn,
             double * ecrs, double * eczet, double * alfc);
void corgga (double * rs, double * zet, double * t, double * uu, double * vv,
             double * ww, double * h, double * dvcup, double * dvcdn,
             double * fk, double * sk, double * g, double * ec, double * ecrs,
             double * eczet);
void corpbe (double * rs, double * zet, double * t, double * uu, double * vv, double * ww,
             int *lgga, int *lpot, double * ec, double * vcup, double * vcdn,
             double * h, double * dvcup, double * DVCDN);
void cpotrf (char *uplo, int *n, double * a, int *lda, int *info);
void exch (double * d, double * s, double * u, double * v, double * ex, double * vx);
void exchpbe (double * d, double * s, double * u, double * v, int *lgga, int *lpot,
              double * ex, double * vx);
void symrho (double * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
             int *irg, int *ftau);
void ssyev (char *jobz, char *uplo, int *numm, double * ss, int *numn,
            double * work1, double * work2, int *lwork, int *info);
void cheev (char *jobz, char *uplo, int *numm, double * ss, int *numn,
            double * work1, double * work2, int *lwork, double *, int *info);
void spotri (char *, int *, double *, int *, int *);
void dgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void zgemm(const char *, const char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void dgetrf( int *, int *, double *, int *, int *, int *);
void dgetri(int *, double *, int *, int *, double *, int *, int *);
void zgetrf( int *, int *, double *, int *, int *, int *);
void zgetri(int *, double *, int *, int *, double *, int *, int *);
void zhegst(int *, const char *, int *, double *, int *, double *, int *, int *);
void zheevd(const char *, const char *, int *, double *, int *, double *, double *, int *, double *, int *, int *, int *, int *);
void dgesv (int *, int*, double *, int *, int *, double *, int *, int *);
void zgesv (int *, int*, double *, int *, int *, double *, int *, int *);
void dgemv ( char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
void dsygvx (int *itype, char *jobz, char *range, char *uplo, int *n, double *A, int *lda, double *B, int *ldb, double *vl, double *vu, 
             int *il, int *iu, double *tol, int *m, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *ifail, int *info);
void zhegvx (int *itype, char *jobz, char *range, char *uplo, int *n, double *A, int *lda, double *B, int *ldb, double *vl, double *vu, 
             int *il, int *iu, double *tol, int *m, double *w, double *z, int *ldz, double *work, int *lwork, double *rwork, int *iwork, int *ifail, int *info);
void dsyev (char *jobz, char *uplo, int *numm, double * ss, int *numn,
            double * work1, double * work2, int *lwork, int *info);
void dsyevd (char *jobz, char *uplo, int *numm, double * ss, int *numn,
            double * work1, double * work2, int *lwork, int *iwork, int *liwork, int *info);
void zheev (char *jobz, char *uplo, int *numm, double * ss, int *numn,
            double * work1, double * work2, int *lwork, double *rwork, int *info);
void dsyevx (char *, char *, char *, int *, double *, int *, double *, double *, int *, int *, double *, int *, double *, double *, int *, double *, int *, int *, int *, int *);
void dsyevr (char *, char *, char *, int *, double *, int *, double *, double *, int *, int *, double *, int *, double *, double *, int *, int *, double *, int *, int *, int *, int *);




#if __cplusplus
}
#endif
/******/
