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
#include "rmgtypedefs.h"
#include "rmg_mangling.h"


#define		sger		RMG_FC_GLOBAL(sger, SGER)
#define		dger		RMG_FC_GLOBAL(dger, DGER)
#define		daxpy		RMG_FC_GLOBAL(daxpy, DAXPY)
#define		zaxpy		RMG_FC_GLOBAL(zaxpy, ZAXPY)
#define		dzasum		RMG_FC_GLOBAL(dzasum, DZASUM)
#define		dcopy		RMG_FC_GLOBAL(dcopy, DCOPY)
#define		zcopy		RMG_FC_GLOBAL(zcopy, ZCOPY)
#define		ddot		RMG_FC_GLOBAL(ddot, DDOT)
#define		dscal		RMG_FC_GLOBAL(dscal, DSCAL)
#define		dsyrk		RMG_FC_GLOBAL(dsyrk, DSYRK)
#define		dsyr2k		RMG_FC_GLOBAL(dsyr2k, DSYR2K)
#define		zsyr2k		RMG_FC_GLOBAL(zsyr2k, ZSYR2K)
#define		dpotrf		RMG_FC_GLOBAL(dpotrf, DPOTRF)
#define		zpotrf		RMG_FC_GLOBAL(zpotrf, ZPOTRF)
#define		dswap		RMG_FC_GLOBAL(dswap, DSWAP)
#define		cheev		RMG_FC_GLOBAL(cheev, CHEEV)
#define		xbecke		RMG_FC_GLOBAL(xbecke, XBECKE)
#define		corlyp_f90	RMG_FC_GLOBAL(corlyp_f90, CORLYP_F90)
#define		corpbe		RMG_FC_GLOBAL(corpbe, CORPBE)
#define		exchpbe		RMG_FC_GLOBAL(exchpbe, EXCHPBE)
#define		exch		RMG_FC_GLOBAL(exch, EXCH)
#define		corlsd		RMG_FC_GLOBAL(corlsd, CORLSD)
#define		corgga		RMG_FC_GLOBAL(corgga, CORGGA)
#define		symrho		RMG_FC_GLOBAL(symrho, SYMRHO)
#define		symmetry	RMG_FC_GLOBAL(symmetry, SYMMETRY)
#define		fgram		RMG_FC_GLOBAL(fgram, FGRAM)
#define		dlamch		RMG_FC_GLOBAL(dlamch, DLAMCH)
#define		dgemm		RMG_FC_GLOBAL(dgemm, DGEMM)
#define		zgemm		RMG_FC_GLOBAL(zgemm, ZGEMM)
#define		ZGEMM		RMG_FC_GLOBAL(ZGEMM, ZGEMM)
#define		dgetrf		RMG_FC_GLOBAL(dgetrf, DGETRF)
#define		dgetri		RMG_FC_GLOBAL(dgetri, DGETRI)
#define		zgetrf		RMG_FC_GLOBAL(zgetrf, ZGETRF)
#define		zgetri		RMG_FC_GLOBAL(zgetri, ZGETRI)
#define		dpotri		RMG_FC_GLOBAL(dpotri, DPOTRI)
#define		dgesv		RMG_FC_GLOBAL(dgesv, DGESV)
#define		zgesv		RMG_FC_GLOBAL(zgesv, ZGESV)
#define		dgemv		RMG_FC_GLOBAL(dgemv, DGEMV)
#define		dsygvx		RMG_FC_GLOBAL(dsygvx, DSYGVX)
#define		dsygvd		RMG_FC_GLOBAL(dsygvd, DSYGVD)
#define		zhegvx		RMG_FC_GLOBAL(zhegvx, ZHEGVX)
#define		zhegvd		RMG_FC_GLOBAL(zhegvd, ZHEGVD)
#define		zhegst		RMG_FC_GLOBAL(zhegst, ZHEGST)
#define		zheevd		RMG_FC_GLOBAL(zheevd, ZHEEVD)
#define		dsyev		RMG_FC_GLOBAL(dsyev, DSYEV)
#define		dsyevd		RMG_FC_GLOBAL(dsyevd, DSYEVD)
#define		dsyevx		RMG_FC_GLOBAL(dsyevx, DSYEVX)
#define		dsyevr		RMG_FC_GLOBAL(dsyevr, DSYEVR)
#define		zheev		RMG_FC_GLOBAL(zheev, ZHEEV)
#define		dtrsm		RMG_FC_GLOBAL(dtrsm, DTRSM)
#define		dsygst		RMG_FC_GLOBAL(dsygst, DSYGST)
#define		zgeev		RMG_FC_GLOBAL(zgeev, ZGEEV)
#define		zgemv		RMG_FC_GLOBAL(zgemv, ZGEMV)
#define		zdotc		RMG_FC_GLOBAL(zdotc, ZDOTC)
#define		dgels		RMG_FC_GLOBAL(dgels, DGELS)
#define		dsytrf		RMG_FC_GLOBAL(dsytrf, DSYTRF)
#define		dsytri		RMG_FC_GLOBAL(dsytri, DSYTRI)
#define		dlange		RMG_FC_GLOBAL(dlange, DLANGE)

#if __cplusplus
extern "C" {
#endif

void my_copy(double *in, double *out, int length);
void my_scal(double alpha, double *vect, int length);
void my_axpy(double alpha, double *in, double *out, int length);
void my_swap(double *vec1, double *vec2, int length);


int ilaenv (int *ispec, char *name, char *opts, int *n1, int *n2, int *n3,
            int *n4);
double dlange(const char *norm, int *m, int *n, double *A, int *lda, double *work);
void daxpy (int *n, double * alpha, double * x, int *incx, double * y, int *incy);
void dscal (int *n, double * alpha, double * x, int *incx);
void dcopy (int *n, double * x, int *incx, double * y, int *incy);
double ddot (int *n, double * x, int *incx, double * y, int *incy);
double dnrm2 (int *n, double * x, int *incx);
void dsyrk (const char *uplo, const char *trans, int *n, int *k, double * alpha, double * A,
            int *lda, double * beta, double * c, int *ldc);
void dsyr2k (const char *uplo, const char *trans, int *n, int *k, double * alpha, double * A,
            int *lda, double *b, int *ldb, double * beta, double * c, int *ldc);
void dswap (int *n, double * x, int *incx, double * y, int *incy);
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
void dpotrf (char *uplo, int *n, double * a, int *lda, int *info);
void zpotrf (char *uplo, int *n, double * a, int *lda, int *info);
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
void dpotri (char *, int *, double *, int *, int *);
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

void dsygvj(int *, char *, char *, int *, double *, int *, double *, int *, double *, double *, int*, int *, int*, int*);

void dtrsm(char *side, char *uplo, char *transa, char *diag, int *M, int *N, double *alpha, double *A, int *lda, double *B, int *ldb);
void dsygst( int *itype, char *uplo, int *N, double *A, int *LDA, double *B, int *LDB, int *INFO );
double dzasum(int *, double *A, int *);
void dsygvd(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *eigs, double *work, int *lwork, int *iwork, int *liwork, int *info);
void zhegvd(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, double *eigs, double *work, int *lwork, double
*rwork, int *lrwork, int *iwork, int *liwork, int *info);
void zcopy(int *, DoubleC *, int *, DoubleC *, int *);
void zaxpy(int *, DoubleC *, DoubleC *, int *, DoubleC *, int*);
void zgeev(char *, char *, int *, DoubleC*, int *, DoubleC*, DoubleC*, int *, DoubleC*, int *, DoubleC*, int *, double *, int * );


void zgemv (char *, int *, int *, DoubleC *, DoubleC*, int *, DoubleC*, int *, DoubleC *,DoubleC *,  int *);
void dgels(char *trans, int *M, int *N, int *nrhs, double *A, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);
void dsytrf(char *, int *, double *, int *, int *, double *, int *, int *);
void dsytri(char *, int *, double *, int *, int *, double *, int *);
void dger(int *, int *, double *, double *, int *, double *, int *, double *, int *);


DoubleC zdotc(int*, DoubleC *, int*, DoubleC *, int*);

void openblas_set_num_threads(int nthreads);

#if __cplusplus
}
#endif
/******/
