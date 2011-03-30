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

#  define		global_sums     global_sums_
#  define		saxpy		daxpy_
#  define		scopy   	dcopy_
#  define		dcopy   	dcopy_
#  define		sdot    	ddot_
#  define		sscal   	dscal_
#  define		ssyrk   	dsyrk_
#  define		spotrf  	dpotrf_
#  define		spotri  	dpotri_
#  define		cpotrf  	zpotrf_
#  define		cpotri  	zpotri_
#  define		snrm2   	dnrm2_
#  define		sswap   	dswap_
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
#  define		latgen 		latgen_
#  define		fgram   	dfgram_
#  define		dlamch  	dlamch_
#  define		dgemm  	        dgemm_
#  define		zgemm  	        zgemm_
#  define		dgetrf 	        dgetrf_
#  define		dgetri 	        dgetri_
#  define		dgesv 	        dgesv_
#  define               sgesv           dgesv_
#  define		zgesv 	        zgesv_

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

/* Macro to simplify calling dcopy*/
#define my_copy(in, out, length)  do{ int my_one = 1; int my_length = length;  dcopy(& (my_length),  (in), & (my_one),  (out), & (my_one)); }while(0)

void fsymforces (REAL * force, int *s, int *irg, int *irt,
                 int *nat, int *ibrav, int *nsym,
                 REAL * celldm, int *nr1, int *nr2, int *nr3);
int ilaenv (int *ispec, char *name, char *opts, int *n1, int *n2, int *n3,
            int *n4);
void saxpy (int *n, REAL * alpha, REAL * x, int *incx, REAL * y, int *incy);
void sscal (int *n, REAL * alpha, REAL * x, int *incx);
void scopy (int *n, REAL * x, int *incx, REAL * y, int *incy);
REAL sdot (int *n, REAL * x, int *incx, REAL * y, int *incy);
REAL snrm2 (int *n, REAL * x, int *incx);
void ssyrk (char *uplo, char *trans, int *n, int *k, REAL * alpha, REAL * A,
            int *lda, REAL * beta, REAL * c, int *ldc);
void sswap (int *n, REAL * x, int *incx, REAL * y, int *incy);
void xbecke (REAL * d, REAL * s, REAL * u, REAL * v, REAL * ex, REAL * vx);
void corlyp_f90 (REAL * dp, REAL * dm, REAL * dp1, REAL * dm1, REAL * dp2,
                 REAL * dm2, REAL * ec, REAL * vcp0, REAL * vcm0, int *ndm);

void corlsd (REAL * rs, REAL * zet, REAL * ec, REAL * vcup, REAL * vcdn,
             REAL * ecrs, REAL * eczet, REAL * alfc);
void corgga (REAL * rs, REAL * zet, REAL * t, REAL * uu, REAL * vv,
             REAL * ww, REAL * h, REAL * dvcup, REAL * dvcdn,
             REAL * fk, REAL * sk, REAL * g, REAL * ec, REAL * ecrs,
             REAL * eczet);
void corpbe (REAL * rs, REAL * zet, REAL * t, REAL * uu, REAL * vv, REAL * ww,
             int *lgga, int *lpot, REAL * ec, REAL * vcup, REAL * vcdn,
             REAL * h, REAL * dvcup, REAL * DVCDN);
void cpotrf (char *uplo, int *n, REAL * a, int *lda, int *info);
void exch (REAL * d, REAL * s, REAL * u, REAL * v, REAL * ex, REAL * vx);
void exchpbe (REAL * d, REAL * s, REAL * u, REAL * v, int *lgga, int *lpot,
              REAL * ex, REAL * vx);
void symmetry (int *ibrav, int *s, int *nsym, int *irg, int *irt,
               int *ftau, int *nat, REAL * tau, int *ityp, int *nks,
               REAL * xk, REAL * wk, REAL * celldm, int *nr1, int *nr2,
               int *nr3, int *wflag);
void symrho (REAL * rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
             int *irg, int *ftau);
void ssyev (char *jobz, char *uplo, int *numm, REAL * ss, int *numn,
            REAL * work1, REAL * work2, int *lwork, int *info);
void cheev (char *jobz, char *uplo, int *numm, REAL * ss, int *numn,
            REAL * work1, REAL * work2, int *lwork, REAL *, int *info);
void spotri (char *, int *, REAL *, int *, int *);
void dgemm(char *, char *, int *, int *, int *, REAL *, REAL *, int *, REAL *, int *, REAL *, REAL *, int *);
void zgemm(char *, char *, int *, int *, int *, REAL *, REAL *, int *, REAL *, int *, REAL *, REAL *, int *);
void dgetrf( int *, int *, REAL *, int *, int *, int *);
void dgetri(int *, REAL *, int *, int *, REAL *, int *, int *);
void dgesv (int *, int*, REAL *, int *, int *, REAL *, int *, int *);
void zgesv (int *, int*, REAL *, int *, int *, REAL *, int *, int *);

/******/
