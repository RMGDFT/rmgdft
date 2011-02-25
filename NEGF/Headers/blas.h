/************************** SVN Revision Information **************************
 **    $Id: blas.h 1242 2011-02-02 18:55:23Z luw $    **
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

#if (LINUX || LINUX_SMP || LINUX_MPI || SGI_ORIGIN || M_SGI_ORIGIN_MPI || SOLARIS)

  #define               global_sums     global_sums_
  #define		app6_del2	app6_del2_
  #define		app_del2f	app_del2f_
  #define		sgesv	dgesv_
  #define		ssum	dsum_
  #define		saxpy	daxpy_
  #define		scopy   dcopy_
  #define		zcopy   zcopy_
  #define		dcopy   dcopy_
  #define		DCOPY   dcopy_
  #define		sdot    ddot_
  #define		ddot    ddot_
  #define		sscal   dscal_
  #define		dscal   dscal_
  #define               sgemvx  dgemvx_
  #define               sgemv   dgemv_
  #define               sgemm   dgemm_
  #define               cgemm   zgemm_
  #define               zgemm   zgemm_
  #define               zgesv   zgesv_
  #define               sger    dger_
  #define               ssyrk   dsyrk_
  #define               spotrf  dpotrf_
  #define               spotri  dpotri_
  #define               cpotrf  zpotrf_
  #define               cpotri  zpotri_
  #define               sgerx   dgerx_
  #define               snrm2   dnrm2_
  #define               dznrm2  dznrm2_
  #define               sswap   dswap_
  #define		ssyev	dsyev_
  #define		cheev	zheev_
  #define		xbecke	xbecke_
  #define               corlyp_f90      corlyp_f90_
  #define		exch	exch_
  #define		corlsd	corlsd_
  #define		corgga	corgga_
  #define		exchpbe exchpbe_
  #define		corpbe  corpbe_
  #define               symrho  symrho_
  #define               fsymforces  fsymforces_
  #define               symmetry symmetry_
  #define               latgen latgen_
  #define               fgram   dfgram_
  #define 		ZGESV zgesv_
  #define 		ZCOPY zcopy_
  #define 		ZAXPY zaxpy_
  #define 		ZGEMM zgemm_
  
#endif

#if (AIX || AIX_MPI )

  #define		sgesv	dgesv
  #define		ssum	dsum
  #define		saxpy	daxpy
  #define		scopy   dcopy
  #define		SCOPY   dcopy
  #define		DCOPY   dcopy
  #define		sdot    ddot
  #define		sscal   dscal
  #define               sgemvx  dgemvx
  #define               sgemv   dgemv
  #define               sgemm   dgemm
  #define               cgemm   zgemm
  #define               sger    dger
  #define               ssyrk   dsyrk
  #define               spotrf  dpotrf
  #define               spotri  dpotri
  #define               sgerx   dgerx
  #define               snrm2   dnrm2
  #define               sswap   dswap
  #define		ssyev	dsyev
  #define		xbecke	xbecke
  #define               corlyp_f90      corlyp_f90
  #define		exch	exch
  #define		corlsd	corlsd
  #define		corgga	corgga
  #define               cpotrf  zpotrf
  #define               cpotri  zpotri
  #define               cheev   zheev
  #define               fgram   dfgram
  #define 		ZGESV zgesv
  #define 		ZCOPY zcopy
  #define 		ZAXPY zaxpy
  #define 		ZGEMM zgemm
  #define 		dGEMM dgemm

#endif

#if (CRAY_YMP || CRAY_C90 || CRAY_T3E)

  #define		saxpy	SAXPY
  #define		scopy   SCOPY
  #define		sdot    SDOT
  #define		sscal   SSCAL
  #define               ssum    SSUM
  #define               sgemm   SGEMM
  #define               cgemm   CGEMM
  #define               sgemvx  SGEMVX
  #define               sgemv   SGEMV
  #define               sgesv   SGESV
  #define               sgerx   SGERX
  #define               sger    SGER
  #define               snrm2   SNRM2
  #define               sswap   SSWAP
  #define               spotrf  SPOTRF
  #define               spotri  SPOTRI
  #define               ssyev   SSYEV
  #define               global_sums GLOBAL_SUMS
  #define               ilaenv  ILAENV
  #define               app_del2f  APP_DEL2F
  #define               xbecke    XBECKE
  #define               corlyp_f90      corlyp_f90
  #define               exch      EXCH
  #define               corlsd    CORLSD
  #define               corgga    CORGGA
  #define               corpbe    CORPBE
  #define               exchpbe   EXCHPBE
  #define               symrho    SYMRHO  
  #define               symmetry  SYMMETRY 
  #define               fsymforces FSYMFORCES
  #define 		cpotrf     CPOTRF
  #define 		cpotri     CPOTRI
  #define 		cheev      CHEEV
  #define 		fgram      FGRAM
  #define               latgen 	   LATGEN

#endif

  void fsymforces(REAL *force, int *s, int *irg, int *irt,
                         int *nat, int *ibrav, int *nsym, 
                         REAL *celldm, int *nr1, int *nr2, int *nr3);
  int ilaenv(int *ispec, char *name, char *opts, int *n1, int *n2, int *n3,
             int *n4);
  void saxpy(int *n, REAL *alpha, REAL *x, int *incx, REAL *y, int *incy);
  void sscal(int *n, REAL *alpha, REAL *x, int *incx);
  void scopy(int *n, REAL *x, int *incx, REAL *y, int *incy);
  REAL sdot(int *n, REAL *x, int *incx, REAL *y, int *incy);
  REAL snrm2(int *n, REAL *x, int *incx);
  REAL dznrm2(int *n, REAL *x, int *incx);
  void sgemvx(int *m, int *n, REAL *alpha, REAL *a, int *ax, int *ay,
              REAL *b, int *bx, REAL *beta, REAL *c, int *cx);
  void sger(int *m, int *n, REAL *alpha, REAL *x, int *incx, REAL *y, int *incy, REAL *a, int *lda);       
  REAL ssum(int *n, REAL *x, int *incx);
  void ssyrk(char *uplo, char *trans, int *n, int *k, REAL *alpha, REAL *A,
             int *lda, REAL *beta, REAL *c, int *ldc);
  void sswap(int *n, REAL *x, int *incx, REAL *y, int *incy);
  void app_del2f(REAL *a, REAL *b, int *dimx, int *dimy, int *dimz,
                REAL *gridhx, REAL *gridhy, REAL *gridhz);
  void xbecke(REAL *d, REAL *s, REAL *u, REAL *v, REAL *ex, REAL *vx);
  void corlyp_f90(REAL *dp,REAL *dm, REAL *dp1, REAL *dm1, REAL *dp2,REAL *dm2,REAL *ec,REAL *vcp0,REAL *vcm0,int *ndm);

  void corlsd(REAL *rs, REAL *zet, REAL *ec, REAL *vcup, REAL *vcdn, 
              REAL *ecrs, REAL *eczet, REAL *alfc);
  void corgga(REAL *rs, REAL *zet, REAL *t, REAL *uu, REAL *vv,
              REAL *ww, REAL *h, REAL *dvcup, REAL *dvcdn,
              REAL *fk, REAL *sk, REAL *g, REAL *ec, REAL *ecrs, REAL *eczet);
  void corpbe(REAL *rs, REAL *zet, REAL *t, REAL *uu, REAL *vv,
              REAL *ww, int *lgga, int *lpot, REAL *ec, REAL *vcup,
              REAL *vcdn, REAL *h, REAL *dvcup, REAL *DVCDN);
  void exch(REAL *d, REAL *s, REAL *u, REAL *v, REAL *ex, REAL *vx);
  void exchpbe(REAL *d, REAL *s, REAL *u, REAL *v, int *lgga, int *lpot, REAL *ex, REAL *vx);
  void sgesv(int *N, int *NRHS, REAL *A, int *LDA, int *IPIV, REAL *b, 
             int *LDB, int *INFO);
  void symmetry(int *ibrav, int *s, int *nsym, int *irg, int *irt,
	        int *ftau, int *nat, REAL *tau, int *ityp, int *nks,
	        REAL *xk, REAL *wk, REAL *celldm, int *nr1, int *nr2, int *nr3, 
                int *wflag);
  void symrho(REAL *rho, int *nr1, int *nr2, int *nr3, int *nsym, int *s,
              int *irg, int *ftau);
  void ssyev(char *jobz, char *uplo, int *numm, REAL *ss, int *numn, REAL *work1,
		  REAL *work2, int *lwork, int *info);
  void spotri(char*,int*,REAL*,int*,int*);
  void zgesv(int*, int*, doublecomplex *, int*, int*, doublecomplex*, int*, int*);
#if !(CRAY_YMP || CRAY_C90 || CRAY_T3E)   
  void cpotrf(char *uplo, int *n, REAL *a, int *lda, int *info); 
#endif
/******/
