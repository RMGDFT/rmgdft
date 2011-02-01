/************************** SVN Revision Information **************************
 **    $Id: blas3_def.h 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 

#if CRAY_T3E

#include <fortran.h>
#define CCHAR    _fcd

#else

#define CCHAR    char*

#endif




#if (AIXFLAG || AIX_MPI)

#define               ssyrk   dsyrk
#define               sgemv   dgemv
#define               ssymm   dsymm
#define               sgemm   dgemm
#define               ssymv   dsymv
#define               dsymv   dsymv
#define               dsymm   dsymm
#define               dgemm   dgemm
#define               dgemv   dgemv
#define               dsyrk   dsyrk
#define               dtrmm   dtrmm

#endif

#if ( SGIFLAG || LINUX || LINUX_MPI )

#define               sgemm   dgemm_
#define               ssymm   dsymm_
#define               ssyrk   dsyrk_
#define               ssymv   dsymv_
#define               sgemv   dgemv_
#define	        dsymv	dsymv_
#define	        dsymm   dsymm_
#define	        dgemm   dgemm_
#define	        dgemv   dgemv_
#define	        dsyrk   dsyrk_
#define	        dtrmm   dtrmm_

#endif

#if CRAY_T3E

#define               sgemm   SGEMM
#define               sgemvx  SGEMVX
#define               ssymv   SSYMV
#define               sgemv   SGEMV
#define               sgesv   SGESV

#define               dgemm   SGEMM
#define               dgemvx  SGEMVX
#define               dsymv   SSYMV
#define               dgemv   SGEMV
#define               dgesv   SGESV

#endif


void ssymv (CCHAR, int *, double *, double *, int *lda,
            double *x, int *incx, double *beta, double *y, int *incy);
void sgemv (CCHAR, int *, int *, double *, double *, int *lda,
            double *x, int *incx, double *beta, double *y, int *incy);
void ssyrk (CCHAR, CCHAR, int *, int *,
            double *, double *, int *, double *, double *, int *);
void sgemm (CCHAR, CCHAR, int *, int *, int *,
            double *, double *, int *, double *, int *,
            double *, double *, int *);
void ssymm (CCHAR, CCHAR, int *, int *,
            double *, double *, int *, double *, int *,
            double *, double *, int *);
void dgemv (CCHAR, int *, int *, double *, double *, int *lda,
            double *x, int *incx, double *beta, double *y, int *incy);
void dsyrk (CCHAR, CCHAR, int *, int *,
            double *, double *, int *, double *, double *, int *);
void dgemm (CCHAR, CCHAR, int *, int *, int *,
            double *, double *, int *, double *, int *,
            double *, double *, int *);
void dsymm (CCHAR, CCHAR, int *, int *,
            double *, double *, int *, double *, int *,
            double *, double *, int *);
void dtrmm (char *, char *, char *, char *, int *, int *, double *, double *,
            int *, double *, int *);


void dsymm_c (char *, char *, int *, int *, double *, double *, int *,
              double *, int *, double *, double *, int *);
void dgemm_c (char *, char *, int *, int *, int *,
              double *, double *, int *, double *,
              int *, double *, double *, int *);
void dgemv_c (char *, int *, int *,
              double *, double *, int *, double *,
              int *, double *, double *, int *);
void ssyrk_c (char *, char *, int *, int *, double *, double *,
              int *, double *, double *, int *);
void dtrmm_c (char *, char *, char *, char *, int *, int *, double *,
              double *, int *, double *, int *);
