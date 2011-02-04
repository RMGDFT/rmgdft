/************************** SVN Revision Information **************************
 **    $Id: lapack_def.h 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 

#if CRAY_T3E

#include <fortran.h>
#define CCHAR    _fcd
#define exit	globalexit

  #define               dposv   SPOSV
  #define		dsyev	SSYEV
  #define               dgeev   SGEEV
  #define		dpotri	SPOTRI
  #define		dpotrf	SPOTRF
  #define		dsygst  SSYGST
  #define		dtrtrs	STRTRS
  #define		dpocon	SPOCON
  #define		dsygv	SSYGV



#else

#define CCHAR    char*

#endif


/*  #if LINUX */
#if  (LINUX || LINUX_MPI ||M_SGI_ORIGIN_MPI) 


  #define               dposv   dposv_
  #define		dsyev	dsyev_
  #define               dgeev   dgeev_
  #define		dpotri	dpotri_
  #define		dpotrf	dpotrf_
  #define               dpocon  dpocon_
  #define		dsygst  dsygst_
  #define		dtrtrs	dtrtrs_
  #define		dpocon	dpocon_
  #define		dsygv	dsygv_

#endif




/*  #if AIXFLAG  */
#if (AIX || AIX_MPI || AIX_SMP ) 
#define dpotrf_ dpotrf
#define dpotri_ dpotri
#define dgelss_ dgelss
#define dgelsx_ dgelsx
#define dgels_  dgels
#define dposv_  dposv
#define dsyev_  dsyev
#define dgeev_  dgeev
#define dsygv_  dsygv
#define dtrsv_  dtrsv
#define dgesv_  dgesv

#endif





void sgelss(int *,int *,int *,double *,int *,
             double *,int *,
             double*,double*,int *,double*,int*,int*);
void sgelsx(int *,int *,int *,double *,int *,double *,int *,int *,
             double *,int *,double *,int *);


void dgsls(CCHAR, int*, int*, int*, double *, int *,
            double *, int*, double*, int*, int *);
void dposv(CCHAR,int *,int *,double *,int *,
                                        double *,int *,int *);
void dsyev(CCHAR, CCHAR, int *, double *, int *, 
              double *, double *, int *, int *);
void dgeev(CCHAR, CCHAR, int*, double*, int*, double*,
    double*, double*, int*, double*, int*, double*, int*, int*);
void dsygv(int*,CCHAR,CCHAR,int*,double*,int*,double*,
    int*,double*,double*,int*,int*);
void dpotri(CCHAR, int*, double*, int*, int*);
void dpotrf(CCHAR, int*, double*, int*, int*);
void dpocon(CCHAR, int *, double *, int *, double *,
              double *, double *, int *,
             int *);
void dtrtrs(CCHAR, CCHAR, CCHAR, int*, int*, double*, int*, double*, int*, int*);
void dsygst(int*, CCHAR, int*, double*, int*, double*, int*, int*);






void dposv_c(char*, int*, int*, double*, int*, 
             double*, int*, int*);
void dsyev_c(char*,char*,int*,double*,int*,double*,double*,int*,int*);
void dgees_c_(char*, char*, char*, int*, double*, int*, int*,
              double*, double*, double*, int*, double*, double*,
              double*, int*, int*);
void dsygst_c(int *, char *, int *,
             double *, int *, double *,
             int *l, int *);
void dtrtrs_c(char *, char *, char *,
             int *, int *, double *, int *,
             double *, int *, int *);
void dpotrf_c(char *, int *, double *, int *,
             int *);
void dpotri_c(char *, int *, double *, int *,
             int *);
void dpocon_c(char *, int *, double *, int *, double *,
              double *, double *, int *,
             int *);
void dsygv_c(int *itype, char *jobz, char *uplo, int *n, double *a, 
             int *lda, double *b, int *ldb, double *work1, double*work2, 
             int *lwork, int *info);
