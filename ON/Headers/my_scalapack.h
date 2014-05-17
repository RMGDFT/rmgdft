/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include "params.h"
#include "blacs.h"

#if CRAY_T3E

#include <fortran.h>

#else

#define _fcd char*
/*
char *_cptofcd(char*, int);
*/
#endif


#if (AIX || IBMSP || AIX_SMP || AIX_MPI )

#define NUMROC numroc
#define INDXG2P indxg2p
#define DESCINIT descinit
#define PSSYEV pdsyev
#define PSPOCON pdpocon
#define PSPOTRF pdpotrf
#define PSPOTRI pdpotri
#define PSSYGST pdsygst
#define PSTRTRS pdtrtrs
#define PSGEMM  pdgemm
#define PSSYMM pdsymm
#define PDTRAN pdtran

void Cpdgemr2d(int m, int n,
                double* a, int ia, int ja, int* desca,
                double* b, int ib, int jb, int* descb,
                int gcontext);

void pdsyev (_fcd, _fcd, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, int *);
void pdpocon (_fcd, int *, double *, int *, int *, int *, double *, double *,
              double *, int *, int *, int *, int *);
void pdpotrf (_fcd, int *, double *, int *, int *, int *, int *);
void pdpotri (_fcd, int *, double *, int *, int *, int *, int *);
void pdsygst (int *, _fcd, int *, double *, int *, int *, int *, double *,
              int *, int *, int *, double *, int *);
void pdtrtrs (_fcd, _fcd, _fcd, int *, int *, double *, int *, int *, int *,
              double *, int *, int *, int *, int *);
void pdgemm (_fcd, _fcd, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void pdsymm (_fcd, _fcd, int *, int *, double *, double *, int *, int *,
             int *, double *, int *, int *, int *, double *, double *, int *,
             int *, int *);


/* Function prototypes */
void init_scalapack (void);
void psubdiag (char *, char *, int, double *, int, double *, int *);
void matinit (double *, int *, double *, int);
void matgather (double *, int *, double *, int);
void DESCINIT (int*, int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void sl_init (int *, int, int);
void distribute_mat (double *, double *);
void dsymm_dis (char *, char *, int *, double *, double *, double *);

int numroc (int *, int *, int *, int *, int *);
int indxg2p (int *, int *, int *, int *, int *);



#endif

#if  ( LINUX)

#define  	NUMROC  	numroc_
#define		INDXG2P  	indxg2p_
#define 	SL_INIT 	sl_init_
#define 	DESCINIT 	descinit_
#define 	PSSYEV 		pdsyev_
#define 	PSPOCON		pdpocon_
#define 	PSPOTRF		pdpotrf_
#define 	PSPOTRI		pdpotri_
#define 	PSSYGST		pdsygst_
#define 	PSTRTRS		pdtrtrs_
#define 	PSGEMM 		pdgemm_
#define 	PSSYMM 		pdsymm_
#define 	PDTRAN 		pdtran_

#define 	PSUBDIAG 	psubdiag_

void Cpdgemr2d(int m, int n,
                double* a, int ia, int ja, int* desca,
                double* b, int ib, int jb, int* descb,
                int gcontext);


int NUMROC (int *, int *, int *, int *, int *);
int INDXG2P (int *, int *, int *, int *, int *);
void SL_INIT (int *, int, int);
void DESCINIT (int[], int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void PSSYEV (_fcd, _fcd, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, int *);
void PSPOCON (_fcd, int *, double *, int *, int *, int *, double *, double *,
              double *, int *, int *, int *, int *);
void PSPOTRF (_fcd, int *, double *, int *, int *, int *, int *);
void PSPOTRI (_fcd, int *, double *, int *, int *, int *, int *);
void PSSYGST (int *, _fcd, int *, double *, int *, int *, int *, double *,
              int *, int *, int *, double *, int *);
void PSTRTRS (_fcd, _fcd, _fcd, int *, int *, double *, int *, int *, int *,
              double *, int *, int *, int *, int *);
void PSGEMM (_fcd, _fcd, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void PSSYMM (_fcd, _fcd, int *, int *, double *, double *, int *, int *,
             int *, double *, int *, int *, int *, double *, double *, int *,
             int *, int *);

void PSUBDIAG (char *, char *, int, double *, int, double *, int *);


#endif

void matinit (double *, int *, double *, int);
void matgather (double *, int *, double *, int);
void distribute_mat (double *, double *);
void dsymm_dis (char *, char *, int *, double *, double *, double *);


/* Blacs dimension */
#define DLEN    9

/* bolck size, it is defined in init_dimensions.c */
int NB;
