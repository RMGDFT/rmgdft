/************************** SVN Revision Information **************************
 **    $Id: my_scalapack.h 2073 2013-11-02 18:19:59Z ebriggs $    **
******************************************************************************/
 
#include "blacs.h"



#define _fcd char*


#ifdef AIX
#define  	NUMROC  	numroc
#define		INDXG2P  	indxg2p
#define 	DESCINIT 	descinit
#define 	PDSYEV 		pdsyev
#define 	PCHEEV 		pzheev
#define 	PSPOCON		pdpocon
#define 	PSPOTRF		pdpotrf
#define 	PSPOTRI		pdpotri
#define 	PSSYGST		pdsygst
#define 	PSSYEV		pdsyev
#define 	PSTRTRS		pdtrtrs
#define 	PDGEMM 		pdgemm
#define 	PSGEMM 		pdgemm
#define 	PZGEMM 		pzgemm
#define 	PSSYMM 		pdsymm
#define 	PDGESV 		pdgesv
#define 	PZGESV 		pzgesv
#define 	PSUBDIAG 	psubdiag
#define 	PDSYGVX 	pdsygvx
#define         PDSYEVX         pdsyevx
#define         PDPOTRF         pdpotrf
#define         PDSYNGST        pdsyngst
#define         PDTRSM          pdtrsm
#define 	PDTRAN  	pdtran




#else

#ifdef  LINUX

#define  	NUMROC  	numroc_
#define		INDXG2P  	indxg2p_
#define 	DESCINIT 	descinit_
#define 	PDSYEV 		pdsyev_
#define 	PCHEEV 		pzheev_
#define 	PSPOCON		pdpocon_
#define 	PSPOTRF		pdpotrf_
#define 	PSPOTRI		pdpotri_
#define 	PSSYGST		pdsygst_
#define 	PSTRTRS		pdtrtrs_
#define 	PDGEMM 		pdgemm_
#define 	PSGEMM 		pdgemm_
#define 	PZGEMM 		pzgemm_
#define 	PSSYMM 		pdsymm_
#define 	PSSYEV		pdsyev_
#define 	PDGESV 		pdgesv_
#define 	PZGESV 		pzgesv_
#define 	PSUBDIAG 	psubdiag_
#define 	PDSYGVX 	pdsygvx_
#define 	PDSYEVX 	pdsyevx_
#define 	PDPOTRF 	pdpotrf_
#define 	PDSYNGST 	pdsyngst_
#define 	PDTRSM  	pdtrsm_
#define 	PDTRAN  	pdtran_

#endif
#endif


int NUMROC (int *, int *, int *, int *, int *);
int INDXG2P (int *, int *, int *, int *, int *);
void DESCINIT (int[], int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void PDGESV(int *, int *, rmg_double_t *, int * , int *, int *, int *, rmg_double_t *,
	int *, int *, int *, int *);
void PZGESV(int *, int *, rmg_double_t *, int * , int *, int *, int *, rmg_double_t *,
	int *, int *, int *, int *);
void PDGEMM (_fcd, _fcd, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void PZGEMM (_fcd, _fcd, int *, int *, int *, double *, double *, int *,
             int *, int *, double *, int *, int *, int *, double *, double *,
             int *, int *, int *);
void PDSYEV (_fcd, _fcd, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, int *);
void PCHEEV (_fcd, _fcd, int *, double *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, int *, rmg_double_t *, int *, int *);
void PSPOCON (_fcd, int *, double *, int *, int *, int *, double *, double *,
              double *, int *, int *, int *, int *);
void PSPOTRF (_fcd, int *, double *, int *, int *, int *, int *);
void PSPOTRI (_fcd, int *, double *, int *, int *, int *, int *);
void PSSYGST (int *, _fcd, int *, double *, int *, int *, int *, double *,
              int *, int *, int *, double *, int *);
void PSTRTRS (_fcd, _fcd, _fcd, int *, int *, double *, int *, int *, int *,
              double *, int *, int *, int *, int *);
void PSSYMM (_fcd, _fcd, int *, int *, double *, double *, int *, int *,
             int *, double *, int *, int *, int *, double *, double *, int *,
             int *, int *);

void PSUBDIAG (char *, char *, int, double *, int, double *, int *);
void PDSYGVX(int *, char*, char*, char*, int*, rmg_double_t *, int*, int*, int*, rmg_double_t*, int*, int*, 
       int*, rmg_double_t*, rmg_double_t *, int*, int*, rmg_double_t*, int*, int*, rmg_double_t*, rmg_double_t*, rmg_double_t*, int*, 
       int*, int*, rmg_double_t*, int*, int*, int*, int*, int*, rmg_double_t*, int*);       
void PDSYEVX(char*, char*, char*, int*, rmg_double_t *, int*, int*, int*, rmg_double_t*, rmg_double_t*, int*,
       int*, rmg_double_t*, int*, int*, rmg_double_t*, rmg_double_t*, rmg_double_t*, int*,                                           
       int*, int*, rmg_double_t*, int*, int*, int*, int*, int*, rmg_double_t*, int*);



void matgather (double *, int *, double *, int);
void dsymm_dis (char *, char *, int *, double *, double *, double *);
void proc_gridsetup(int nproc, int *nprow, int *npcol);


/* Blacs dimension */
#define DLEN    9
