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

#if  LINUX || __CYGWIN__

#define  	NUMROC  	numroc_
#define		INDXG2P  	indxg2p_
#define 	DESCINIT 	descinit_
#define 	PDSYEV 		pdsyev_
#define 	PCHEEV 		pzheev_
#define 	PSPOCON		pdpocon_
#define 	PSPOTRF		pdpotrf_
#define 	PSPOTRI		pdpotri_
#define 	PSTRTRS		pdtrtrs_
#define 	PDGEMM 		pdgemm_
#define 	PSGEMM 		pdgemm_
#define 	PZGEMM 		pzgemm_
#define 	pzgemm 		pzgemm_
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
#define 	PZTRANC  	pztranc_
#define 	PZHEGVX         pzhegvx_
#define 	pdgetrf     pdgetrf_
#define 	pdgetrs     pdgetrs_
#define 	pdsygvx     pdsygvx_
#define     pdlmach     pdlmach_
#define     dlmach      dlmach_
#define     idamax      idamax_

#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif


int idamax_(int *, double *, int*);
double pdlamch(int *, char *);
int NUMROC (int *, int *, int *, int *, int *);
int INDXG2P (int *, int *, int *, int *, int *);
void DESCINIT (int[], int *, int *, int *, int *, int *, int *, int *, int *,
               int *);
void PDGESV(int *, int *, double *, int * , int *, int *, int *, double *,
	int *, int *, int *, int *);
void PZGESV(int *, int *, double *, int * , int *, int *, int *, double *,
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
             double *, int *, int *, int *, double *, int *, double *, int *, int *);
void PSPOCON (_fcd, int *, double *, int *, int *, int *, double *, double *,
              double *, int *, int *, int *, int *);
void PSPOTRF (_fcd, int *, double *, int *, int *, int *, int *);
void PSPOTRI (_fcd, int *, double *, int *, int *, int *, int *);
void PSTRTRS (_fcd, _fcd, _fcd, int *, int *, double *, int *, int *, int *,
              double *, int *, int *, int *, int *);
void PSSYMM (_fcd, _fcd, int *, int *, double *, double *, int *, int *,
             int *, double *, int *, int *, int *, double *, double *, int *,
             int *, int *);

void PSUBDIAG (char *, char *, int, double *, int, double *, int *);
void PDSYGVX(int *, char*, char*, char*, int*, double *, int*, int*, int*, double*, int*, int*, 
       int*, double*, double *, int*, int*, double*, int*, int*, double*, double*, double*, int*, 
       int*, int*, double*, int*, int*, int*, int*, int*, double*, int*);       
void PZHEGVX(int *, char*, char*, char*, int*, double *, int*, int*, int*, double*, int*, int*, 
       int*, double*, double *, int*, int*, double*, int*, int*, double*, double*, double*, int*, 
       int*, int*, double*, int *, double *, int*, int*, int*, int*, int*, double*, int*);       
void PDSYEVX(char*, char*, char*, int*, double *, int*, int*, int*, double*, double*, int*,
       int*, double*, int*, int*, double*, double*, double*, int*,                                           
       int*, int*, double*, int*, int*, int*, int*, int*, double*, int*);
void pzgetri(int *, double *, int *, int*, int*, int*, double*, int*, int*, int*,int*);
void pzgetrf(int *, int*,  double *, int *, int*, int*, int*,int*);
void pdgetrf(int *, int*,  double *, int *, int*, int*, int*,int*);
void pdgetrs(char *, int *, int*,  double *, int *, int*, int*, int*, double *, int*, int*, int*, int*);
void pztranu_(int *, int*,  double *, double *, int *, int*, int*, double *, double *, int *, int*,int*);
void pdsygvx( int *, char *, char *, char *, int *, double *, int *, int *,
                 int *, double *, int *, int *, int *, double *, double *, int *, int *,
                 double *, int *, int *, double *, double *, double *, int *, int *, int *,
                 double *, int *, int *, int *, int *, int *,
                 double *, int * );

void pdtran_(int *, int*,  double *, double *, int *, int*, int*, double *, double *, int *, int*,int*);

void matgather (double *, int *, double *, int);
void dsymm_dis (char *, char *, int *, double *, double *, double *);
void proc_gridsetup(int nproc, int *nprow, int *npcol); 

void Cpdgemr2d(int n, int m, double *A, int ia, int ja, int *desca, 
        double *B, int ib, int jb, int *descb, int ictxt);
int indxl2g(int *, int *, int*, int *, int*);

#ifdef __cplusplus
}
#endif

/* Blacs dimension */
#define DLEN    9
