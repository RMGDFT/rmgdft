/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#ifndef INCLUDE_METHOD
#define INCLUDE_METHOD

#define		DEBUG		0


/* Use nonorthogonal functions (1) or Ritz functions (0) */
#define		NONORTHO	1

/* Mehrstellenverfahren FD with positive definite r.h.s. operator */
#define		MPD		0

/* Use regular 4th order FD scheme instead of MFD */
#define         FD4             0

/* Use distributed matrices (PBLAS and ScaLapack) */
#define		USE_DIS_MAT	( 1 && NONORTHO )


/* line minimization for wavefunctions corrections */
#define         LM              ( 1 && NONORTHO )

/* conjugate gradient directions for wavefunctions corrections */
#define         CG              ( 1 && LM )



/* Update Hij after correction of the eigenvectors (1) or not (0)*/
/* Note: should be at 1 to reproduce NONORTHO=0 when NONORTHO=1 */
#define	UPDATE_H	(1 && NONORTHO)


#endif
