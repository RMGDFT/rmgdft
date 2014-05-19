/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/* 
 *     SD (steepest decent) method subroutine 
 *
 *
 * 	Input: 	step-- iteration step for KAIN,  MUST START from ZERO 
 * 		N  ---	Dinmention of array  xm and fm 	
 * 		xm ---	current sollution 
 * 		fm ---	current residual 
 *
 * 	Output: xm ----- updated 
 * 
 *
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "prototypes_on.h"



void precond(double *x);


void sd(int step, int N, double *xm, double *fm)
{
    int ione = 1;
    double gamma = -0.5;

    precond(fm);
    saxpy(&N, &gamma, fm, &ione, xm, &ione);
}
