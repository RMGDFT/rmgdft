/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/rmg_fastrelax.c *****
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
 *   Calculates constants alpha and beta for wave function and charge extrapolation
 *   between MD steps. See PRB 45, 1538 (1992) and Comp. Phys. Comm. 118. 31 (1999).
 * INPUTS
 *   no explict input
 * OUTPUT
 *   atomic coordinates are updated
 * PARENTS
 *   
 * CHILDREN
 *   
 * SOURCE
 */

#define VERBOSE 0

#include "main.h"
void get_extrapolation_constants (rmg_double_t *alpha, rmg_double_t *beta)
{
    int ion;
    ION *iptr;
    rmg_double_t d01x, d01y, d01z, d12x, d12y, d12z, d23x, d23y, d23z, detA;
    rmg_double_t a11, a12, a22, b1, b2, a21;


    a11 = 0.0;
    a12 = 0.0;
    a21 = 0.0;
    a22 = 0.0;
    b1 = 0.0;
    b2 = 0.0;
    
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        
	/* Get ion pointer */
        iptr = &ct.ions[ion];
	
	d01x = iptr->crds[0] - iptr->ocrds1[0];
	d01y = iptr->crds[1] - iptr->ocrds1[1];
	d01z = iptr->crds[2] - iptr->ocrds1[2];

	d12x = iptr->ocrds1[0] - iptr->ocrds2[0];
	d12y = iptr->ocrds1[1] - iptr->ocrds2[1];
	d12z = iptr->ocrds1[2] - iptr->ocrds2[2];
	
	d23x = iptr->ocrds2[0] - iptr->ocrds3[0];
	d23y = iptr->ocrds2[1] - iptr->ocrds3[1];
	d23z = iptr->ocrds2[2] - iptr->ocrds3[2];


	a11 += d12x*d12x + d12y*d12y + d12z*d12z;

	a12 += d12x*d23x + d12y*d23y + d12z*d23z;

	a22 += d23x*d23x + d23y*d23y + d23z*d23z;

	b1 += d01x*d12x + d01y*d12y + d01z*d12z;
	
	b2 += d01x*d23x + d01y*d23y + d01z*d23z;
    }

    a21 = a12;
    detA = a11*a22 - a12*a21;

    *alpha = (b1*a22 - b2*a12) / detA;
    *beta =  (b2*a11 - b1*a21) / detA;

#if VERBOSE
    printf("\n alpha: %8.5f  beta:%8.5f a11:%8.5e a12:%8.5e a21:%8.5e a22:%8.5e b1:%8.5e b2:%8.5e detA:%8.5e", *alpha, *beta, a11, a12, a21, a22, b1, b2, detA);


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
	printf ("\n Extrapolated position for ion %d is %10.7f  %10.7f  %10.7f", ion,
		iptr->ocrds1[0] + *alpha *(iptr->ocrds1[0] - iptr->ocrds2[0]) + *beta * (iptr->ocrds2[0] - iptr->ocrds3[0]), 
		iptr->ocrds1[1] + *alpha *(iptr->ocrds1[1] - iptr->ocrds2[1]) + *beta * (iptr->ocrds2[1] - iptr->ocrds3[1]), 
		iptr->ocrds1[2] + *alpha *(iptr->ocrds1[2] - iptr->ocrds2[2]) + *beta * (iptr->ocrds2[2] - iptr->ocrds3[2]));
	printf("\n Previous position0:  %10.7f  %10.7f  %10.7f", iptr->crds[0], iptr->crds[1], iptr->crds[2]);
	printf("\n Previous position1:  %10.7f  %10.7f  %10.7f", iptr->ocrds1[0], iptr->ocrds1[1], iptr->ocrds1[2]);
	printf("\n Previous position2:  %10.7f  %10.7f  %10.7f", iptr->ocrds2[0], iptr->ocrds2[1], iptr->ocrds2[2]);
	printf("\n Previous position3:  %10.7f  %10.7f  %10.7f", iptr->ocrds3[0], iptr->ocrds3[1], iptr->ocrds3[2]);
	printf("\n Previous position:   d1: %5.5e %5.5e %5.5e d2: %5.5e %5.5e %5.5e", 
		iptr->ocrds1[0] - iptr->ocrds2[0], iptr->ocrds1[1] - iptr->ocrds2[1], iptr->ocrds1[2] - iptr->ocrds2[2],
		iptr->ocrds2[0] - iptr->ocrds3[0], iptr->ocrds2[1] - iptr->ocrds3[1], iptr->ocrds2[2] - iptr->ocrds3[2]);

	printf("\n 1st order extrapolated position (%8.5e): %10.7f  %10.7f  %10.7f", b1/a11, 
		iptr->ocrds1[0] + b1/a11 * (iptr->ocrds1[0] - iptr->ocrds2[0]),
		iptr->ocrds1[1] + b1/a11 * (iptr->ocrds1[1] - iptr->ocrds2[1]),
		iptr->ocrds1[2] + b1/a11 * (iptr->ocrds1[2] - iptr->ocrds2[2]));
    }
#endif 


}
