/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/app_nl.c *****
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
 * INPUTS
 * OUTPUT
 * PARENTS
 * CHILDREN
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void mix_betaxpsi (int mix)
{

    int size;
    REAL *newsintR, *oldsintR, *newsintI, *oldsintI;


      newsintR = ct.ions[0].newsintR;
      oldsintR = ct.ions[0].oldsintR;
#if !GAMMA_PT
      newsintI = ct.ions[0].newsintI;
      oldsintI = ct.ions[0].oldsintI;
#endif

      size = ct.num_kpts * ct.num_ions * ct.num_states * ct.max_nl;

    
    if (mix)
    {
	my_scal( 1.0 - ct.prjmix, oldsintR, size);
	my_axpy(ct.prjmix, newsintR, oldsintR, size); 

#if !GAMMA_PT
	my_scal( 1.0 - ct.prjmix, oldsintI, size);
	my_axpy(ct.prjmix, newsintI, oldsintI, size); 
#endif
    }


    else
    {
	my_copy (newsintR, oldsintR, ct.num_kpts * ct.num_ions * ct.num_states * ct.max_nl);
#if !GAMMA_PT
	my_copy (newsintI, oldsintI, ct.num_kpts * ct.num_ions * ct.num_states * ct.max_nl);
#endif
    }

    

}


