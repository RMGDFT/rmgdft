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
    double *newsintR, *oldsintR;

    /*Local version*/
    newsintR = pct.newsintR_local;
    oldsintR = pct.oldsintR_local;
    
    size = ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl;
    if(!ct.is_gamma) size *=2;
    
      
    if (mix)
    {
        my_scal( 1.0 - ct.prjmix, oldsintR, size);
        my_axpy(ct.prjmix, newsintR, oldsintR, size); 

    }
    else
    {
        my_copy (newsintR, oldsintR, size);
    }

}


