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

// Mixes projectors for a single orbital
// using a power function.
void mix_betaxpsi1 (STATE *sp)
{

    int size, koffset, loffset, idx1, kpt, istate;
    rmg_double_t scale;

    kpt = sp->kidx;
    istate = sp->istate;

    scale = pow(1.0 - ct.prjmix, (rmg_double_t)istate);
    if(istate == 0) scale = 1.0 - ct.prjmix;
    size = ct.max_nl;
    koffset = kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl;
      
    for (idx1 = 0; idx1 < pct.num_nonloc_ions; idx1++) {

        /* For localized <beta|psi>, there is offset due to both k-point and ion*/
        loffset = koffset + idx1 * ct.num_states * ct.max_nl;
        my_scal( scale, &pct.oldsintR_local[loffset + istate * ct.max_nl], size);
        my_axpy(1.0 - scale, &pct.newsintR_local[loffset + istate * ct.max_nl],
                &pct.oldsintR_local[loffset + istate * ct.max_nl], size);


#if !GAMMA_PT

        my_scal( scale, &pct.oldsintI_local[loffset + istate * ct.max_nl], size);
        my_axpy(1.0 - scale, &pct.newsintI_local[loffset + istate * ct.max_nl],
                &pct.oldsintI_local[loffset + istate * ct.max_nl], size);

#endif

    }
}


