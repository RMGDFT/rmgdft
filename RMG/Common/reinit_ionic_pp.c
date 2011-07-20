/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/***** RMG/Common/reinit_ionic_pp.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2010  : Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void reinit_ionic_pp (STATE * states, REAL * vnuc, REAL * rhocore, REAL * rhoc)
 *   drive routine for updating systemic potentials due to motion of ions.
 * INPUTS
 *   states: all wave functions (see main.h)
 *   vnuc: pseudopotential
 *   rhocore: charge density of core electrons, only useful when we 
 *            include non-linear core correction for pseudopotential.
 *   rhoc:    compensating charge density
 * OUTPUT
 *   all the above inputs are updated
 * PARENTS
 *   fastrlx.c, moldyn.c
 * CHILDREN
 *   init_nuc,get_QI,get_nlop,get_weight,get_qqq,betaxpsi,mix_betaxpsi
 * SOURCE
 */

#include "main.h"



void reinit_ionic_pp (STATE * states, REAL * vnuc, REAL * rhocore, REAL * rhoc)
{

    /* Update items that change when the ionic coordinates change */
    init_nuc (vnuc, rhoc, rhocore);
    get_QI ();
    get_nlop ();

    /*Other things that need to be recalculated when ionic positions change */
    get_weight ();
    get_qqq ();

    if (!verify ("calculation_mode", "Band Structure Only"))
    {
        betaxpsi (states);
        mix_betaxpsi(0);
    }

}                               /* end reinit_ionic_pp */


/******/
