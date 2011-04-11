/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/get_ke.c *****
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
 *   REAL get_ke(STATE *sp, int tid)
 *   Computes the kinetic energy of a given orbital. 
 * INPUTS
 *   sp: points to orbital structure (see main.h)
 *   tid: thread ID
 * OUTPUT
 *   kinetic energy is returned
 * PARENTS
 *   quench.c
 * CHILDREN
 *   gather_psi.c pack_ptos.c app6_del2.c
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"



#if GAMMA_PT
REAL get_ke (STATE * sp, int tid)
{

    int pbasis, sbasis;
    REAL *work1, *work2;
    REAL *tmp_psi, KE;
    int dimx, dimy, dimz;

#if 1
    REAL time1, time2;
    time1 = my_crtc ();
#endif


    dimx = sp->dimx;
    dimy = sp->dimy;
    dimz = sp->dimz;
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;


    /* Grab some memory */
    my_malloc (tmp_psi, 2 * (sbasis), REAL);
    work2 = tmp_psi + sbasis;

    gather_psi (tmp_psi, NULL, sp, tid);

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psi, tmp_psi, dimx, dimy, dimz);

    app6_del2 (tmp_psi, (P0_GRID *) work2);

    KE = -0.5 * ct.vel * QMD_sdot (pbasis, tmp_psi, 1, work2, 1);
    KE = real_sum_all (KE, pct.grid_comm);

    /* Release our memory */
    my_free (tmp_psi);

#if 1
    time2 = my_crtc ();
    rmg_timings (MG_EIGTIME, (time2 - time1), tid);
#endif

    return KE;

}                               /* end get_ke */

#else


/* Complex version */
REAL get_ke (STATE * sp, int tid)
{

    int pbasis, sbasis;
    REAL *work1, *work2;
    REAL *tmp_psiR, *tmp_psiI, KE;
    int dimx, dimy, dimz;

#if 1
    REAL time1, time2;
    time1 = my_crtc ();
#endif


    dimx = sp->dimx;
    dimy = sp->dimy;
    dimz = sp->dimz;
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;


    /* Grab some memory */
    my_malloc (tmp_psiR, 2 * (sbasis), REAL);
    work1 = tmp_psiR + sbasis;

    my_malloc (tmp_psiI, 2 * (sbasis), REAL);
    work2 = tmp_psiI + sbasis;


    gather_psi (tmp_psiR, tmp_psiI, sp, tid);

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psiR, tmp_psiR, dimx, dimy, dimz);

    //pack_ptos (sg_psiI, tmp_psiI, dimx, dimy, dimz);

    app6_del2 (tmp_psiR, (P0_GRID *) work1);
    app6_del2 (tmp_psiI, (P0_GRID *) work2);

    KE = -0.5 * ct.vel * QMD_sdot (pbasis, tmp_psiR, 1, work1, 1);
    KE += -0.5 * ct.vel * QMD_sdot (pbasis, tmp_psiI, 1, work2, 1);
    KE = real_sum_all (KE, pct.grid_comm);

    /* Release our memory */
    my_free (tmp_psiR);
    my_free (tmp_psiI);

#if 1
    time2 = my_crtc ();
    rmg_timings (MG_EIGTIME, (time2 - time1), tid);
#endif

    return KE;


}                               /* end get_ke */


#endif




/******/
