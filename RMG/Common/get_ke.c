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
 *   rmg_double_t get_ke(STATE *sp, int tid)
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
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"



#if GAMMA_PT
rmg_double_t get_ke (STATE * sp, int tid)
{

    int pbasis, sbasis;
    rmg_double_t *work2;
    rmg_double_t *tmp_psi, KE;
    int dimx, dimy, dimz;

#if 1
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif


    dimx = get_PX0_GRID();
    dimy = get_PY0_GRID();
    dimz = get_PZ0_GRID();
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;


    /* Grab some memory */
    my_malloc (tmp_psi, 2 * (sbasis), rmg_double_t);
    work2 = tmp_psi + sbasis;

    gather_psi (tmp_psi, NULL, sp, tid);

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psi, tmp_psi, dimx, dimy, dimz);

    app6_del2 (tmp_psi, work2, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_hxgrid(), get_hygrid(), get_hzgrid());

    KE = -0.5 * get_vel() * QMD_ddot (pbasis, tmp_psi, 1, work2, 1);
    KE = real_sum_all (KE, pct.grid_comm);

    /* Release our memory */
    my_free (tmp_psi);

#if 1
    time2 = my_crtc ();
#endif

    return KE;

}                               /* end get_ke */

#else


/* Complex version */
rmg_double_t get_ke (STATE * sp, int tid)
{

    int pbasis, sbasis;
    rmg_double_t *work1, *work2;
    rmg_double_t *tmp_psiR, *tmp_psiI, KE;
    int dimx, dimy, dimz;

#if 1
    rmg_double_t time1, time2;
    time1 = my_crtc ();
#endif


    dimx = sp->dimx;
    dimy = sp->dimy;
    dimz = sp->dimz;
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;


    /* Grab some memory */
    my_malloc (tmp_psiR, 2 * (sbasis), rmg_double_t);
    work1 = tmp_psiR + sbasis;

    my_malloc (tmp_psiI, 2 * (sbasis), rmg_double_t);
    work2 = tmp_psiI + sbasis;


    gather_psi (tmp_psiR, tmp_psiI, sp, tid);

    /* Pack psi into smoothing array */
    //pack_ptos (sg_psiR, tmp_psiR, dimx, dimy, dimz);

    //pack_ptos (sg_psiI, tmp_psiI, dimx, dimy, dimz);

    app6_del2 (tmp_psiR, (P0_GRID *) work1);
    app6_del2 (tmp_psiI, (P0_GRID *) work2);

    KE = -0.5 * get_vel() * QMD_ddot (pbasis, tmp_psiR, 1, work1, 1);
    KE += -0.5 * get_vel() * QMD_ddot (pbasis, tmp_psiI, 1, work2, 1);
    KE = real_sum_all (KE, pct.grid_comm);

    /* Release our memory */
    my_free (tmp_psiR);
    my_free (tmp_psiI);

#if 1
    time2 = my_crtc ();
#endif

    return KE;


}                               /* end get_ke */


#endif




/******/
