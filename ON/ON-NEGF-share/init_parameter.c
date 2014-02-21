/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/init_parameter.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void init_parameter()
 *   Basic initialization of some parameters.
 * INPUTS
 *
 * OUTPUT
 *   
 * PARENTS
 *   init.c
 * CHILDREN
 *
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main_on.h"

void init_parameter(STATE * states)
{

    int kpt, kst1, state, ion, st1;
    ION *iptr;
    rmg_double_t v1, v2, v3;
    int ispin;

    ct.psi_nbasis = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();
    ct.psi_nxgrid = get_NX_GRID();
    ct.psi_nygrid = get_NY_GRID();
    ct.psi_nzgrid = get_NZ_GRID();


    /* Set hartree boundary condition stuff */
    ct.vh_pxgrid = get_FPX0_GRID();
    ct.vh_pygrid = get_FPY0_GRID();
    ct.vh_pzgrid = get_FPZ0_GRID();

    ct.vh_pbasis = ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid;
    my_malloc_init( ct.vh_ext, ct.vh_pbasis, rmg_double_t );


    ct.vh_nbasis = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();
    ct.vh_nxgrid = get_FNX_GRID();
    ct.vh_nygrid = get_FNY_GRID();
    ct.vh_nzgrid = get_FNZ_GRID();
    /* Initialize some k-point stuff */
    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        v1 = TWO * PI * ct.kp[kpt].kpt[0] / ct.xside;
        v2 = TWO * PI * ct.kp[kpt].kpt[1] / ct.yside;
        v3 = TWO * PI * ct.kp[kpt].kpt[2] / ct.zside;

        ct.kp[kpt].kvec[0] = v1;
        ct.kp[kpt].kvec[1] = v2;
        ct.kp[kpt].kvec[2] = v3;

        ct.kp[kpt].kmag = v1 * v1 + v2 * v2 + v3 * v3;

        ct.kp[kpt].kstate = &states[kpt * ct.num_states];
        ct.kp[kpt + ct.spin * ct.num_states].kstate =
            &states[kpt * ct.num_states + ct.spin * ct.num_kpts * ct.num_states];
        ct.kp[kpt].kidx = kpt;

    }                           /* end for */

    /* Count up the total number of electrons */
    ct.ionic_charge = ZERO;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        ct.ionic_charge += ct.sp[iptr->species].zvalence;

    }                           /* end for */

    ct.nel = ct.ionic_charge + ct.background_charge;


    ct.hmaxgrid = ct.xside * ct.hxgrid;
    if (ct.yside * ct.hygrid > ct.hmaxgrid)
        ct.hmaxgrid = ct.yside * ct.hygrid;
    if (ct.zside * ct.hzgrid > ct.hmaxgrid)
        ct.hmaxgrid = ct.zside * ct.hzgrid;

    ct.hmingrid = ct.xside * ct.hxgrid;
    if (ct.yside * ct.hygrid < ct.hmingrid)
        ct.hmingrid = ct.yside * ct.hygrid;
    if (ct.zside * ct.hzgrid < ct.hmingrid)
        ct.hmingrid = ct.zside * ct.hzgrid;

    ct.anisotropy = ct.hmaxgrid / ct.hmingrid;

    if (ct.anisotropy > 1.1)
        error_handler(" Anisotropy too large");

    /* Set discretization array */
    ct.xcstart = ZERO;
    ct.ycstart = ZERO;
    ct.zcstart = ZERO;


    /* Some multigrid parameters */
    ct.poi_parm.sb_step = 1.0;
    ct.eig_parm.sb_step = 1.0;
    double t1, t2, part_occ;
    int full_occ;
    t1 = modf(ct.nel/2.0, &t2);
    full_occ = (int)(t2);
     
    part_occ = ct.nel - full_occ *2.0;

    for (ispin = 0; ispin <= ct.spin; ispin++)
    {
        for (kpt = pct.kstart; kpt < pct.kend; kpt++)
        {
            for (st1 = 0; st1 < ct.num_states; st1++)
            {
                kst1 = ispin * ct.num_states * ct.num_kpts + kpt * ct.num_states + st1;
                states[kst1].kidx = kpt;
                states[kst1].istate = st1;
                states[kst1].firstflag = 0;
                states[kst1].occupation[0] = 0.0;
                if(st1 < full_occ)  states[kst1].occupation[0] = 2.0;
                if(st1 == full_occ)  states[kst1].occupation[0] = part_occ;
            }                   /* end for */

        }                       /* end for */
    }                           /* end for ispin */



}                               /* end init */

/******/
