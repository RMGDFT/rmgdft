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
#include "main.h"
#include "prototypes_on.h"

void init_parameter(STATE * states)
{

    int kpt, kst1, ion, st1;
    ION *iptr;
    double v1, v2, v3;

    ct.psi_nbasis = get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();


    /* Set hartree boundary condition stuff */
    ct.vh_pxgrid = get_FPX0_GRID();
    ct.vh_pygrid = get_FPY0_GRID();
    ct.vh_pzgrid = get_FPZ0_GRID();

    ct.vh_pbasis = ct.vh_pxgrid * ct.vh_pygrid * ct.vh_pzgrid;
    my_malloc_init( ct.vh_ext, ct.vh_pbasis, double );


    ct.vh_nbasis = get_FNX_GRID() * get_FNY_GRID() * get_FNZ_GRID();
    /* Initialize some k-point stuff */


    ct.is_gamma = true;

    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        v1 = TWO * PI * ct.kp[kpt].kpt[0] / get_xside();
        v2 = TWO * PI * ct.kp[kpt].kpt[1] / get_yside();
        v3 = TWO * PI * ct.kp[kpt].kpt[2] / get_zside();

        ct.kp[kpt].kvec[0] = v1;
        ct.kp[kpt].kvec[1] = v2;
        ct.kp[kpt].kvec[2] = v3;

        ct.kp[kpt].kpt[0] *= TWO * PI;
        ct.kp[kpt].kpt[1] *= TWO * PI;
        ct.kp[kpt].kpt[2] *= TWO * PI;

        ct.kp[kpt].kmag = v1 * v1 + v2 * v2 + v3 * v3;

        if(ct.kp[kpt].kmag != 0.0) ct.is_gamma = false;

        ct.kp[kpt].kstate = &states[kpt * ct.num_states];
        ct.kp[kpt + (int)ct.spin_polarization * ct.num_states].kstate =
            &states[kpt * ct.num_states + (int)ct.spin_polarization * ct.num_kpts * ct.num_states];
        ct.kp[kpt].kidx = kpt;

    }                           /* end for */

    /* Count up the total number of electrons */
    ct.ionic_charge = ZERO;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &Atoms[ion];
        ct.ionic_charge += Species[iptr->species].zvalence;

    }                           /* end for */

    ct.nel = ct.ionic_charge + ct.background_charge;


    ct.hmaxgrid = get_xside() * get_hxgrid();
    if (get_yside() * get_hygrid() > ct.hmaxgrid)
        ct.hmaxgrid = get_yside() * get_hygrid();
    if (get_zside() * get_hzgrid() > ct.hmaxgrid)
        ct.hmaxgrid = get_zside() * get_hzgrid();

    ct.hmingrid = get_xside() * get_hxgrid();
    if (get_yside() * get_hygrid() < ct.hmingrid)
        ct.hmingrid = get_yside() * get_hygrid();
    if (get_zside() * get_hzgrid() < ct.hmingrid)
        ct.hmingrid = get_zside() * get_hzgrid();



    if (get_anisotropy() > 1.1)
    {
        dprintf("\n ct.hmaxgrid  %f %f ", ct.hmaxgrid, ct.hmingrid);
        error_handler(" Anisotropy too large");
    }
    /* Set discretization array */
    ct.xcstart = ZERO;
    ct.ycstart = ZERO;
    ct.zcstart = ZERO;


    /* Some multigrid parameters */
    ct.poi_parm.sb_step = 1.0;
    ct.eig_parm.sb_step = 1.0;
    double t2, part_occ;
    int full_occ;
    modf(ct.nel/2.0, &t2);
    full_occ = (int)(t2);
     
    part_occ = ct.nel - full_occ *2.0;

    for (kpt = 0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            kst1 = kpt * ct.num_states + st1;
            states[kst1].kidx = kpt;
            states[kst1].istate = st1;
            states[kst1].firstflag = 0;
            states[kst1].occupation[0] = 0.0;
            if(st1 < full_occ)  states[kst1].occupation[0] = 2.0/(1.0 + ct.spin_flag);
            if(st1 == full_occ)  states[kst1].occupation[0] = part_occ/(1.0+ ct.spin_flag);
            states[kst1].occupation[1] = 0.0;
            if(st1 < full_occ)  states[kst1].occupation[1] = 2.0/(1.0 + ct.spin_flag);
            if(st1 == full_occ)  states[kst1].occupation[1] = part_occ/(1.0+ ct.spin_flag);
        }                   /* end for */

    }                       /* end for */



}                               /* end init */

/******/
