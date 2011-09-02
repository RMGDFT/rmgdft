/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "md.h"


void precond_mg(double *res, double *work1, double *work2, int istate)
{
    int cycles, ione = 1, nits, stopp0;
    double diag, t1, time1, time2, d1, one = 1.;
    STATE *sp;

    int idx;
    REAL tem, tem1, tem2;
    int ion, ixx, iyy, izz;
    /* Pre and post smoothings on each level */


    int eig_pre[6] = { 0, 2, 2, 2, 2, 2 };
    int eig_post[6] = { 0, 2, 2, 2, 2, 2 };


    time1 = my_crtc();

    diag = -1. / ct.Ac;
#if FD4
    diag *= 1.0;
#endif
    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    /* Get state pointer */
    sp = &states[istate];

    ion = state_to_ion[istate];
    ixx = states[istate].ixmax - states[istate].ixmin + 1;
    iyy = states[istate].iymax - states[istate].iymin + 1;
    izz = states[istate].izmax - states[istate].izmin + 1;

    stopp0 = ixx * iyy * izz;




    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {


        pack_ptos(sg_orbit, work1, ixx, iyy, izz);

        app_cil_orbital(sg_orbit, work2, ixx, iyy, izz, ct.hxgrid, ct.hygrid, ct.hzgrid);

        saxpy(&stopp0, &one, res, &ione, work2, &ione);

        /*app_mask(istate, work2, 0); */
        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {

            /* Pack the residual data into multigrid array */
            pack_ptos(sg_orbit_res, work2, ixx, iyy, izz);


            /* Do multigrid step with solution in sg_twovpsi */


            mgrid_solv(sg_orbit, sg_orbit_res, work2, ixx, iyy, izz,
                       ct.hxgrid, ct.hygrid, ct.hzgrid, 0, pct.neighbors,
                       ct.eig_parm.levels, eig_pre, eig_post, 1, istate, &sp->inum, 1);


            pack_stop(sg_orbit, work2, ixx, iyy, izz);


            t1 = -1.;

        }
        else
        {

            t1 = ct.eig_parm.gl_step * diag;    /*shuchun wang */

        }                       /* end if cycles == ct.eig_parm.gl_pre */

        /* Update correction for wavefuntion */
        saxpy(&stopp0, &t1, work2, &ione, work1, &ione);

        app_mask(istate, work1, 0);

    }                           /* end for Smoothing cycles */




    time2 = my_crtc();
    d1 = time2 - time1;
    rmg_timings(PRECOND_TIME, d1);


}
