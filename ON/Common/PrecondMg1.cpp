/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "transition.h"

void PrecondMg1(double *res, double *work1, double *work2, int istate)
{
    int cycles, ione = 1, nits, stopp0;
    double diag, t1, one = 1.;
    STATE *sp;
    double Zfac = 2.0 * ct.max_zvalence;
    int levels = ct.eig_parm.levels;
    FiniteDiff FD(&Rmg_L);
    if(ct.scf_steps < 2) levels = 0;

    int idx;
    int ixx, iyy, izz;
    /* Pre and post smoothings on each level */


    int eig_pre[6] = { 6, 6, 6, 6, 6, 6 };
    int eig_post[6] = { 3, 3, 3, 3, 3, 3 };

    double *work3;

//    eig_pre[ct.eig_parm.levels] = 50;
    double hxgrid = get_hxgrid();
    double hygrid = get_hygrid();
    double hzgrid = get_hzgrid();


    diag = -1. / ct.Ac;
#if FD4
    diag *= 1.0;
#endif
    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    /* Get state pointer */
    sp = &states[istate];

    ixx = states[istate].ixmax - states[istate].ixmin + 1;
    iyy = states[istate].iymax - states[istate].iymin + 1;
    izz = states[istate].izmax - states[istate].izmin + 1;

    stopp0 = ixx * iyy * izz;


    idx = (ixx + 8) * (iyy +8) * (izz+8);
    work3 = new double[idx];

    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {

        Rmg_T->trade_imagesx_central_local(work1, work3, ixx, iyy, izz, 1);
        //diag = FD.app_cil_sixth (work3, work2, ixx, iyy, izz, hxgrid, hygrid, hzgrid);
        diag = FD.app_cil_fourth (work3, work2, ixx, iyy, izz, hxgrid, hygrid, hzgrid);
        //diag = FD.app8_del2 (work3, work2, ixx, iyy, izz, hxgrid, hygrid, hzgrid);


        daxpy(&stopp0, &one, res, &ione, work2, &ione);
        /*app_mask(istate, work2, 0);*/

        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {

            /* Pack the residual data into multigrid array */
            //pack_ptos(sg_orbit_res, work2, ixx, iyy, izz);
            Rmg_T->trade_imagesx_central_local(sg_orbit_res, work3, ixx, iyy, izz, 1);
            CPP_app_smooth (work3, sg_orbit_res, ixx, iyy, izz);


            /* Do multigrid step with solution in sg_twovpsi */


//            mgrid_solv_local(sg_orbit, sg_orbit_res, work2, ixx, iyy, izz,
//                       get_hxgrid(), get_hygrid(), get_hzgrid(), 0, get_neighbors(),
//                       ct.eig_parm.levels, eig_pre, eig_post, ct.eig_parm.sb_step,
//                       1, istate, &sp->inum, 1, Zfac);
            MgridSolvLocal(sg_orbit, sg_orbit_res, work2, ixx, iyy, izz,
                       get_hxgrid(), get_hygrid(), get_hzgrid(), 0, get_neighbors(),
                       levels, eig_pre, eig_post, ct.eig_parm.sb_step,
                       1, sp->istate, &sp->inum, Zfac);


            pack_stop(sg_orbit, work2, ixx, iyy, izz);


            t1 = -1.;

        }
        else
        {
            double t5 = diag - Zfac;
            t5 = -1.0 / t5;
            t1 = ct.eig_parm.gl_step * t5;

        }                       /* end if cycles == ct.eig_parm.gl_pre */

        /* Update correction for wavefuntion */
        daxpy(&stopp0, &t1, work2, &ione, work1, &ione);

        app_mask(istate, work1, 0);

    }                           /* end for Smoothing cycles */

    delete [] work3;
}
