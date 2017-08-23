/************************** SVN Revision Information **************************
 **    $Id: precond_mg.c 3394 2016-03-05 15:31:45Z ebriggs $    **
******************************************************************************/
 
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "transition.h"
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"

extern BaseGrid *OG;
FiniteDiff *MGFD;

void PrecondMg(double *psiR, double *work1, STATE *sp)
{
    int cycles, ione = 1, nits, stopp0;
    double diag, t1, one = 1.;
    double Zfac = 2.0 * ct.max_zvalence;
    int levels = ct.eig_parm.levels;
    if(ct.scf_steps < 2) levels = 0;


    int idx;
    int ixx, iyy, izz;
    /* Pre and post smoothings on each level */


    int eig_pre[6] = { 0, 6, 6, 6, 6, 6 };
    int eig_post[6] = { 0, 3, 3, 3, 3, 3 };

//    eig_pre[ct.eig_parm.levels] = 50;

    diag = -1. / ct.Ac;
#if FD4
    diag *= 1.0;
#endif
    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    ixx = sp->ixmax - sp->ixmin + 1;
    iyy = sp->iymax - sp->iymin + 1;
    izz = sp->izmax - sp->izmin + 1;

    double hxgrid = get_hxgrid();
    double hygrid = get_hygrid();
    double hzgrid = get_hzgrid();

    if(OG == NULL) OG = new BaseGrid(ixx, iyy, izz, 1, 1, 1, 0, 1);
    if(MGFD == NULL) MGFD = new FiniteDiff(&Rmg_L, OG, CLUSTER, CLUSTER, CLUSTER, 1, 2);
    stopp0 = ixx * iyy * izz;


    idx = (ixx + 8) * (iyy + 8) * (izz + 8);
    double *work2 = new double[4*idx];
    double *work3 = new double[4*idx];
    double *work4 = new double[4*idx];

    ZeroBoundary(work1, ixx, iyy, izz);

    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {


        diag = MGFD->app_del2_np(work1, work2, hxgrid, hygrid, hzgrid);
        daxpy(&stopp0, &one, psiR, &ione, work2, &ione);


        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {

            /* Pack the residual data into multigrid array */
            pack_ptos(work3, work2, ixx, iyy, izz);


            /* Do multigrid step with solution in sg_twovpsi */


            mgrid_solv_local(work4, work3, work2, ixx, iyy, izz,
                       hxgrid, hygrid, hzgrid, 0, get_neighbors(),
                       levels, eig_pre, eig_post, ct.eig_parm.sb_step,
                       1, sp->istate, &sp->inum, 1, Zfac);


            pack_stop(work4, work2, ixx, iyy, izz);
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
        ZeroBoundary(work1, ixx, iyy, izz);

    }                           /* end for Smoothing cycles */



    delete [] work4;
    delete [] work3;
    delete [] work2;


}
