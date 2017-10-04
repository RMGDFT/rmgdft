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
    double step = ct.eig_parm.gl_step;
    int levels = ct.eig_parm.levels;
    if(ct.scf_steps < 2) levels = 0;


    int idx;
    /* Pre and post smoothings on each level */
    int eig_pre[MAX_MG_LEVELS] = { 8, 8, 8, 20, 20, 20, 20, 20 };
    int eig_post[MAX_MG_LEVELS] = { 2, 2, 2, 2, 2, 2, 2, 2 };

    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;

    int ixx = sp->ixmax - sp->ixmin + 1;
    int iyy = sp->iymax - sp->iymin + 1;
    int izz = sp->izmax - sp->izmin + 1;
    int dx2 = ixx / 2 + 1;
    int dy2 = iyy / 2 + 1;
    int dz2 = izz / 2 + 1;
    int ixoff = 0;
    int iyoff = 0;
    int izoff = 0;

    double hxgrid = get_hxgrid();
    double hygrid = get_hygrid();
    double hzgrid = get_hzgrid();

    if(OG == NULL) OG = new BaseGrid(ixx, iyy, izz, 1, 1, 1, 0, 1);
    if(MGFD == NULL) MGFD = new FiniteDiff(&Rmg_L, OG, CLUSTER, CLUSTER, CLUSTER, 1, 2);
    stopp0 = ixx * iyy * izz;

    int offsiz = ct.kohn_sham_fd_order + 2;
    idx = (ixx + offsiz) * (iyy + offsiz) * (izz + offsiz);
    offsiz /= 2;
    double *work2 = new double[4*idx];
    double *work3 = new double[4*idx];
    double *work4 = new double[4*idx];

//    ZeroBoundary(work1, ixx, iyy, izz);

    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {


        diag = MGFD->app_del2_np(work1, work2, hxgrid, hygrid, hzgrid);
        daxpy(&stopp0, &one, psiR, &ione, work2, &ione);


        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {

            /* Pack the residual data into multigrid array */
            ZeroBoundary(work2, ixx, iyy, izz, 1);
            mg_restrict(work2, work4, ixx-2, iyy-2, izz-2, dx2-2, dy2-2, dz2-2, ixoff, iyoff, izoff);
            ZeroBoundary(work4, dx2,dy2,dz2, 1);
        
            /* Do multigrid step with solution in sg_twovpsi */
            MgridSolvLocal(work2, work4, work3, dx2, dy2, dz2,
                       2.0*hxgrid, 2.0*hygrid, 2.0*hzgrid, 0, get_neighbors(),
                       levels, eig_pre, eig_post, step,
                       1, sp->istate, &sp->inum, 2.0*Zfac);

            ZeroBoundary(work2, dx2,dy2,dz2, 1);
            mg_prolong(work4, work2, ixx-2, iyy-2, izz-2, dx2-2, dy2-2, dz2-2, ixoff, iyoff, izoff);
            ZeroBoundary(work4, ixx, iyy, izz, 1);

            t1 = -1.;
            /* Update correction for wavefuntion */
            daxpy(&stopp0, &t1, work4, &ione, work1, &ione);

        }
        else
        {
            double t5 = diag - Zfac;
            t5 = -1.0 / t5;
            t1 = ct.eig_parm.gl_step * t5;
            /* Update correction for wavefuntion */
            daxpy(&stopp0, &t1, work2, &ione, work1, &ione);
        }                       /* end if cycles == ct.eig_parm.gl_pre */

        ZeroBoundary(work1, ixx, iyy, izz);

    }                           /* end for Smoothing cycles */



    delete [] work4;
    delete [] work3;
    delete [] work2;


}
