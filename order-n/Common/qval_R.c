/************************** SVN Revision Information **************************
 **    $Id: qval_R.c 587 2006-08-18 22:57:56Z miro $    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "md.h"

void qval_R(int ih, int jh, REAL r, REAL * x, REAL * qlig, REAL * drqlig, REAL invdr, int *nhtol,
            int *nhtom, int *indv, REAL * ylm, REAL * ylm_x, REAL * ylm_y, REAL * ylm_z,
            REAL ap[][9][9], int lpx[][9], int lpl[][9][9], REAL * Q_x, REAL * Q_y, REAL * Q_z,
            SPECIES * sp)
{
    int ivl, jvl;
    int nb, mb, nmb, lm, lp, l;
    REAL qrad, drqrad, *q_tpr, *drq_tpr;
    REAL x_hat, y_hat, z_hat;

    x_hat = x[0] / (r + 1.0e-9);
    y_hat = x[1] / (r + 1.0e-9);
    z_hat = x[2] / (r + 1.0e-9);

    nb = indv[ih];
    mb = indv[jh];
    if (nb < mb)
        nmb = mb * (mb + 1) / 2 + nb;
    else
        nmb = nb * (nb + 1) / 2 + mb;


    ivl = nhtol[ih] * nhtol[ih] + nhtom[ih];
    jvl = nhtol[jh] * nhtol[jh] + nhtom[jh];

    *Q_x = 0.0;
    *Q_y = 0.0;
    *Q_z = 0.0;

    for (lm = 1; lm <= lpx[ivl][jvl]; lm++)
    {
        lp = lpl[ivl][jvl][lm - 1];
        if (lp == 0)
            l = 0;
        else if ((lp >= 1) && (lp < 4))
            l = 1;
        else if ((lp >= 4) && (lp < 9))
            l = 2;
        else if ((lp >= 9) && (lp < 16))
            l = 3;
        else if ((lp >= 16) && (lp < 25))
            l = 4;
        else
            error_handler("L>4");

        q_tpr = qlig + (nmb * sp->nlc + l) * MAX_QLIG;
        drq_tpr = drqlig + (nmb * sp->nlc + l) * MAX_QLIG;
        qrad = linint(q_tpr, r, invdr);
        drqrad = linint(drq_tpr, r, invdr);

        *Q_x += (qrad * ap[lp][ivl][jvl] * ylm_x[lp] + drqrad * ap[lp][ivl][jvl] * ylm[lp] * x_hat);
        *Q_y += (qrad * ap[lp][ivl][jvl] * ylm_y[lp] + drqrad * ap[lp][ivl][jvl] * ylm[lp] * y_hat);
        *Q_z += (qrad * ap[lp][ivl][jvl] * ylm_z[lp] + drqrad * ap[lp][ivl][jvl] * ylm[lp] * z_hat);
/*if(r<1.0e-8) printf("r=%f  qrad=%f  drqrd=%f  lp=%d  ylm=%f ylm_x=%f  ylm_y=%f  ylm_z=%f\n",r,qrad,drqrad,lp,ylm[lp],ylm_x[lp],ylm_y[lp],ylm_z[lp]);*/
    }
}
