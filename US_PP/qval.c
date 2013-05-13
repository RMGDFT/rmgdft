/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include <float.h>

rmg_double_t qval (int ih, int jh, rmg_double_t r, rmg_double_t invdr, rmg_double_t * ptpr, int *nhtol, int *nhtom,
           int *indv, rmg_double_t * ylm, rmg_double_t ap[][9][9], int lpx[][9], int lpl[][9][9], SPECIES * sp)
{
    int ivl, jvl;
    int nb, mb, nmb, lm, lp, l;
    rmg_double_t qrad, sum, *ptpr1;

    nb = indv[ih];
    mb = indv[jh];
    if (nb < mb)
        nmb = mb * (mb + 1) / 2 + nb;
    else
        nmb = nb * (nb + 1) / 2 + mb;

    ivl = nhtol[ih] * nhtol[ih] + nhtom[ih];
    jvl = nhtol[jh] * nhtol[jh] + nhtom[jh];

    sum = 0;
    for (lm = 0; lm < lpx[ivl][jvl]; lm++)
    {
        lp = lpl[ivl][jvl][lm];
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
            error_handler ("L>4");

        ptpr1 = ptpr + (nmb * sp->nlc + l) * MAX_QLIG;
        /*shuchun wang */
        qrad = linint (ptpr1, r, invdr);
        sum += qrad * ap[lp][ivl][jvl] * ylm[lp];
    }
    return (sum);
}
