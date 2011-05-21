/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include <float.h>

void pack_rho_ctof (REAL * rho, REAL * rho_f)
{
    int i, j, k, ii, jj, kk, basis1, basis2, basis3;
    int pbasis, sbasis, dimx, dimy, dimz, in, jn, kn;
    REAL tmp1, tmp2, tmp3, frac, cc[10][4];
    REAL sum_rho, sum_rhof, coef;
    SS0_GRID rho_c;
    FP0_GRID *tpr_rho_f;

    int num = 0;

    dimx = PX0_GRID;
    dimy = PY0_GRID;
    dimz = PZ0_GRID;
    pbasis = P0_BASIS;
    sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);

    sum_rho = 0.0;
    sum_rhof = 0.0;
    for (i = 0; i < pbasis; i++)
        sum_rho += rho[i];
    sum_rho = real_sum_all (sum_rho, pct.grid_comm);
    sum_rho *= ct.vel;

    for (i = 0; i < FG_NX; i++)
    {
        frac = (REAL) i / (REAL) FG_NX;
        cc[i][0] = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
        cc[i][1] = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
        cc[i][2] = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
        cc[i][3] = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;
    }

    trade_imagesx (rho, &rho_c.b[0][0][0], dimx, dimy, dimz, 2);

    tpr_rho_f = (FP0_GRID *) rho_f;

    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {
                tpr_rho_f->s1.b[FG_NX * (i - 2)][FG_NX * (j - 2)][FG_NX * (k - 2)] =
                    rho_c.b[i][j][k];
                ++num;
            }
        }
    }

    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_NX; in++)
                {
                    tmp1 = 0.0;
                    tmp2 = 0.0;
                    tmp3 = 0.0;
                    basis1 = -1;

                    for (ii = 0; ii < 4; ii++)
                    {
                        tmp1 += cc[in][ii] * rho_c.b[i + basis1][j][k];
                        tmp2 += cc[in][ii] * rho_c.b[i][j + basis1][k];
                        tmp3 += cc[in][ii] * rho_c.b[i][j][k + basis1];
                        ++basis1;
                    }

                    tpr_rho_f->s1.b[FG_NX * (i - 2) + in][FG_NX * (j - 2)][FG_NX * (k - 2)] = tmp1;
                    tpr_rho_f->s1.b[FG_NX * (i - 2)][FG_NX * (j - 2) + in][FG_NX * (k - 2)] = tmp2;
                    tpr_rho_f->s1.b[FG_NX * (i - 2)][FG_NX * (j - 2)][FG_NX * (k - 2) + in] = tmp3;
                    num += 3;
                }

            }
        }
    }

    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_NX; in++)
                {
                    for (jn = 1; jn < FG_NX; jn++)
                    {

                        tmp1 = 0.0;
                        tmp2 = 0.0;
                        tmp3 = 0.0;
                        basis1 = -1;
                        for (ii = 0; ii < 4; ii++)
                        {
                            basis2 = -1;
                            for (jj = 0; jj < 4; jj++)
                            {
                                tmp1 +=
                                    cc[in][ii] * cc[jn][jj] * rho_c.b[i + basis1][j + basis2][k];
                                tmp2 +=
                                    cc[in][ii] * cc[jn][jj] * rho_c.b[i + basis1][j][k + basis2];
                                tmp3 +=
                                    cc[in][ii] * cc[jn][jj] * rho_c.b[i][j + basis1][k + basis2];
                                ++basis2;
                            }
                            ++basis1;
                        }
                        tpr_rho_f->s1.b[FG_NX * (i - 2) + in][FG_NX * (j - 2) +
                                                              jn][FG_NX * (k - 2)] = tmp1;
                        tpr_rho_f->s1.b[FG_NX * (i - 2) + in][FG_NX * (j - 2)][FG_NX * (k - 2) +
                                                                               jn] = tmp2;
                        tpr_rho_f->s1.b[FG_NX * (i - 2)][FG_NX * (j - 2) + in][FG_NX * (k - 2) +
                                                                               jn] = tmp3;
                        num += 3;
                    }
                }
            }
        }
    }


    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {

                for (in = 1; in < FG_NX; in++)
                {
                    for (jn = 1; jn < FG_NX; jn++)
                    {
                        for (kn = 1; kn < FG_NX; kn++)
                        {

                            tmp1 = 0.0;
                            basis1 = -1;
                            for (ii = 0; ii < 4; ii++)
                            {
                                basis2 = -1;
                                for (jj = 0; jj < 4; jj++)
                                {
                                    basis3 = -1;
                                    for (kk = 0; kk < 4; kk++)
                                    {
                                        tmp1 +=
                                            cc[in][ii] * cc[jn][jj] * cc[kn][kk] * rho_c.b[i +
                                                                                           basis1][j
                                                                                                   +
                                                                                                   basis2]
                                            [k + basis3];
                                        ++basis3;
                                    }
                                    ++basis2;
                                }
                                ++basis1;
                            }
                            tpr_rho_f->s1.b[FG_NX * (i - 2) + in][FG_NX * (j - 2) +
                                                                  jn][FG_NX * (k - 2) + kn] = tmp1;
                            num += 1;

                        }
                    }
                }
            }
        }
    }

    if (num != FP0_BASIS)
        error_handler ("there is something wrong here");

    for (i = 0; i < num; i++)
        sum_rhof += rho_f[i];
    sum_rhof = real_sum_all (sum_rhof, pct.grid_comm);
    sum_rhof *= ct.vel_f;
    coef = sum_rho / sum_rhof;
    for (i = 0; i < num; i++)
        rho_f[i] *= coef;

}
