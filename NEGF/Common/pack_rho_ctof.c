/************************** SVN Revision Information **************************
 **    $Id: pack_rho_ctof.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include "md.h"
#include <float.h>

void pack_rho_ctof (REAL * rho_c, REAL * rho_f)
{
    int i, j, k, ii, jj, kk, basis1, basis2, basis3;
    int dimx, dimy, dimz, in, jn, kn;
    REAL tmp1, tmp2, tmp3, frac, cc[6][4];
    REAL sum_rho, sum_rhof, coef;
    REAL *ssrho_c;
    int ixs, iys, ix1, iy1;

    int num = 0;

    dimx = PX0_GRID;
    dimy = PY0_GRID;
    dimz = PZ0_GRID;

    ixs = (dimy + 4) * (dimz + 4);
    iys = (dimz + 4);
    ix1 = FPZ0_GRID * FPY0_GRID;
    iy1 = FPZ0_GRID;

    my_malloc_init( ssrho_c, (dimx + 4) * (dimy + 4) * (dimz + 4), REAL );

    sum_rho = 0.0;
    sum_rhof = 0.0;
    for (i = 0; i < P0_BASIS; i++)
        sum_rho += rho_c[i];
    sum_rho = real_sum_all (sum_rho);
    sum_rho *= ct.vel;

    for (i = 0; i < RHO_NX; i++)
    {
        frac = (REAL) i / (REAL) RHO_NX;
        cc[i][0] = -frac * (1.0 - frac) * (2.0 - frac) / 6.0;
        cc[i][1] = (1.0 + frac) * (1.0 - frac) * (2.0 - frac) / 2.0;
        cc[i][2] = (1.0 + frac) * frac * (2.0 - frac) / 2.0;
        cc[i][3] = -(1.0 + frac) * frac * (1.0 - frac) / 6.0;
    }

    trade_images2 (rho_c, ssrho_c, dimx, dimy, dimz);


    for (i = 2; i < PX0_GRID + 2; i++)
    {
        for (j = 2; j < PY0_GRID + 2; j++)
        {
            for (k = 2; k < PZ0_GRID + 2; k++)
            {
                rho_f[(RHO_NX * (i - 2)) * ix1 + (RHO_NX * (j - 2)) * iy1 + (RHO_NX * (k - 2))] =
                    ssrho_c[i * ixs + j * iys + k];
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

                for (in = 1; in < RHO_NX; in++)
                {
                    tmp1 = 0.0;
                    tmp2 = 0.0;
                    tmp3 = 0.0;
                    basis1 = -1;

                    for (ii = 0; ii < 4; ii++)
                    {
                        tmp1 += cc[in][ii] * ssrho_c[(i + basis1) * ixs + j * iys + k];
                        tmp2 += cc[in][ii] * ssrho_c[i * ixs + (j + basis1) * iys + k];
                        tmp3 += cc[in][ii] * ssrho_c[i * ixs + j * iys + (k + basis1)];
                        ++basis1;
                    }

                    rho_f[(RHO_NX * (i - 2) + in) * ix1 + (RHO_NX * (j - 2)) * iy1 +
                          (RHO_NX * (k - 2))] = tmp1;
                    rho_f[(RHO_NX * (i - 2)) * ix1 + (RHO_NX * (j - 2) + in) * iy1 +
                          (RHO_NX * (k - 2))] = tmp2;
                    rho_f[(RHO_NX * (i - 2)) * ix1 + (RHO_NX * (j - 2)) * iy1 +
                          (RHO_NX * (k - 2) + in)] = tmp3;
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

                for (in = 1; in < RHO_NX; in++)
                {
                    for (jn = 1; jn < RHO_NX; jn++)
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
                                    cc[in][ii] * cc[jn][jj] * ssrho_c[(i + basis1) * ixs +
                                                                      (j + basis2) * iys + k];
                                tmp2 +=
                                    cc[in][ii] * cc[jn][jj] * ssrho_c[(i + basis1) * ixs + j * iys +
                                                                      (k + basis2)];
                                tmp3 +=
                                    cc[in][ii] * cc[jn][jj] * ssrho_c[i * ixs + (j + basis1) * iys +
                                                                      (k + basis2)];
                                ++basis2;
                            }
                            ++basis1;
                        }
                        rho_f[(RHO_NX * (i - 2) + in) * ix1 + (RHO_NX * (j - 2) + jn) * iy1 +
                              (RHO_NX * (k - 2))] = tmp1;
                        rho_f[(RHO_NX * (i - 2) + in) * ix1 + (RHO_NX * (j - 2)) * iy1 +
                              (RHO_NX * (k - 2) + jn)] = tmp2;
                        rho_f[(RHO_NX * (i - 2)) * ix1 + (RHO_NX * (j - 2) + in) * iy1 +
                              (RHO_NX * (k - 2) + jn)] = tmp3;
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

                for (in = 1; in < RHO_NX; in++)
                {
                    for (jn = 1; jn < RHO_NX; jn++)
                    {
                        for (kn = 1; kn < RHO_NX; kn++)
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
                                            cc[in][ii] * cc[jn][jj] * cc[kn][kk] *
                                            ssrho_c[(i + basis1) * ixs + (j + basis2) * iys +
                                                    (k + basis3)];
                                        ++basis3;
                                    }
                                    ++basis2;
                                }
                                ++basis1;
                            }
                            rho_f[(RHO_NX * (i - 2) + in) * ix1 + (RHO_NX * (j - 2) + jn) * iy1 +
                                  (RHO_NX * (k - 2) + kn)] = tmp1;
                            num += 1;

                        }
                    }
                }
            }
        }
    }

    if (num != FP0_BASIS)
        error_handler ("there is something wrong here");
    my_free(ssrho_c);

    for (i = 0; i < num; i++)
        sum_rhof += rho_f[i];
    sum_rhof = real_sum_all (sum_rhof);
    sum_rhof *= ct.vel_f;
    /*printf("\n\nsum_rho=%f   sum_rhof=%f \n\n",sum_rho,sum_rhof); */
    coef = sum_rho / sum_rhof;
    for (i = 0; i < num; i++)
        rho_f[i] *= coef;

    coef = 1.0 / (REAL) (RHO_NX * RHO_NY * RHO_NZ);
    for (i = 0; i < num; i++)
        rho_f[i] *= coef;

}
