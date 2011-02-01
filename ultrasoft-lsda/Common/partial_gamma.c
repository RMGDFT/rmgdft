/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "main.h"

void partial_gamma (int ion, REAL * par_gammaR, REAL * par_omegaR, ION * iptr, int nh,
                    REAL * newsintR_x, REAL * newsintR_y, REAL * newsintR_z,
                    REAL * newsintI_x, REAL * newsintI_y, REAL * newsintI_z)
{
    int i, j, idx, kidx, istate, size, index;
    REAL *gamma_x, *gamma_y, *gamma_z;
    REAL *omega_x, *omega_y, *omega_z;
    REAL t1, betaxpsiNR, betaxpsiMR;
#if !GAMMA_PT
    REAL betaxpsiNI, betaxpsiMI;
    REAL *psixbetaI_x, *psixbetaI_y, *psixbetaI_z;
#endif
    REAL *psixbetaR_x, *psixbetaR_y, *psixbetaR_z;
    STATE *sta;


    size = nh * (nh + 1) / 2;
    gamma_x = par_gammaR;
    gamma_y = gamma_x + size;
    gamma_z = gamma_y + size;

    omega_x = par_omegaR;
    omega_y = omega_x + size;
    omega_z = omega_y + size;

    for (idx = 0; idx < size; idx++)
    {
        gamma_x[idx] = 0.0;
        gamma_y[idx] = 0.0;
        gamma_z[idx] = 0.0;
        omega_x[idx] = 0.0;
        omega_y[idx] = 0.0;
        omega_z[idx] = 0.0;
    }


    for (kidx = 0; kidx < ct.num_kpts; kidx++)
    {
        sta = ct.kp[kidx].kstate;
        for (istate = 0; istate < ct.num_states; istate++)
        {
            t1 = sta->occupation * ct.kp[kidx].kweight;


            index =
                ion * ct.num_kpts * ct.num_states * ct.max_nl + kidx * ct.num_states * ct.max_nl +
                istate * ct.max_nl;

            psixbetaR_x = &newsintR_x[index];
            psixbetaR_y = &newsintR_y[index];
            psixbetaR_z = &newsintR_z[index];

#if !GAMMA_PT
            psixbetaI_x = &newsintI_x[index];
            psixbetaI_y = &newsintI_y[index];
            psixbetaI_z = &newsintI_z[index];
#endif


            idx = 0;
            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    betaxpsiNR = iptr->newsintR[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                                istate * ct.max_nl + i];
                    betaxpsiMR = iptr->newsintR[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                                istate * ct.max_nl + j];


#if GAMMA_PT
                    gamma_x[idx] +=
                        t1 * (psixbetaR_x[i] * betaxpsiMR + betaxpsiNR * psixbetaR_x[j]);
                    gamma_y[idx] +=
                        t1 * (psixbetaR_y[i] * betaxpsiMR + betaxpsiNR * psixbetaR_y[j]);
                    gamma_z[idx] +=
                        t1 * (psixbetaR_z[i] * betaxpsiMR + betaxpsiNR * psixbetaR_z[j]);
                    omega_x[idx] +=
                        t1 * sta->eig * (psixbetaR_x[i] * betaxpsiMR + betaxpsiNR * psixbetaR_x[j]);
                    omega_y[idx] +=
                        t1 * sta->eig * (psixbetaR_y[i] * betaxpsiMR + betaxpsiNR * psixbetaR_y[j]);
                    omega_z[idx] +=
                        t1 * sta->eig * (psixbetaR_z[i] * betaxpsiMR + betaxpsiNR * psixbetaR_z[j]);
#else
                    betaxpsiNI = iptr->newsintI[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                                istate * ct.max_nl + i];
                    betaxpsiMI = iptr->newsintI[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                                istate * ct.max_nl + j];



                    gamma_x[idx] += t1 * (psixbetaR_x[i] * betaxpsiMR + psixbetaI_x[i] * betaxpsiMI
                                          + betaxpsiNR * psixbetaR_x[j] +
                                          betaxpsiNI * psixbetaI_x[j]);
                    gamma_y[idx] +=
                        t1 * (psixbetaR_y[i] * betaxpsiMR + psixbetaI_y[i] * betaxpsiMI +
                              betaxpsiNR * psixbetaR_y[j] + betaxpsiNI * psixbetaI_y[j]);
                    gamma_z[idx] +=
                        t1 * (psixbetaR_z[i] * betaxpsiMR + psixbetaI_z[i] * betaxpsiMI +
                              betaxpsiNR * psixbetaR_z[j] + betaxpsiNI * psixbetaI_z[j]);

                    omega_x[idx] +=
                        t1 * sta->eig * (psixbetaR_x[i] * betaxpsiMR + psixbetaI_x[i] * betaxpsiMI +
                                         betaxpsiNR * psixbetaR_x[j] + betaxpsiNI * psixbetaI_x[j]);
                    omega_y[idx] +=
                        t1 * sta->eig * (psixbetaR_y[i] * betaxpsiMR + psixbetaI_y[i] * betaxpsiMI +
                                         betaxpsiNR * psixbetaR_y[j] + betaxpsiNI * psixbetaI_y[j]);
                    omega_z[idx] +=
                        t1 * sta->eig * (psixbetaR_z[i] * betaxpsiMR + psixbetaI_z[i] * betaxpsiMI +
                                         betaxpsiNR * psixbetaR_z[j] + betaxpsiNI * psixbetaI_z[j]);
#endif
                    ++idx;
                }               /* end for j */
            }                   /*end for i */
            sta++;
        }                       /*end for istate */
    }                           /*end for kidx */

}
