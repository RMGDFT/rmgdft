/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "main.h"
#include "common_prototypes.h"

void partial_betaxpsi (int ion, fftwnd_plan p2, rmg_double_t * newsintR_x, rmg_double_t * newsintR_y,
                       rmg_double_t * newsintR_z, rmg_double_t * newsintI_x, rmg_double_t * newsintI_y, rmg_double_t * newsintI_z,
                       ION * iptr)
{

    int idx, kidx, istate, size, nh, index;
    int alloc, prjcount, count;
    int incx = 1, *pidx, ip;
    rmg_double_t *workR;
    rmg_double_t *beta_x, *beta_y, *beta_z;
    rmg_double_t *temp_psiR;
    STATE *sta;
#if !GAMMA_PT
    rmg_double_t *workI, *pI, *temp_psiI, *pR;
#endif

    nh = ct.sp[ct.ions[ion].species].nh;
    alloc = get_P0_BASIS();
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;

    count = pct.idxptrlen[ion];
    pidx = pct.nlindex[ion];

    if (count)
    {
#if GAMMA_PT
        my_malloc (workR, alloc, rmg_double_t);
#else
        my_malloc (workR, 2 * alloc, rmg_double_t);
        workI = workR + alloc;
#endif
    }


    size = nh * pct.idxptrlen[ion];
    if (size)
    {
#if !FDIFF_BETA
        my_malloc (beta_x, 3 * size, rmg_double_t);
        beta_y = beta_x + size;
        beta_z = beta_y + size;

        for (idx = 0; idx < size; idx++)
        {
            beta_x[idx] = ZERO;
            beta_y[idx] = ZERO;
            beta_z[idx] = ZERO;
        }

        /*partial_beta(ion, beta_x, beta_y,beta_z, iptr, p1, p2); */
        get_derweight (ion, beta_x, beta_y, beta_z, iptr, p2);
#else

        beta_x = pct.weight_derx[ion];
        beta_y = pct.weight_dery[ion];
        beta_z = pct.weight_derz[ion];

#endif
    }
    else
    {
        beta_x = NULL;
        beta_y = NULL;
        beta_z = NULL;
    }

    for (kidx = 0; kidx < ct.num_kpts; kidx++)
    {

        sta = ct.kp[kidx].kstate;

#if !GAMMA_PT
        pR = pct.phaseptr[ion];
        pR += 2 * kidx * count;
        pI = pR + count;
#endif

        for (istate = 0; istate < ct.num_states; istate++)
        {

            if (count)
            {

                /*Gather_psi is not necessary, getting pointer should be enough
                   gather_psi(temp_psiR, temp_psiI, sta, 0); */
                temp_psiR = sta->psiR;
#if !GAMMA_PT
                temp_psiI = sta->psiI;
#endif


                for (idx = 0; idx < count; idx++)
                {
#if GAMMA_PT
                    workR[idx] = temp_psiR[pidx[idx]];
#else
                    workR[idx] = temp_psiR[pidx[idx]] * pR[idx] - temp_psiI[pidx[idx]] * pI[idx];
                    workI[idx] = temp_psiI[pidx[idx]] * pR[idx] + temp_psiR[pidx[idx]] * pI[idx];
#endif
                }

                prjcount = 0;
                index =
                    ion * ct.num_kpts * ct.num_states * ct.max_nl +
                    kidx * ct.num_states * ct.max_nl + istate * ct.max_nl;
                for (ip = 0; ip < nh; ip++)
                {

                    newsintR_x[index] =
                        get_vel() * QMD_ddot (count, workR, incx, &beta_x[prjcount], incx);
                    newsintR_y[index] =
                        get_vel() * QMD_ddot (count, workR, incx, &beta_y[prjcount], incx);
                    newsintR_z[index] =
                        get_vel() * QMD_ddot (count, workR, incx, &beta_z[prjcount], incx);
#if !GAMMA_PT
                    newsintI_x[index] =
                        get_vel() * QMD_ddot (count, workI, incx, &beta_x[prjcount], incx);
                    newsintI_y[index] =
                        get_vel() * QMD_ddot (count, workI, incx, &beta_y[prjcount], incx);
                    newsintI_z[index] =
                        get_vel() * QMD_ddot (count, workI, incx, &beta_z[prjcount], incx);
#endif
                    prjcount += count;
                    index++;

                }               /*end for ip */
            }                   /*end if count */
            sta++;
        }                       /* end for istate */
    }

    if (count)
        my_free (workR);
    /*my_free(temp_psiR); */
# if !FDIFF_BETA
    if (size)
        my_free (beta_x);
#endif
}
