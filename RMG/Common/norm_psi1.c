/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* This function is used to normalize the wavefunctions */
void norm_psi1 (STATE * sp, int istate, int kidx)
{

    int idx, ion, nh, i, j, size, incx = 1, inh, sidx, sidx_local;
    REAL sumbetaR, sumpsi, sum, t1;
    REAL *tmp_psiR, *tmp_psiI, *qqq, *sintR, *sintI, *ptr;
    ION *iptr;
#if !GAMMA_PT
    REAL sumbetaI;
#endif

    sidx = kidx * ct.num_ions * ct.num_states * ct.max_nl + istate * ct.max_nl;

    size = sp->pbasis;

    tmp_psiR = sp->psiR;
#if !GAMMA_PT
    tmp_psiI = sp->psiI;
#endif

    sumpsi = 0.0;
    sumbetaR = 0.0;
#if !GAMMA_PT
    sumbetaI = 0.0;
#endif

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        qqq = pct.qqq[ion];
        nh = pct.prj_per_ion[ion];
        iptr = &ct.ions[ion];


        sintR = &iptr->newsintR[sidx];
#if !GAMMA_PT
        sintI = &iptr->newsintI[sidx];
#endif


        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            for (j = 0; j < nh; j++)
            {
                if (qqq[inh + j] != 0.0)
                {
#if GAMMA_PT
                    sumbetaR += qqq[inh + j] * sintR[i] * sintR[j];
#else
                    sumbetaR += qqq[inh + j] * (sintR[i] * sintR[j] + sintI[i] * sintI[j]);
                    sumbetaI += qqq[inh + j] * (sintR[i] * sintI[j] - sintI[i] * sintR[j]);
#endif

                }
            }
        }
    }

#if !GAMMA_PT
    if (fabs (sumbetaI / sumbetaR) > 1.0e-10)
    {
        printf ("<psi|s|psi> = %e + i %e\n", sumbetaR, sumbetaI);
        error_handler ("<psi|s|psi> can't be complex number");
    }
#endif

    for (idx = 0; idx < P0_BASIS; idx++)
    {
        sumpsi += tmp_psiR[idx] * tmp_psiR[idx];
#if !GAMMA_PT
        sumpsi += tmp_psiI[idx] * tmp_psiI[idx];
#endif

    }

    sumpsi = real_sum_all (sumpsi, pct.grid_comm);
    sum = sumpsi * ct.vel + sumbetaR;
    sum = 1.0 / sum;
    if (sum < 0.0)
    {
        printf ("the %dth state is wrong\n", istate);
        error_handler ("<psi|S|psi> cann't be negative");
    }

    t1 = sqrt (sum);
    QMD_sscal (size, t1, tmp_psiR, incx);
#if !GAMMA_PT
    QMD_sscal (size, t1, tmp_psiI, incx);
#endif


    /* update <beta|psi> */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];

        ptr = &iptr->newsintR[sidx];
        QMD_sscal (ct.max_nl, t1, ptr, incx);
#if !GAMMA_PT
        ptr = &iptr->newsintI[sidx];
        QMD_sscal (ct.max_nl, t1, ptr, incx);
#endif
    }
    
    
    /* update <beta|psi> - Local version*/
    sidx_local = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + istate * ct.max_nl;
    
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

	ptr = &pct.newsintR_local[ion * ct.num_states * ct.max_nl];
        ptr += sidx_local;
        
	QMD_sscal (ct.max_nl, t1, ptr, incx);
#if !GAMMA_PT
	ptr = &pct.newsintI_local[ion * ct.num_states * ct.max_nl];
        ptr += sidx_local;
        QMD_sscal (ct.max_nl, t1, ptr, incx);
#endif
    }


}                               /* end norm_psi1 */
