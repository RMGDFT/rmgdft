/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"
#include "common_prototypes.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* This function is used to normalize the wavefunctions */
void norm_psi1 (STATE * sp, int istate, int kidx)
{

    int idx, ion, nh, i, j, size, incx = 1, inh, sidx, sidx_local, nidx, oion;
    rmg_double_t sumbeta, sumpsi, sum, t1;
    rmg_double_t *tmp_psiR, *tmp_psiI, *qqq, *sintR, *sintI, *ptr;
    ION *iptr;

    sidx = kidx * ct.num_ions * ct.num_states * ct.max_nl + istate * ct.max_nl;
    sidx_local = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + istate * ct.max_nl;

    size = sp->pbasis;

    tmp_psiR = sp->psiR;
#if !GAMMA_PT
    tmp_psiI = sp->psiI;
#endif

    sumpsi = 0.0;
    sumbeta = 0.0;

    nidx = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {

	oion = pct.owned_ions_list[ion];
        
        iptr = &ct.ions[oion];
       
        nh = ct.sp[iptr->species].nh;
	
	/* Figure out index of owned ion in nonloc_ions_list array*/
	do {
	    
	    nidx++;
	    if (nidx >= pct.num_nonloc_ions)
		error_handler("Could not find matching entry in pct.nonloc_ions_list for owned ion %d", oion);
	
	} while (pct.nonloc_ions_list[nidx] != oion);

        qqq = pct.qqq[oion];


	/*nidx adds offset due to current ion*/
        sintR = &pct.newsintR_local[sidx_local + nidx * ct.num_states * ct.max_nl];
#if !GAMMA_PT
        sintI = &pct.newsintI_local[sidx_local+ nidx * ct.num_states * ct.max_nl];
#endif


        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            for (j = 0; j < nh; j++)
            {
                if (qqq[inh + j] != 0.0)
                {
#if GAMMA_PT
                    sumbeta += qqq[inh + j] * sintR[i] * sintR[j];
#else
                    sumbeta += qqq[inh + j] * (sintR[i] * sintR[j] + sintI[i] * sintI[j]);
                    sumbeta += qqq[inh + j] * (sintR[i] * sintI[j] - sintI[i] * sintR[j]);
#endif

                }
            }
        }
    }


    for (idx = 0; idx < get_P0_BASIS(); idx++)
    {
        sumpsi += tmp_psiR[idx] * tmp_psiR[idx];
#if !GAMMA_PT
        sumpsi += tmp_psiI[idx] * tmp_psiI[idx];
#endif
    }

    sum = real_sum_all (ct.vel * sumpsi + sumbeta, pct.grid_comm);
    sum = 1.0 / sum;
    if (sum < 0.0)
    {
        printf ("the %dth state is wrong\n", istate);
        error_handler ("<psi|S|psi> cann't be negative");
    }

    t1 = sqrt (sum);
    QMD_dscal (size, t1, tmp_psiR, incx);
#if !GAMMA_PT
    QMD_dscal (size, t1, tmp_psiI, incx);
#endif


    /* update <beta|psi> - Local version*/
    
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {

	ptr = &pct.newsintR_local[ion * ct.num_states * ct.max_nl];
        ptr += sidx_local;
        
	QMD_dscal (ct.max_nl, t1, ptr, incx);
#if !GAMMA_PT
	ptr = &pct.newsintI_local[ion * ct.num_states * ct.max_nl];
        ptr += sidx_local;
        QMD_dscal (ct.max_nl, t1, ptr, incx);
#endif
    }


}                               /* end norm_psi1 */
