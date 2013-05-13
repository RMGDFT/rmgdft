/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"

void get_gamma (rmg_double_t * gammaR, int ion, int nh)
{
    int i, j, idx, kidx, istate;
    rmg_double_t t1, sintNR, sintMR;
#if !GAMMA_PT
    rmg_double_t sintNI, sintMI;
#endif
    STATE *sta;


    idx = 0;
    for (i = 0; i < nh; i++)
    {
        for (j = i; j < nh; j++)
        {
            gammaR[idx] = 0.0;
            for (kidx = 0; kidx < ct.num_kpts; kidx++)
            {
                sta = ct.kp[kidx].kstate;
                for (istate = 0; istate < ct.num_states; istate++)
                {
                    t1 = sta->occupation[0] * ct.kp[kidx].kweight;
                    sintNR = pct.newsintR_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + i];
                    sintMR = pct.newsintR_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + j];
#if GAMMA_PT
                    gammaR[idx] += t1 * sintNR * sintMR;
#else
                    sintNI = pct.newsintI_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + i];
                    sintMI = pct.newsintI_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + j];
                    gammaR[idx] += t1 * (sintNR * sintMR + sintNI * sintMI);
#endif
                    sta++;
                }               /*end for istate */
            }                   /*end for kidx */
            ++idx;
        }                       /*end for j */
    }                           /*end for i */
}
