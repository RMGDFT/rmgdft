/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"

void get_gamma (REAL * gammaR, ION * iptr, int nh)
{
    int i, j, idx, kidx, istate;
    REAL t1, sintNR, sintMR;
#if !GAMMA_PT
    REAL sintNI, sintMI;
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
                    t1 = sta->occupation * ct.kp[kidx].kweight;
                    sintNR = iptr->newsintR[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + i];
                    sintMR = iptr->newsintR[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + j];
#if GAMMA_PT
                    gammaR[idx] += t1 * sintNR * sintMR;
#else
                    sintNI = iptr->newsintI[kidx * ct.num_ions * ct.num_states * ct.max_nl +
                                            istate * ct.max_nl + i];
                    sintMI = iptr->newsintI[kidx * ct.num_ions * ct.num_states * ct.max_nl +
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
