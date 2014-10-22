/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"


template void GetGamma<double> (double * gammaR, int ion, int nh, Kpoint<double> **Kptr);
template void GetGamma<std::complex<double>> (double * gammaR, int ion, int nh, Kpoint<std::complex<double>> **Kptr);


template <typename OrbitalType> void GetGamma (double * gammaR, int ion, int nh , Kpoint<OrbitalType> **Kptr)
{
    int i, j, idx, kidx, istate;
    rmg_double_t t1, sintNR, sintMR;
    rmg_double_t sintNI, sintMI;


    idx = 0;
    for (i = 0; i < nh; i++)
    {
        for (j = i; j < nh; j++)
        {
            gammaR[idx] = 0.0;
            for (kidx = 0; kidx < ct.num_kpts; kidx++)
            {
                for (istate = 0; istate < ct.num_states; istate++)
                {
                    t1 = Kptr[kidx]->Kstates[istate].occupation[0] * Kptr[kidx]->kweight;
                    sintNR = pct.newsintR_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                        istate * ct.max_nl + i];
                    sintMR = pct.newsintR_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                        istate * ct.max_nl + j];
                    if(ct.is_gamma)
                    {
                        gammaR[idx] += t1 * sintNR * sintMR;
                    }
                    else
                    {
                        sintNI = pct.newsintI_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                            istate * ct.max_nl + i];
                        sintMI = pct.newsintI_local[kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ion *  ct.num_states * ct.max_nl +
                            istate * ct.max_nl + j];
                        gammaR[idx] += t1 * (sintNR * sintMR + sintNI * sintMI);
                    }
                }               /*end for istate */
            }                   /*end for kidx */
            ++idx;
        }                       /*end for j */
    }                           /*end for i */
}
