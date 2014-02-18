/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
get_qnmpsi:

modified by Wenchang Lu, 09-24-2005

 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "main_on.h"


/*shuchun wang*/
void get_dnmpsi(STATE *sp, double *kbpsi_one_state, double *work)
{
    int ion, ipindex, idx, ip1, ip2;
    double *prjptr, time1, time2;
    int PROJECTOR_SPACE;
    int ion2, nh;
    double *qqq, *kbpsi_one_ion, *qnm_weight;
    double *prj_sum;
    double qsum, qsum_all;


    /*
     *  dnm_weght[i] = sum_j dnm[i, j] *<beta_j|phi> 
     * prj_sum(r) = sum_i nm_weight[i] *beta_i(r)  
     */

    my_malloc_init( qnm_weight, ct.max_nl, REAL );
    my_malloc_init( prj_sum, ct.max_nlpoints, REAL );


    /* Zero out the work array */
    for (idx = 0; idx < sp->size; idx++)
        work[idx] = 0.;

    prjptr = projectors;

    for (ion2 = 0; ion2 < pct.n_ion_center; ion2++)
    {
        ion = pct.ionidx[ion2];
        ipindex = ion2 * ct.max_nl;
        qqq = pct.dnmI[ion];
        kbpsi_one_ion = &kbpsi_one_state[ipindex];
        nh = pct.prj_per_ion[ion];

        /* calculate qnm * <kb|phi> */

        qsum_all = 0.0;
        for (ip1 = 0; ip1 < nh; ip1++)
        {

            qsum = 0.0;

            for (ip2 = 0; ip2 < nh; ip2++)
            {

                qsum += qqq[ip1 * nh + ip2] * kbpsi_one_ion[ip2];
            }

            qnm_weight[ip1] = qsum;
            qsum_all += qsum;

        }

        if (fabs(qsum_all) > 0.0)
        {

            for (idx = 0; idx < ct.max_nlpoints; idx++)
            {
                prj_sum[idx] = 0.0;
            }


            for (ip1 = 0; ip1 < nh; ip1++)
            {
                for (idx = 0; idx < ct.max_nlpoints; idx++)
                {

                    prj_sum[idx] += qnm_weight[ip1] * prjptr[ip1 * ct.max_nlpoints + idx];
                }
            }

            /*
             *  project the prj_sum to the orbital |phi>  and stored in work 
             */



            qnm_beta_betapsi(sp, ion, prj_sum, work);

        }                       /* end if(fabs(...  */


        prjptr += pct.prj_per_ion[ion] * ct.max_nlpoints;       /*shuchu wang */

    }                           /* end for ion */

    my_free(qnm_weight);
    my_free(prj_sum);


}
