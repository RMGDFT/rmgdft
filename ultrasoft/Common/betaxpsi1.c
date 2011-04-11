/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


void betaxpsi1 (STATE * states, int kpt)
{

    int idx, ion, stop, ip, ipindex, alloc, index_global, index_local;
    int istate;
    int id1, incx = 1, *pidx;
    REAL *nlarrayR, *nlarrayI, *sintR;
    REAL *weiptr, *psiR;
    ION *iptr;
    STATE *st;
#if !GAMMA_PT
    REAL *pR, *pI, *sintI, *psiI;
#endif


    alloc = P0_BASIS;

    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;
    my_calloc (nlarrayR, 2 * alloc, REAL);
    nlarrayI = nlarrayR + alloc;


    sintR = &ct.ions[0].newsintR[kpt * ct.num_ions * ct.num_states * ct.max_nl];
#if !GAMMA_PT
    sintI = &ct.ions[0].newsintI[kpt * ct.num_ions * ct.num_states * ct.max_nl];
#endif

    /*Zero whole array first */
    for (idx = 0; idx < ct.num_ions * ct.num_states * ct.max_nl; idx++)
    {
        sintR[idx] = 0.0;
#if !GAMMA_PT
        sintI[idx] = 0.0;
#endif
    }




    /* Loop over ions on this processor */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];
        stop = pct.idxptrlen[ion];

        if (stop)
        {

            pidx = pct.nlindex[ion];
#if !GAMMA_PT
            pR = pct.phaseptr[ion];
            pR += 2 * kpt * stop;
            pI = pR + stop;
#endif


            for (istate = 0; istate < ct.num_states; istate++)
            {

                st = &states[istate];
                psiR = st->psiR;
#if !GAMMA_PT
                psiI = st->psiI;
#endif


#if GAMMA_PT
                /* Copy wavefunction into temporary array */
                for (idx = 0; idx < stop; idx++)
                    nlarrayR[idx] = psiR[pidx[idx]];
#else
                for (idx = 0; idx < stop; idx++)
                    nlarrayR[idx] = psiR[pidx[idx]] * pR[idx] - psiI[pidx[idx]] * pI[idx];

                for (idx = 0; idx < stop; idx++)
                    nlarrayI[idx] = psiI[pidx[idx]] * pR[idx] + psiR[pidx[idx]] * pI[idx];
#endif

                /* <Beta|psi>                                       */

                weiptr = pct.weight[ion];
                ipindex = istate * ct.max_nl;

                for (ip = 0; ip < pct.prj_per_ion[ion]; ip++)
                {

                    sintR[ipindex] = ct.vel * sdot (&stop, nlarrayR, &incx, weiptr, &incx);
#if !GAMMA_PT
                    sintI[ipindex] = ct.vel * sdot (&stop, nlarrayI, &incx, weiptr, &incx);
#endif

                    weiptr += pct.idxptrlen[ion];
                    ipindex++;

                }



            }

        }

        /*Advance pointers to the next ion */
        sintR += ct.num_states * ct.max_nl;
#if !GAMMA_PT
        sintI += ct.num_states * ct.max_nl;
#endif

    }                           /*end for (istate = 0; istate < ct.num_states; istate++) */



    /*Reset pointers to the beginning */
    sintR = &ct.ions[0].newsintR[kpt * ct.num_ions * ct.num_states * ct.max_nl];
#if !GAMMA_PT
    sintI = &ct.ions[0].newsintI[kpt * ct.num_ions * ct.num_states * ct.max_nl];
#endif

    /* Sum the sint array over all processors */
    id1 = ct.num_states * ct.num_ions * ct.max_nl;
    global_sums (sintR, &id1, pct.grid_comm);
#if !GAMMA_PT
    global_sums (sintI, &id1, pct.grid_comm);
#endif



    /*Pick up only the data relevant for this processor*/
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {
	index_global = kpt * pct.nonloc_ions_list[ion] * ct.num_states * ct.max_nl;
	index_local = kpt * ion *  ct.num_states * ct.max_nl;

	my_copy(&sintR[index_global], &pct.newsintR_local[index_local], ct.num_states * ct.max_nl);

    }


    my_free (nlarrayR);

}
