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

    int idx, nion, ion, stop, ip, ipindex, alloc, index_global, index_local;
    int istate, goffset, loffset;
    int id1, incx = 1, *pidx;
    REAL *nlarrayR, *nlarrayI, *sintR, *sintR_global;
    REAL *weiptr, *psiR;
    ION *iptr;
    SPECIES *sp;
    STATE *st;
#if !GAMMA_PT
    REAL *pR, *pI, *sintI,  *sintI_global, *psiI;
#endif


    alloc = P0_BASIS;

    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;
    my_calloc (nlarrayR, 2 * alloc, REAL);
    nlarrayI = nlarrayR + alloc;

    my_calloc (sintR_global, ct.num_ions * ct.num_states * ct.max_nl, REAL);
#if !GAMMA_PT
    my_calloc (sintI_global, ct.num_ions * ct.num_states * ct.max_nl, REAL);
#endif


    /*Zero whole array first */
    for (idx = 0; idx < ct.num_ions * ct.num_states * ct.max_nl; idx++)
    {
        sintR_global[idx] = 0.0;
#if !GAMMA_PT
        sintI_global[idx] = 0.0;
#endif
    }




    /* Loop over ions on this processor */
    for (nion = 0; nion < pct.num_nonloc_ions; nion++)
    {
        
        /*Actual index of the ion under consideration*/
        ion = pct.nonloc_ions_list[nion];

        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];
        stop = pct.idxptrlen[ion];

        if (stop)
        {
        
            sintR = &sintR_global[ion * ct.num_states * ct.max_nl];
#if !GAMMA_PT
            sintI = &sintI_global[ion * ct.num_states * ct.max_nl];
#endif

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

                for (ip = 0; ip < sp->nh; ip++)
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


    }                           /*end for (istate = 0; istate < ct.num_states; istate++) */



    /*Reset pointers to the beginning */
    sintR = sintR_global;
#if !GAMMA_PT
    sintI = sintI_global;
#endif

    /* Sum the sint array over all processors */
    id1 = ct.num_states * ct.num_ions * ct.max_nl;
    global_sums (sintR, &id1, pct.grid_comm);
#if !GAMMA_PT
    global_sums (sintI, &id1, pct.grid_comm);
#endif



    /*Pick up only the data relevant for this processor*/

    /*Offset due to a kpoint*/
    goffset = kpt * ct.num_ions * ct.num_states * ct.max_nl;
    loffset = kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl;

    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {
	index_global = goffset + pct.nonloc_ions_list[ion] * ct.num_states * ct.max_nl;
	index_local = loffset + ion *  ct.num_states * ct.max_nl;

	my_copy(&sintR_global[index_global], &pct.newsintR_local[index_local], ct.num_states * ct.max_nl);
#if !GAMMA_PT
	my_copy(&sintI_global[index_global], &pct.newsintI_local[index_local], ct.num_states * ct.max_nl);
#endif

    }

    my_free (nlarrayR);
    my_free (sintR_global);
#if !GAMMA_PT
    my_free (sintI_global);
#endif

}
