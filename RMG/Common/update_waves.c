/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include "main.h"
#include <float.h>
#include <stdlib.h>
#include <math.h>


/* Updates wavefunction of state sp2 so that it is orthogonal to the one from state sp1 
 * (at least that is what I think it is doing)
 * Also newsint (beta|psi) is updated*/

/*Loop over ions is parallelized in this function only certain parts of newsintR
 * are updated when wavefunctions are updated, so newsintR (projectors applied on wavefunctions)
 * should be recalculated (using betaxpsi) once new_psi is finished. This is done 
 * at the end of ortho_full.c*/

#if GAMMA_PT
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t cR, rmg_double_t cI)
{
    int ion, sidx1, sidx2, lsidx1, lsidx2;
    int size, incx = 1;
    rmg_double_t *tmp_psi2, *tmp_psi1, *ptr1, *ptr2;
    rmg_double_t t1;
    ION *iptr;

    size =pct.P0_BASIS;

    sidx1 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;
    
    lsidx1 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    lsidx2 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1 = sp1->psiR;
    tmp_psi2 = sp2->psiR;


    t1 = cR;


// For FAST_ORTHO we update in a block in the main loop
#if !FAST_ORTHO
    /*update the wavefunction psi2 */
    QMD_daxpy (size, t1, tmp_psi1, incx, tmp_psi2, incx);
#endif

    /* update localized  <beta|psi2> */
    for (ion = 0; ion < pct.num_nonloc_ions; ion ++)
    {

        ptr1 = &pct.newsintR_local[lsidx1 + ion * ct.num_states * ct.max_nl];
        ptr2 = &pct.newsintR_local[lsidx2 + ion * ct.num_states * ct.max_nl];
        QMD_daxpy (ct.max_nl, t1, ptr1, incx, ptr2, incx);
    }

}
#else
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t cR, rmg_double_t cI)
{
    int ion;
    int incx = 1;
    int sidx1, sidx2, lsidx1, lsidx2, idx;
    rmg_double_t *tmp_psi2R, *tmp_psi2I, *tmp_psi1R, *tmp_psi1I, *ptr1R, *ptr1I, *ptr2R, *ptr2I;
    rmg_double_t sumpsiR, sumpsiI, sumbetaR, sumbetaI, sri, sii;
    ION *iptr;


    sidx1 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;
    
    lsidx1 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    lsidx2 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1R = sp1->psiR;
    tmp_psi1I = sp1->psiI;
    tmp_psi2R = sp2->psiR;
    tmp_psi2I = sp2->psiI;



    /*update the wavefunction psi2 */
    for (idx = 0; idx <pct.P0_BASIS; idx++)
    {
        tmp_psi2R[idx] += -cR * tmp_psi1R[idx] + cI * tmp_psi1I[idx];
        tmp_psi2I[idx] += -cR * tmp_psi1I[idx] - cI * tmp_psi1R[idx];
    }

    /* update localized <beta|psi2> */
    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {


        ptr1R = &pct.newsintR_local[lsidx1 + ion * ct.num_states * ct.max_nl];
        ptr1I = &pct.newsintI_local[lsidx1 + ion * ct.num_states * ct.max_nl];
        ptr2R = &pct.newsintR_local[lsidx2 + ion * ct.num_states * ct.max_nl];
        ptr2I = &pct.newsintI_local[lsidx2 + ion * ct.num_states * ct.max_nl];

        QMD_daxpy (ct.max_nl, -cR, ptr1R, incx, ptr2R, incx);
        QMD_daxpy (ct.max_nl, cI, ptr1I, incx, ptr2R, incx);

        QMD_daxpy (ct.max_nl, -cR, ptr1I, incx, ptr2I, incx);
        QMD_daxpy (ct.max_nl, -cI, ptr1R, incx, ptr2I, incx);
    }

}
#endif
