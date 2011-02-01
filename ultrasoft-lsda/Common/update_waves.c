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
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL cR, REAL cI)
{
    int ion, sidx1, sidx2;
    int size, incx = 1;
    REAL *tmp_psi2, *tmp_psi1, *ptr1, *ptr2;
    REAL t1;
    ION *iptr;

    size = P0_BASIS;

    sidx1 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1 = sp1->psiR;
    tmp_psi2 = sp2->psiR;


    t1 = cR;



    /*update the wavefunction psi2 */
    QMD_saxpy (size, t1, tmp_psi1, incx, tmp_psi2, incx);

    /* update <beta|psi2> */
    /*Parallelization over ions */
    for (ion = pct.thispe; ion < ct.num_ions; ion += NPES)
    {

        iptr = &ct.ions[ion];

        ptr1 = &iptr->newsintR[sidx1];
        ptr2 = &iptr->newsintR[sidx2];
        QMD_saxpy (ct.max_nl, t1, ptr1, incx, ptr2, incx);
    }

}
#else
void update_waves (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL cR, REAL cI)
{
    int ion;
    int incx = 1;
    int sidx1, sidx2, idx;
    REAL *tmp_psi2R, *tmp_psi2I, *tmp_psi1R, *tmp_psi1I, *ptr1R, *ptr1I, *ptr2R, *ptr2I;
    REAL sumpsiR, sumpsiI, sumbetaR, sumbetaI, sri, sii;
    ION *iptr;


    sidx1 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1R = sp1->psiR;
    tmp_psi1I = sp1->psiI;
    tmp_psi2R = sp2->psiR;
    tmp_psi2I = sp2->psiI;



    /*update the wavefunction psi2 */
    for (idx = 0; idx < P0_BASIS; idx++)
    {
        tmp_psi2R[idx] += -cR * tmp_psi1R[idx] + cI * tmp_psi1I[idx];
        tmp_psi2I[idx] += -cR * tmp_psi1I[idx] - cI * tmp_psi1R[idx];
    }


    /* update <beta|psi2> */
    /*Parallelization over ions */
    for (ion = pct.thispe; ion < ct.num_ions; ion += NPES)
    {

        iptr = &ct.ions[ion];

        ptr1R = &iptr->newsintR[sidx1];
        ptr1I = &iptr->newsintI[sidx1];
        ptr2R = &iptr->newsintR[sidx2];
        ptr2I = &iptr->newsintI[sidx2];

        QMD_saxpy (ct.max_nl, -cR, ptr1R, incx, ptr2R, incx);
        QMD_saxpy (ct.max_nl, cI, ptr1I, incx, ptr2R, incx);

        QMD_saxpy (ct.max_nl, -cR, ptr1I, incx, ptr2I, incx);
        QMD_saxpy (ct.max_nl, -cI, ptr1R, incx, ptr2I, incx);
    }


}
#endif
