/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include "main.h"
#include <float.h>
#include <stdlib.h>
#include <math.h>

/* Calculates coefficients that will be used in orthogonalizing wavefunctions*/

/*Loop over ions is parallelized in this function only certain parts of newsintR
 * are updated when wavefunctions are updated, so newsintR (projectors applied on wavefunctions)
 * should be recalculated (using betaxpsi) once new_psi is finished. This is done 
 * at the end of ortho_full.c*/

#if GAMMA_PT
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL * cR, REAL * cI)
{
    int ion, i, j, inh, sidx1, sidx2;
    int idx, nh;
    REAL *tmp_psi2, *tmp_psi1;
    REAL *sint1R, *sint2R, *qqq;
    REAL sumpsi, sumbeta, sri;
    ION *iptr;


    sidx1 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1 = sp1->psiR;
    tmp_psi2 = sp2->psiR;

    sumpsi = 0.0;
    sumbeta = 0.0;

    /*This loop should be the same for all ions - qqq and newsintR should have 
     * the same values everywhere*/
    /*This is parallelized over ions */
    for (ion = pct.gridpe; ion < ct.num_ions; ion += NPES)
    {
        qqq = pct.qqq[ion];
        nh = pct.prj_per_ion[ion];
        iptr = &ct.ions[ion];

        /* get<beta|psi1> and <beta|psi2> */
        sint1R = &iptr->newsintR[sidx1];
        sint2R = &iptr->newsintR[sidx2];

        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            sri = sint1R[i];

            for (j = 0; j < nh; j++)
                sumbeta += qqq[inh + j] * sri * sint2R[j];
        }                       /*end for i */

    }                           /*end for ion */


    for (idx = 0; idx < P0_BASIS; idx++)
        sumpsi += tmp_psi2[idx] * tmp_psi1[idx];



    *cR = -1.0 * (ct.vel * sumpsi + sumbeta);

}



#else
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, REAL * cR, REAL * cI)
{
    int ion, i, j;
    int idx, nh;
    int sidx1, sidx2, inh;
    REAL *tmp_psi2R, *tmp_psi2I, *tmp_psi1R, *tmp_psi1I;
    REAL *sint1R, *sint1I, *sint2R, *sint2I, *qqq;
    REAL sumpsiR, sumpsiI, sumbetaR, sumbetaI, sri, sii;
    ION *iptr;


    sidx1 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * ct.num_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1R = sp1->psiR;
    tmp_psi1I = sp1->psiI;
    tmp_psi2R = sp2->psiR;
    tmp_psi2I = sp2->psiI;

    sumpsiR = 0.0;
    sumpsiI = 0.0;
    sumbetaR = 0.0;
    sumbetaI = 0.0;

    /*This is parallelized over ions */
    for (ion = pct.gridpe; ion < ct.num_ions; ion += NPES)
    {
        qqq = pct.qqq[ion];
        nh = pct.prj_per_ion[ion];
        iptr = &ct.ions[ion];

        /* get<beta|psi1> and <beta|psi2> */
        sint1R = &iptr->newsintR[sidx1];
        sint1I = &iptr->newsintI[sidx1];
        sint2R = &iptr->newsintR[sidx2];
        sint2I = &iptr->newsintI[sidx2];


        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            sri = sint1R[i];
            sii = sint1I[i];

            for (j = 0; j < nh; j++)
            {
                sumbetaR += qqq[inh + j] * (sri * sint2R[j] + sii * sint2I[j]);
                sumbetaI += qqq[inh + j] * (sri * sint2I[j] - sii * sint2R[j]);
            }                   /*end for j */
        }                       /*end for i */
    }                           /*end for ion */

    for (idx = 0; idx < P0_BASIS; idx++)
    {
        sumpsiR += (tmp_psi2R[idx] * tmp_psi1R[idx] + tmp_psi2I[idx] * tmp_psi1I[idx]);
        sumpsiI += (tmp_psi2I[idx] * tmp_psi1R[idx] - tmp_psi2R[idx] * tmp_psi1I[idx]);
    }



    *cR = ct.vel * sumpsiR + sumbetaR;
    *cI = ct.vel * sumpsiI + sumbetaI;

}

#endif
