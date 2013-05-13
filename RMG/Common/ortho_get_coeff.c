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

#if FAST_ORTHO
#if GAMMA_PT
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t * cR, rmg_double_t * cI, rmg_double_t *Oij)
{
    int ion, i, j, inh, sidx1, sidx2, nidx, oion;
    int idx, nh;
    rmg_double_t *tmp_psi2, *tmp_psi1;
    rmg_double_t *sint1R, *sint2R, *qqq;
    rmg_double_t sumpsi, sumbeta, sri;
    ION *iptr;
    SPECIES *sp;


    sidx1 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1 = sp1->psiR;
    tmp_psi2 = sp2->psiR;

    sumpsi = 0.0;
    sumbeta = 0.0;

    /*This loop should be the same for all ions - qqq and newsintR should have 
     * the same values everywhere*/
    /*This is parallelized over ions */
    nidx = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {
	oion = pct.owned_ions_list[ion];
        
        iptr = &ct.ions[oion];
        sp = &ct.sp[iptr->species];
       
        nh = sp->nh;
	
	/* Figure out index of owned ion in nonloc_ions_list array*/
	do {
	    
	    nidx++;
	    if (nidx >= pct.num_nonloc_ions)
		error_handler("Could not find matching entry in pct.nonloc_ions_list for owned ion %d", oion);
	
	} while (pct.nonloc_ions_list[nidx] != oion);
        
        qqq = pct.qqq[oion];

        /* get<beta|psi1> and <beta|psi2> */
        sint1R = &pct.newsintR_local[sidx1 + nidx * ct.num_states * ct.max_nl];
        sint2R = &pct.newsintR_local[sidx2 + nidx * ct.num_states * ct.max_nl];

        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            sri = sint1R[i];

            for (j = 0; j < nh; j++)
                sumbeta += qqq[inh + j] * sri * sint2R[j];
        }                       /*end for i */

    }                           /*end for ion */

    sumpsi = Oij[ist1 * ct.num_states + ist2];

    *cR = -1.0 * (ct.vel * sumpsi + sumbeta);

}



#else
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t * cR, rmg_double_t * cI)
{
    int ion, i, j, nidx, oion;
    int idx, nh;
    int sidx1, sidx2, inh;
    rmg_double_t *tmp_psi2R, *tmp_psi2I, *tmp_psi1R, *tmp_psi1I;
    rmg_double_t *sint1R, *sint1I, *sint2R, *sint2I, *qqq;
    rmg_double_t sumpsiR, sumpsiI, sumbetaR, sumbetaI, sri, sii;
    ION *iptr;
    SPECIES *sp;


    sidx1 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1R = sp1->psiR;
    tmp_psi1I = sp1->psiI;
    tmp_psi2R = sp2->psiR;
    tmp_psi2I = sp2->psiI;

    sumpsiR = 0.0;
    sumpsiI = 0.0;
    sumbetaR = 0.0;
    sumbetaI = 0.0;

    /*This is parallelized over ions */
    nidx = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {
	oion = pct.owned_ions_list[ion];
        
        iptr = &ct.ions[oion];
        sp = &ct.sp[iptr->species];
       
        nh = sp->nh;
	
	/* Figure out index of owned ion in nonloc_ions_list array*/
	do {
	    
	    nidx++;
	    if (nidx >= pct.num_nonloc_ions)
		error_handler("Could not find matching entry in pct.nonloc_ions_list for owned ion %d", oion);
	
	} while (pct.nonloc_ions_list[nidx] != oion);
        
        qqq = pct.qqq[oion];

        /* get<beta|psi1> and <beta|psi2> */
        sint1R = &pct.newsintR_local[sidx1 + nidx * ct.num_states * ct.max_nl];
        sint1I = &pct.newsintI_local[sidx1 + nidx * ct.num_states * ct.max_nl];
        sint2R = &pct.newsintR_local[sidx2 + nidx * ct.num_states * ct.max_nl];
        sint2I = &pct.newsintI_local[sidx2 + nidx * ct.num_states * ct.max_nl];


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

    for (idx = 0; idx <pct.P0_BASIS; idx++)
    {
        sumpsiR += (tmp_psi2R[idx] * tmp_psi1R[idx] + tmp_psi2I[idx] * tmp_psi1I[idx]);
        sumpsiI += (tmp_psi2I[idx] * tmp_psi1R[idx] - tmp_psi2R[idx] * tmp_psi1I[idx]);
    }



    *cR = ct.vel * sumpsiR + sumbetaR;
    *cI = ct.vel * sumpsiI + sumbetaI;

}

#endif

#else

/* Calculates coefficients that will be used in orthogonalizing wavefunctions*/

/*Loop over ions is parallelized in this function only certain parts of newsintR
 * are updated when wavefunctions are updated, so newsintR (projectors applied on wavefunctions)
 * should be recalculated (using betaxpsi) once new_psi is finished. This is done 
 * at the end of ortho_full.c*/

#if GAMMA_PT
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t * cR, rmg_double_t * cI)
{
    int ion, i, j, inh, sidx1, sidx2, nidx, oion;
    int idx, nh;
    rmg_double_t *tmp_psi2, *tmp_psi1;
    rmg_double_t *sint1R, *sint2R, *qqq;
    rmg_double_t sumpsi, sumbeta, sri;
    ION *iptr;
    SPECIES *sp;


    sidx1 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1 = sp1->psiR;
    tmp_psi2 = sp2->psiR;

    sumpsi = 0.0;
    sumbeta = 0.0;

    /*This loop should be the same for all ions - qqq and newsintR should have 
     * the same values everywhere*/
    /*This is parallelized over ions */
    nidx = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {
	oion = pct.owned_ions_list[ion];
        
        iptr = &ct.ions[oion];
        sp = &ct.sp[iptr->species];
       
        nh = sp->nh;
	
	/* Figure out index of owned ion in nonloc_ions_list array*/
	do {
	    
	    nidx++;
	    if (nidx >= pct.num_nonloc_ions)
		error_handler("Could not find matching entry in pct.nonloc_ions_list for owned ion %d", oion);
	
	} while (pct.nonloc_ions_list[nidx] != oion);
        
        qqq = pct.qqq[oion];

        /* get<beta|psi1> and <beta|psi2> */
        sint1R = &pct.newsintR_local[sidx1 + nidx * ct.num_states * ct.max_nl];
        sint2R = &pct.newsintR_local[sidx2 + nidx * ct.num_states * ct.max_nl];

        for (i = 0; i < nh; i++)
        {
            inh = i * nh;
            sri = sint1R[i];

            for (j = 0; j < nh; j++)
                sumbeta += qqq[inh + j] * sri * sint2R[j];
        }                       /*end for i */

    }                           /*end for ion */


    for (idx = 0; idx <pct.P0_BASIS; idx++)
        sumpsi += tmp_psi2[idx] * tmp_psi1[idx];



    *cR = -1.0 * (ct.vel * sumpsi + sumbeta);

}



#else
void ortho_get_coeff (STATE * sp1, STATE * sp2, int ist1, int ist2, int kidx, rmg_double_t * cR, rmg_double_t * cI)
{
    int ion, i, j, nidx, oion;
    int idx, nh;
    int sidx1, sidx2, inh;
    rmg_double_t *tmp_psi2R, *tmp_psi2I, *tmp_psi1R, *tmp_psi1I;
    rmg_double_t *sint1R, *sint1I, *sint2R, *sint2I, *qqq;
    rmg_double_t sumpsiR, sumpsiI, sumbetaR, sumbetaI, sri, sii;
    ION *iptr;
    SPECIES *sp;


    sidx1 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist1 * ct.max_nl;
    sidx2 = kidx * pct.num_nonloc_ions * ct.num_states * ct.max_nl + ist2 * ct.max_nl;


    tmp_psi1R = sp1->psiR;
    tmp_psi1I = sp1->psiI;
    tmp_psi2R = sp2->psiR;
    tmp_psi2I = sp2->psiI;

    sumpsiR = 0.0;
    sumpsiI = 0.0;
    sumbetaR = 0.0;
    sumbetaI = 0.0;

    /*This is parallelized over ions */
    nidx = -1;
    for (ion = 0; ion < pct.num_owned_ions; ion++)
    {
	oion = pct.owned_ions_list[ion];
        
        iptr = &ct.ions[oion];
        sp = &ct.sp[iptr->species];
       
        nh = sp->nh;
	
	/* Figure out index of owned ion in nonloc_ions_list array*/
	do {
	    
	    nidx++;
	    if (nidx >= pct.num_nonloc_ions)
		error_handler("Could not find matching entry in pct.nonloc_ions_list for owned ion %d", oion);
	
	} while (pct.nonloc_ions_list[nidx] != oion);
        
        qqq = pct.qqq[oion];

        /* get<beta|psi1> and <beta|psi2> */
        sint1R = &pct.newsintR_local[sidx1 + nidx * ct.num_states * ct.max_nl];
        sint1I = &pct.newsintI_local[sidx1 + nidx * ct.num_states * ct.max_nl];
        sint2R = &pct.newsintR_local[sidx2 + nidx * ct.num_states * ct.max_nl];
        sint2I = &pct.newsintI_local[sidx2 + nidx * ct.num_states * ct.max_nl];


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

    for (idx = 0; idx <pct.P0_BASIS; idx++)
    {
        sumpsiR += (tmp_psi2R[idx] * tmp_psi1R[idx] + tmp_psi2I[idx] * tmp_psi1I[idx]);
        sumpsiI += (tmp_psi2I[idx] * tmp_psi1R[idx] - tmp_psi2R[idx] * tmp_psi1I[idx]);
    }



    *cR = ct.vel * sumpsiR + sumbetaR;
    *cI = ct.vel * sumpsiI + sumbetaI;

}

#endif
#endif
