/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/***** QMD-MGDFT/get_rho.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 *   Version: 2.1.5
 * COPYRIGHT
 *   Copyright (C) 1995  Emil Briggs
 *   Copyright (C) 1998  Emil Briggs, Charles Brabec, Mark Wensell, 
 *                       Dan Sullivan, Chris Rapcewicz, Jerzy Bernholc
 *   Copyright (C) 2001  Emil Briggs, Wenchang Lu,
 *                       Marco Buongiorno Nardelli,Charles Brabec, 
 *                       Mark Wensell,Dan Sullivan, Chris Rapcewicz,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void get_new_rho(STATE *states, REAL *rho)
 *   Generates new charge density
 * INPUTS
 *   states:  point to orbital structure (see main.h)
 *   rho:     old charge density
 * OUTPUT
 *   rho:    updated charge density
 * PARENTS
 *   scf.c
 * CHILDREN
 *   gather_psi.c symmetrize_rho.f
 * SOURCE
 */


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"





void get_new_rho (STATE * states, REAL * rho)
{

    int istate, kpt, n, incx, idx;
    int *ivec;
    int nh, icount, ncount, i, j, ion, gion;
    REAL *qnmI, *sintR, *qtpr;
    REAL t1, *work, *product;
    REAL time1;
#if !GAMMA_PT
    REAL *sintI;
#endif
    STATE *sp;
    ION *iptr;

    my_calloc (work, P0_BASIS, REAL);

#if GAMMA_PT
    my_malloc (sintR, ct.max_nl, REAL);
#else
    my_malloc (sintR, 2 * ct.max_nl, REAL);
    sintI = sintR + ct.max_nl;
#endif

    /* scale charge accumulator */
    n = FP0_BASIS;
    incx = 1;


    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        sp = ct.kp[kpt].kstate;

        /* Loop over states and accumulate charge */
        for (istate = 0; istate < ct.num_states; istate++)
        {

            t1 = sp->occupation[0] * ct.kp[kpt].kweight;

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                work[idx] += t1 * sp->psiR[idx] * sp->psiR[idx];
#if !GAMMA_PT
                work[idx] += t1 * sp->psiI[idx] * sp->psiI[idx];
#endif
            }                   /* end for */

            sp++;
        }                       /*end for istate */

    }                           /*end for kpt */


    /* Interpolate onto fine grid, result will be stored in rho*/
    time1 = my_crtc ();

    switch (ct.interp_flag)
    {
        case 0:
            pack_rho_ctof (work, rho);
            break;
        case 1:
            bspline_interp_full (work, rho);
            break;
        case 2:
            mg_prolong_MAX10 (rho, work, FPX0_GRID, FPY0_GRID, FPZ0_GRID, PX0_GRID, PY0_GRID, PZ0_GRID, FG_NX, 6);
            break;

        default:

            Dprintf ("charge interpolation is set to %d", ct.interp_flag);
            error_handler ("ct.interp_flag is set to %d. The valid values are 0, 1, 2", ct.interp_flag);


    } 

    rmg_timings (INTERPOLATION_TIME, my_crtc () - time1);



    for (ion = 0; ion < pct.num_nonloc_ions; ion++)
    {
        gion = pct.nonloc_ions_list[ion];
        
        if (pct.Qidxptrlen[gion])
        {
            
            iptr = &ct.ions[gion];
       
            nh = ct.sp[iptr->species].nh;
            
            ivec = pct.Qindex[gion];
            ncount = pct.Qidxptrlen[gion];
            qnmI = pct.augfunc[gion];

            my_calloc (product, (nh + 1) * nh / 2, REAL);

            for (kpt = 0; kpt < ct.num_kpts; kpt++)
            {

                sp = ct.kp[kpt].kstate;
                /* Loop over states and accumulate charge */
                for (istate = 0; istate < ct.num_states; istate++)
                {
                    t1 = sp->occupation[0] * ct.kp[kpt].kweight;

                    for (i = 0; i < ct.max_nl; i++)
                    {
                        sintR[i] =
                            pct.newsintR_local[kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl 
                            + ion * ct.num_states * ct.max_nl + istate * ct.max_nl + i];
#if !GAMMA_PT
                        sintI[i] =
                            pct.newsintI_local[kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl 
                            + ion * ct.num_states * ct.max_nl + istate * ct.max_nl + i];
#endif
                    }               /*end for i */

                    idx = 0;
                    for (i = 0; i < nh; i++)
                    {
                        for (j = i; j < nh; j++)
                        {
#if GAMMA_PT
                            if (i == j)
                                product[idx] += t1 * sintR[i] * sintR[j];
                            else
                                product[idx] += 2 * t1 * sintR[i] * sintR[j];
#else
                            if (i == j)
                                product[idx] += t1 * (sintR[i] * sintR[j] + sintI[i] * sintI[j]);
                            else
                                product[idx] += 2 * t1 * (sintR[i] * sintR[j] + sintI[i] * sintI[j]);
#endif
                            idx++;
                        }           /*end for j */
                    }               /*end for i */
                    sp++;
                }                   /*end for istate */
            }                       /*end for kpt */


            idx = 0;
            for (i = 0; i < nh; i++)
            {
                for (j = i; j < nh; j++)
                {
                    qtpr = qnmI + idx * ncount;
                    for (icount = 0; icount < ncount; icount++)
                    {
                        rho[ivec[icount]] += qtpr[icount] * product[idx];
                    }           /*end for icount */
                    idx++;
                }               /*end for j */
            }                   /*end for i */

            my_free (product);

        }                       /*end if */

    }                           /*end for ion */



#if !GAMMA_PT
    /* Symmetrize the density */
    symmetrize_rho ((FP0_GRID *) rho);
#endif

    /* Check total charge. */
    ct.tcharge = ZERO;
    for (idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /* ct.tcharge = real_sum_all (ct.tcharge); */
    ct.tcharge = real_sum_all (ct.tcharge, pct.img_comm);  
    ct.tcharge = ct.tcharge * ct.vel_f;

    /* Renormalize charge, there could be some discrpancy because of interpolation */
    t1 = ct.nel / ct.tcharge;
    if (pct.imgpe == 0)
        printf ("\n get_new_rho: Normalization constant for new charge is %f", t1);
    QMD_sscal (n, t1, rho, incx);

    /*Update ct.tcharge, do not really recalculate it, just mutltiply it by normalization constant */
    ct.tcharge *= t1;

    /* release our memory */
    my_free (work);
    my_free (sintR);

}                               /* end get_rho */


/******/
