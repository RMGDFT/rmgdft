/************************** SVN Revision Information **************************
 **    $Id: get_rho.c 1251 2011-02-07 14:40:15Z yanli $    **
******************************************************************************/

/****f* QMD-MGDFT/get_rho.c *****
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
 *   void get_rho(STATE *states, REAL *rho)
 *   Generates new charge density and mix linearly with old one.
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


/*

                         get_rho.c




*/




#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"




#if MPI

void get_rho (STATE * states, REAL * rho, REAL * rhocore)
{

    int istate, kpt, n, incx, idx;
    int *ivec, size;
    int nh, icount, ncount, i, j, ion;
    REAL *qnmI, *sintR, *qtpr;
    REAL t1, *work, *tmp_psiR, *sum, *tmp, *product;
    REAL time1, min, min2;
#if !GAMMA_PT
    REAL *sintI, *tmp_psiI;
#endif
    STATE *sp;
    ION *iptr;

    my_malloc (work, P0_BASIS, REAL);
    for (idx = 0; idx < P0_BASIS; idx++)
        work[idx] = 0.0;

#if GAMMA_PT
    my_malloc (tmp_psiR, P0_BASIS, REAL);
    my_malloc (sintR, ct.max_nl, REAL);
#else
    my_malloc (sintR, 2 * ct.max_nl, REAL);
    sintI = sintR + ct.max_nl;
    my_malloc (tmp_psiR, 2 * P0_BASIS, REAL);
    tmp_psiI = tmp_psiR + P0_BASIS;
#endif

    size = (FPX0_GRID + 2) * (FPY0_GRID + 2) * (FPZ0_GRID + 2);
    my_malloc (sum, 2 * size, REAL);
    tmp = sum + size;

    /* scale charge accumulator */
    n = FP0_BASIS;
    incx = 1;

    t1 = ONE - ct.mix;
    QMD_sscal (n, t1, rho, incx);


    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        sp = ct.kp[kpt].kstate;

        /* Loop over states and accumulate charge */
        for (istate = 0; istate < ct.num_states; istate++)
        {

            t1 = ct.mix * sp->occupation * ct.kp[kpt].kweight;
            gather_psi (tmp_psiR, NULL, sp, 0);

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                work[idx] += t1 * tmp_psiR[idx] * tmp_psiR[idx];
            }                   /* end for */
#if !GAMMA_PT
            gather_psi (NULL, tmp_psiI, sp, 0);
            for (idx = 0; idx < P0_BASIS; idx++)
            {
                work[idx] += t1 * tmp_psiI[idx] * tmp_psiI[idx];
            }                   /* end for */
#endif
            sp++;
        }                       /*end for istate */

    }                           /*end for kpt */

    for (i = 0; i < size; i++)
    {
        sum[i] = 0.0;
        tmp[i] = 0.0;
    }

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        ivec = pct.Qindex[ion];
        nh = pct.prj_per_ion[ion];
        ncount = pct.Qidxptrlen[ion];
        qnmI = pct.augfunc[ion];
        iptr = &ct.ions[ion];

        my_malloc (product, (nh + 1) * nh / 2, REAL);
        for (idx = 0; idx < (nh + 1) * nh / 2; idx++)
            product[idx] = 0.0;

        for (kpt = 0; kpt < ct.num_kpts; kpt++)
        {

            sp = ct.kp[kpt].kstate;
            /* Loop over states and accumulate charge */
            for (istate = 0; istate < ct.num_states; istate++)
            {
                t1 = ct.mix * sp->occupation * ct.kp[kpt].kweight;

                for (i = 0; i < ct.max_nl; i++)
                {
                    sintR[i] =
                        iptr->newsintR[kpt * ct.num_ions * ct.num_states * ct.max_nl +
                                       istate * ct.max_nl + i];
#if !GAMMA_PT
                    sintI[i] =
                        iptr->newsintI[kpt * ct.num_ions * ct.num_states * ct.max_nl +
                                       istate * ct.max_nl + i];
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

        if (pct.Qidxptrlen[ion])
        {

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

        }                       /*end if */
        my_free (product);

    }                           /*end for ion */

    time1 = my_crtc ();

    ct.tcharge = ZERO;
    for (idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];

    /*Do required type of interpolation */
    if (!ct.interp_flag)
        pack_rho_ctof (work, sum);
    else
        bspline_interp_full (work, sum);

    rmg_timings (INTERPOLATION_TIME, my_crtc () - time1, 0);
/*    pack_rho_ctof(work,tmp);
    pack_ptos(sum,tmp,FPX0_GRID,FPY0_GRID,FPZ0_GRID);
    trade_images(sum,FPX0_GRID,FPY0_GRID,FPZ0_GRID,pct.neighbors);
    app_smooth1((FS0_GRID *)sum,(FS0_GRID *)tmp);
    pack_stop(tmp,sum,FPX0_GRID,FPY0_GRID,FPZ0_GRID);

    print_density_z_direction(sum,FPX0_GRID,FPY0_GRID,FPZ0_GRID,20);
    print_density_z_direction(rho,FPX0_GRID,FPY0_GRID,FPZ0_GRID,20);
*/

    QMD_saxpy (n, ONE, sum, incx, rho, incx);

#if !GAMMA_PT
    /* Symmetrize the density */
    symmetrize_rho ((FP0_GRID *) rho);
#endif

    /*Find charge minimum */
    min = ZERO;
    min2 = ZERO;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        if (rho[idx] < min)
            min = rho[idx];
        
	/*Here we check charge density with rhocore added*/
	if (pct.spin_flag)
	{
		if ((rho[idx] + rhocore[idx]/2.0) < min2)
            		min2 = rho[idx] + rhocore[idx]/2.0;
	}
	else
	{
		if ((rho[idx] + rhocore[idx]) < min2)
            		min2 = rho[idx] + rhocore[idx];
	}
    }

    /*Find absolute minimum from all PEs in MPI_COMM_WORLD*/
    if (pct.spin_flag)
    {
    	min = real_min_all_spin (min);
    	min2 = real_min_all_spin (min2);
    }
    else
    {
    	min = real_min_all (min);
    	min2 = real_min_all (min2);
    }

    if ((pct.imgpe == 0) && (min < ZERO)){
         printf ("\n\n Charge density is NEGATIVE after interpolation, minimum is %e", min); 
         printf ("\n Minimum charge density with core charge added is %e", min2); 
    }


    /* Check total charge. */
    ct.tcharge = ZERO;
    for (idx = 0; idx < FP0_BASIS; idx++)
        ct.tcharge += rho[idx];


    if (pct.spin_flag)
    	ct.tcharge = real_sum_all_spin (ct.tcharge);  /* for spin case sum over MPI_COMM_WORLD*/
    else
    	ct.tcharge = real_sum_all (ct.tcharge);  /* for pairwise case sum over pct.grid_comm*/
     
    
    ct.tcharge = ct.tcharge * ct.vel_f;


    /* Renormalize charge */
    t1 = ct.nel / ct.tcharge;    /* should be close to 1.0*/
    if (pct.imgpe == 0)
        printf ("\n get_rho: Normalization constant is %f", t1);
    QMD_sscal (n, t1, rho, incx);

    /*Update ct.tcharge, do not really recalculate it, just mutltiply it by normalization constant */
    ct.tcharge *= t1;

    /* release our memory */
    my_free (tmp_psiR);
    my_free (work);
    my_free (sintR);
    my_free (sum);

}                               /* end get_rho */

#else


void get_rho (STATE * states, REAL * rho)
{
    error_handler ("This functions works only in MPI mode");
}
#endif

/******/
