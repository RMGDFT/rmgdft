/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/***** QMD-MGDFT/get_dos.c *****
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
 *   void get_dos(STATE *states, rmg_double_t *rho)
 *   Generates pdos when CONVERGENCE = TRUE 
 * INPUTS
 *   states:  point to orbital structure (see main.h)
 *   rho:     old charge density
 * OUTPUT
 *   pdos
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
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"





void get_pdos (STATE * states, rmg_double_t Emin, rmg_double_t Emax, int E_POINTS)
{

    int istate, kpt, n, incx, idx, max_product, iene, P0_BASIS, FP0_BASIS;
    int PX0_GRID, PY0_GRID, PZ0_GRID;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int *ivec;
    int nh, icount, ncount, i, j, ion, gion, ix, iy, iz, ii, jj, kk, xoff, yoff, zoff;
    rmg_double_t *qnmI, *sintR, *qtpr;
    rmg_double_t t1, *work, *work_temp, *rho_temp, *rho_energy, *product, de, E;
    rmg_double_t time1;
    FILE *file;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();
    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();

#if !GAMMA_PT
    rmg_double_t *sintI;
#endif
    STATE *sp;
    ION *iptr;
    de = (Emax - Emin) / (E_POINTS - 1);

    P0_BASIS = get_P0_BASIS();
    FP0_BASIS = get_FP0_BASIS();

    my_calloc (work, P0_BASIS, rmg_double_t);

#if GAMMA_PT
    my_malloc (sintR, ct.max_nl, rmg_double_t);
#else
    my_malloc (sintR, 2 * ct.max_nl, rmg_double_t);
    sintI = sintR + ct.max_nl;
#endif
            
    max_product = (ct.max_nl + 1) * ct.max_nl / 2;
    my_malloc (product, max_product, rmg_double_t);

    my_malloc_init( rho_energy, E_POINTS * get_FNX_GRID(), rmg_double_t );
                                                                                              
    pe2xyz (pct.gridpe, &ii, &jj, &kk);
    xoff = ii * get_FPX0_GRID();

    /* scale charge accumulator */
    n = FP0_BASIS;
    incx = 1;


/***************Begin get pdos******************/

for (iene = 0; iene < E_POINTS; iene++)
{

	my_malloc_init (work_temp, P0_BASIS, rmg_double_t);
	my_malloc_init (rho_temp, FP0_BASIS, rmg_double_t);

	sp = ct.kp[0].kstate;
	E =  iene * de  + Emin;

	/* Loop over states and accumulate charge */
	for (istate = 0; istate < ct.num_states; istate++)
	{

		t1 = 5.6418959 * exp( -100.0 * (sp->eig[0] * Ha_eV - E) * (sp->eig[0] * Ha_eV - E) ) ;
                printf("t1 = %f, eig = %f , E = %f  \n", t1, sp->eig[0] * Ha_eV, E);

		for (idx = 0; idx < P0_BASIS; idx++)
		{
			work_temp[idx] += t1 * sp->psiR[idx] * sp->psiR[idx];
#if !GAMMA_PT
			work_temp[idx] += t1 * sp->psiI[idx] * sp->psiI[idx];
#endif
		}                   /* end for */

		sp++;
	}                       /*end for istate */



    /* Interpolate onto fine grid, result will be stored in rho*/
    time1 = my_crtc ();

    switch (ct.interp_flag)
    {
	    case 0:
		    pack_rho_ctof (work_temp, rho_temp);
		    break;
	    case 1:
		    bspline_interp_full (work_temp, rho_temp);
		    break;
	    case 2:
		    mg_prolong_MAX10 (rho_temp, work_temp, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_NX(), 6);
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

		    for (i=0; i < max_product; i++)
			    product[i] = 0.0;

		    for (kpt = 0; kpt < ct.num_kpts; kpt++)
		    {

			    sp = ct.kp[kpt].kstate;
			    /* Loop over states and accumulate charge */
			    for (istate = 0; istate < ct.num_states; istate++)
			    {
				    t1 = 5.6418959 * exp( -100.0 * (sp->eig[0] * Ha_eV - E) * (sp->eig[0] * Ha_eV - E) ) ;

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
					    rho_temp[ivec[icount]] += qtpr[icount] * product[idx];
				    }           /*end for icount */
				    idx++;
			    }               /*end for j */
		    }                   /*end for i */


	    }                       /*end if */

    }                           /*end for ion */


for (ix = 0; ix < get_FPX0_GRID(); ix++)
{
	for (iy = 0; iy < get_FPY0_GRID(); iy++)
	{
		for (iz = 0; iz < get_FPZ0_GRID(); iz++)
		{
			idx = iz + iy * get_FPZ0_GRID() + ix * get_FPZ0_GRID() * get_FPY0_GRID();
			rho_energy[iene * get_FNX_GRID() + ix + xoff] += rho_temp[idx];
		}
	}
}

    my_free (work_temp);
    my_free (rho_temp);

}                           /*end for iene*/
/***************End get pdos******************/



/***************begin plot pdos******************/
    iene = E_POINTS * get_FNX_GRID();
    global_sums (rho_energy, &iene, pct.grid_comm);
    if (pct.gridpe == 0)
    {
        double dx = ct.celldm[0] / get_NX_GRID();
        double x0 = 0.5 * ct.celldm[0];

        file = fopen ("dos.dat", "w");
        fprintf (file, "#     x[a0]      E[eV]          dos\n\n");
        for (iene = 0; iene < E_POINTS; iene++)
        {

            for (ix = 0; ix < get_NX_GRID(); ix++)
            {

                fprintf (file, " %10.6f %10.6f %12.6e\n",
                        ix * dx - x0, Emin+iene*de, rho_energy[iene * get_FNX_GRID() + ix * get_FG_NX()]);
            }
            fprintf (file, "\n");
        }

        fclose (file);
    }
/***************End plot pdos******************/


    /* release our memory */
    my_free (product);
    my_free (work);
    my_free (sintR);
    my_free(rho_energy); 

}                               /* end get_dos */


/******/
