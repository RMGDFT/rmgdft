#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "main.h"
#include "common_prototypes.h"





void get_new_tau (STATE * states, double * tau)
{

    int istate, kpt, idx, P0_BASIS, PX0_GRID, PY0_GRID, PZ0_GRID;
    
    double t1, *work, *wfn; 
    double time1;

    STATE *sp;
    
    double *gx, *gy, *gz;
    double hxgrid, hygrid, hzgrid;


    P0_BASIS = get_P0_BASIS();
    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();

    hxgrid=get_hxgrid();
    hygrid=get_hygrid();
    hzgrid=get_hzgrid();
    
    my_malloc (wfn, P0_BASIS, double);
    my_malloc (gx, P0_BASIS, double);
    my_malloc (gy, P0_BASIS, double);
    my_malloc (gz, P0_BASIS, double);


    my_malloc (work, P0_BASIS, double);

            
   printf("\nPO_BASIS: %d, PX0_GRID: %d, PY0_GRID: %d, PZ0_GRID: %d", P0_BASIS, PX0_GRID, PY0_GRID, PZ0_GRID);
   printf("\nFPO_BASIS: %d, FPX0_GRID: %d, FPY0_GRID: %d, FPZ0_GRID: %d", get_FP0_BASIS(), get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID());

   for (idx=0; idx < P0_BASIS; idx++)
	work[idx] = 0.0;


    for (kpt = 0; kpt < ct.num_kpts; kpt++)
    {

        sp = ct.kp[kpt].kstate;

        /* Loop over states and accumulate charge */
        for (istate = 0; istate < ct.num_states; istate++)
        {

            t1 = sp->occupation[0] * ct.kp[kpt].kweight;
            
	    for (idx = 0; idx < P0_BASIS; idx++)
            {
                wfn[idx] = sp->psiR[idx];
		//printf("idx: %d, wfnR: %16.9e\n", idx, wfn[idx]);
            
	    }			/* end for */
    
	    /* Generate the gradient of the real part of the wavefunction*/
    	    app_grad (wfn, gx, gy, gz, PX0_GRID, PY0_GRID, PZ0_GRID, hxgrid, hygrid, hzgrid);

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                work[idx] += t1 * ( gx[idx] * gx[idx] + gy[idx] * gy[idx] + gz[idx] * gz[idx]);
		if (t1 * ( gx[idx] * gx[idx] + gy[idx] * gy[idx] + gz[idx] * gz[idx] ) < 0)
			printf("istate: %d, idx: %d, t1: %16.9f, gx: %16.9f, gy: %16.9f, gz: %16.9f", istate, idx, t1, gx[idx], gy[idx], gz[idx]);
            
	    }			/* end for */

#if !GAMMA_PT
	    for (idx = 0; idx < P0_BASIS; idx++)
            {
                wfn[idx] = sp->psiI[idx];
		//printf("idx: %d, wfnI: %16.9e\n", idx, wfn[idx]);
            
	    }			/* end for */
	    /* Generate the gradient of the imaginary part of the wavefunction*/
    	    app_grad (wfn, gx, gy, gz, PX0_GRID, PY0_GRID, PZ0_GRID, hxgrid, hygrid, hzgrid);

            for (idx = 0; idx < P0_BASIS; idx++)
            {
                work[idx] += t1 * ( gx[idx] * gx[idx] + gy[idx] * gy[idx] + gz[idx] * gz[idx]);
            
	    }                   /* end for */
#endif

            sp++;
        }                       /*end for istate */

    }                           /*end for kpt */


    for (idx = 0; idx < P0_BASIS; idx++)
	if (work[idx] < 0)
		printf("work: %16.9e\n", work[idx]);


    /* Interpolate onto fine grid, result will be stored in tau*/
    time1 = my_crtc ();

    switch (ct.interp_flag)
    {
        case 0:
            pack_rho_ctof (work, tau);
            break;
        case 1:
            bspline_interp_full (work, tau);
            break;
        case 2:
            mg_prolong_MAX10 (tau, work, get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_FG_RATIO(), 6);
            break;

        default:

            Dprintf ("tau interpolation is set to %d", ct.interp_flag);
            error_handler ("ct.interp_flag is set to %d. The valid values are 0, 1, 2", ct.interp_flag);


    } 

    //rmg_timings (INTERPOLATION_TIME, my_crtc () - time1);


    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
		//if ( (tau[idx] <0) &&  (abs(tau[idx]) < 1.e-10 ))
                //	tau[idx] = 0.0;
                tau[idx] = 0.5 * tau[idx];
		if (tau[idx] < 0)
			printf("after interpolation of work, idx: %d, tau: %16.9e\n", idx, tau[idx]);
            
    }                   /* end for */


    /* release our memory */
    my_free (wfn);
    my_free (work);
    my_free (gx);
    my_free (gy);
    my_free (gz);

}                               /* end get_new_tau */


/******/
