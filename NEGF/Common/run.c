/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/run.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   void run()   
 *   Perform any initializations that are required and then
 *   enters the main driver loop. It also handles checkpointing and the
 *   output of intermediate results.
 * INPUTS
 *   nothing
 * OUTPUT
 *   Print standard output.
 * PARENTS
 *   md.c
 * CHILDREN
 *   
 * SEE ALSO
 *   md.h for structure defination
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "md.h"

extern REAL *vh_old, *vxc_old;

void run (STATE * states, STATE * states1)
{
    REAL time1, time2;
    FILE *file;
    int ix, iy, iz, idx;
    int size, iprobe, idx_delta, j;
    double *vbias;

    time1 = my_crtc ();
    /* initialize processor structure, decompose the processor into subdomains
     * pct.pe_kpoint * ( pct.pe_x, pct.pe_y, pct.pe_z) or
     * pct.pe_kpoint * ( pct.pe_column , pct.pe_row)
     */

    init_pe_on ();


    init_dimension ();

    pmo_init();
    if (ct.runflag == 100)
    {
        time1 = my_crtc ();
        lead_bandstructure ();
        time2 = my_crtc ();
        if (pct.gridpe == 0)
            printf ("\nTIME for lead bandstructure %f\n", time2 - time1);
        if (pct.gridpe == 0)
            printf ("\nband structrue file: band.dat\n");
    }
    else if (ct.runflag == 110)
    {
        time1 = my_crtc ();
        get_cond_frommatrix ();
        time2 = my_crtc ();
        if (pct.gridpe == 0)
            printf ("\nTIME for get_cond_dos %f\n", time2 - time1);
    }

    else 
    {
/*    if(ct.runflag == 100)
 *    {
 *        time1 = my_crtc();
 *        get_cond_lead();
 *        time2 = my_crtc();
 *        if(pct.gridpe ==0) printf("\n TIME for get_cond_dos %f", time2 - time1);
 *        exit(0);
 *    }
 */
        
        /* allocate memory for matrixs  */
        allocate_matrix_soft ();
        
        /* Perform some necessary initializations 
         * no matter localized or not  
         */
/*    my_malloc_init( vxc_old, FP0_BASIS, REAL );
 *   my_malloc_init( vh_old, FP0_BASIS, REAL );
 */
        vxc_old = vxc;
        vh_old = vh;
        
        init_soft (vh, rho, rhocore, rhoc, states, states1, vnuc, vxc, vh_old, vxc_old);
        if (ct.runflag == 200)
        {
            time1 = my_crtc ();
            get_dos (states);
            time2 = my_crtc ();
            if (pct.gridpe == 0)
                printf ("\n TIME for get_cond_dos %f", time2 - time1);
        }
        else 
        {
            
            if (pct.gridpe == 0)
                printf ("init_soft is done\n");

            /* total energy point = # of poles + ncircle + nmax_gq1 */

            size = (int) (cei.DELTA2/(2 * PI * cei.KT)) + cei.ncircle + cei.nmax_gq1 +10;

			for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
			{	
                my_malloc(lcr[iprobe].ene, size, complex double);
                my_malloc(lcr[iprobe].weight, size, complex double);
            }
			
            size = cei.nmax_gq2 + 10;
			for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
			{	
                for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                {	  
                    my_malloc(lcr[iprobe].lcr_ne[idx_delta - 1].ene_ne, size, complex double);
                    my_malloc(lcr[iprobe].lcr_ne[idx_delta - 1].weight_ne, size, complex double);
                }
            }
		   	
			
			for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
			{	
                set_energy_weight (lcr[iprobe].ene, lcr[iprobe].weight,
                    lcr[iprobe].bias + lcr[iprobe].EF_new, &lcr[iprobe].nenergy);
			}



			
			for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
			{
                j = 0;	
                for (idx_delta = 1; idx_delta <= cei.num_probe; idx_delta++)
                {	  
                    if(idx_delta != iprobe)
                    {

                        set_energy_weight_ne (lcr[iprobe].lcr_ne[j].ene_ne, lcr[iprobe].lcr_ne[j].weight_ne, 
                               lcr[idx_delta].bias + lcr[idx_delta].EF_new, 
                               lcr[iprobe].bias + lcr[iprobe].EF_new, &lcr[iprobe].lcr_ne[j].nenergy_ne);

                        j++;
                    }
                }
			} 
            
            
            if (ct.runflag == 1)
                get_all_kbpsi (states, states);
            
            /* Wait until everybody gets here */
            my_barrier ();

/*--------------------------------*/

            my_malloc_init( vbias, FPX0_GRID * FPY0_GRID, double );
            if (ct.runflag == 113) apply_potential_drop( vbias );

 
            printf (" apply_potential_drop is done :-) \n");


/*--------------------------------*/

            /* Dispatch to the correct driver routine */
            
            switch (ct.forceflag)
            {

            case MD_QUENCH:            /* Quench the electrons */
                quench (states, states1, vxc, vh, vnuc, vh_old, vxc_old, rho, rhoc, rhocore, vbias);
                break;

            default:
                error_handler ("Undefined MD method");
                
                
            }                           /* end switch */


            time1 = my_crtc () - time1;
            rmg_timings (TOTAL_TIME, time1);


            if(pct.psi1 != NULL ) my_free( pct.psi1 );
            if(pct.psi2 != NULL ) my_free( pct.psi2 );
            if(projectors != NULL ) my_free( projectors );


            /* Save data to output file */
            write_data (ct.outfile, vh, vxc, vh_old, vxc_old, rho, vbias, &states[0]);

            if (ct.runflag == 111)
                write_data_lead (ct.outfile, vh, vxc, vh_old, vxc_old, rho);


            my_barrier ();
            writeout_matrix_p ();

            my_barrier ();
            if (pct.gridpe == 0)
                printf ("\n Run done...\n");
            fflush (NULL);

            my_barrier();

            write_timings();

        }
    }

}                               /* end run */
