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
 *   main.h for structure defination
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "make_conf.h"
#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"



#include "my_scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"



extern double *vh_old, *vxc_old;

void Run (STATE * states, STATE * states1)
{
    int size, iprobe, idx_delta, i, j;
    double *vbias;

    /* initialize processor structure, decompose the processor into subdomains
     * pct.pe_kpoint * ( pct.pe_x, pct.pe_y, pct.pe_z) or
     * pct.pe_kpoint * ( pct.pe_column , pct.pe_row)
     */

    read_trans(&cei);
    read_LCR();
    init_pe_on ();


    init_dimension (&MXLLDA, &MXLCOL);

    pmo_init();
    if (ct.runflag == 100)
    {

        RmgTimer *RT1 = new RmgTimer("BANDSTRUCTURE Calculation");

        lead_bandstructure ();

        delete(RT1);

        if (pct.gridpe == 0)
            printf ("\nband structrue file: band.dat\n");
    }
    else if (ct.runflag == 110)
    {
#if GPU_ENABLED
        init_gpu();
#endif

        /* allocate memory for matrixs  */
        //get_cond_frommatrix ();
        RmgTimer *RT2 = new RmgTimer("Conductance Calculation");
        get_cond_frommatrix_kyz ();
        delete(RT2);

    }

    else 
    {

        /* allocate memory for matrixs  */
        allocate_matrix_soft ();

        /* Perform some necessary initializations 
         * no matter localized or not  
         */
        int FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);
        vxc_old = new double[FP0_BASIS];
        vh_old = new double[FP0_BASIS];


        RmgTimer *RT3 = new RmgTimer("1-TOTAL: init");
        InitNegf (vh, rho, rhocore, rhoc, rho_tf, states, states1, vnuc, vext, vxc, vh_old, vxc_old, states_distribute);
        delete(RT3);

        size = 1;
        for (i = 0; i < ct.num_blocks; i++) size = std::max(size, ct.block_dim[i] * ct.block_dim[i]);
        size = std::max(size, ct.num_states * (ct.state_end - ct.state_begin));
        size = std::max(size, pct.num_local_orbit * pct.num_local_orbit);

        work_matrix = new double[size];


        if (pct.gridpe == 0)
            printf ("init_soft is done\n");

        if (ct.runflag == 200)
        {
            get_dos(states);
        }
        else
            if (ct.runflag == 300)
            {
                get_cond_frommatrix_kyz ();

                /* it will automatically calculate and plot 3Ddos for each peak */
                printf ("\n peakNum = %d", peakNum);
                for (i = 0; i < peakNum; i++)
                {
                    printf ("\n peaks[%d] = %f", i, peaks[i]);
                    get_3Ddos (states, peaks[i]-0.0002 , peaks[i]+0.0002, 3, i);
                }
            }
            else 
            {


                /* total energy point = # of poles + ncircle + nmax_gq1 */

                size = (int) (cei.DELTA2/(2 * PI * cei.KT)) + cei.ncircle + cei.nmax_gq1 +10;

                for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
                {	
                    lcr[iprobe].ene = new DoubleC[size];
                    lcr[iprobe].weight = new DoubleC[size];
                }

                size = cei.nmax_gq2 + 10;
                for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
                {	
                    for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)
                    {	  
                        lcr[iprobe].lcr_ne[idx_delta - 1].ene_ne = new DoubleC[size];
                        lcr[iprobe].lcr_ne[idx_delta - 1].weight_ne = new DoubleC[size];
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


                /*--------------------------------*/

                int FPX0_GRID = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
                int FPY0_GRID = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);

                vbias = new double[FPX0_GRID * FPY0_GRID];
                for(i = 0; i < FPX0_GRID * FPY0_GRID; i++) vbias[i] = 0.0;
                if (ct.runflag == 113) apply_potential_drop( vbias );


                printf (" apply_potential_drop is done :-) \n");


                /*--------------------------------*/

                /* Dispatch to the correct driver routine */

                RmgTimer *RT4 = new RmgTimer("1-TOTAL: Quench");
                switch (ct.forceflag)
                {

                    case MD_QUENCH:            /* Quench the electrons */
                        QuenchNegf (states, states1, states_distribute, vxc, vh, vnuc, vext, vh_old, vxc_old, rho, rhoc, rhocore, rho_tf, vbias);

                        break;

                    default:

                        exit(0);

                }                           /* end switch */

                delete(RT4);




                /* Save data to output file */
                RmgTimer *RT5 = new RmgTimer("1-TOTAL: Write_data");
                write_data_negf (ct.outfile, vh, vxc, vh_old, vxc_old, rho, vbias, &states[0]);

                if (ct.runflag == 111)
                    write_data_lead (ct.outfile, vh, vxc, vh_old, vxc_old, rho);


                my_barrier ();
                writeout_matrix_p ();
                delete(RT5);

                my_barrier ();
                if (pct.gridpe == 0)
                    printf ("\n Run done...\n");
                fflush (NULL);

                my_barrier();


            }
    }

}                               /* end run */

