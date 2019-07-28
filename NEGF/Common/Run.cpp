#include "negf_prototypes.h"
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

#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "LCR.h"
#include "prototypes_on.h"
#include "prototypes_negf.h"
#include "init_var.h"



#include "Scalapack.h"
#include "blas.h"
#include "Kbpsi.h"
#include "FiniteDiff.h"



extern double *vh_old, *vxc_old;

void Run (STATE * states, STATE * states1, std::unordered_map<std::string, InputKey *>& ControlMap)
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

#if GPU_ENABLED
    init_gpu();
#endif
    if (ct.runflag == 100)
    {

        RmgTimer *RT1 = new RmgTimer("BANDSTRUCTURE Calculation");

        lead_bandstructure ();

        delete(RT1);

        if (pct.gridpe == 0)
            rmg_printf ("\nband structrue file: band.dat\n");
    }
    else if (ct.runflag == 110)
    {

        /* allocate memory for matrixs  */
        //get_cond_frommatrix ();
        allocate_matrix_soft ();
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
        InitNegf (vh, rho, rhocore, rhoc, rho_tf, states, states1, vnuc, vext, vxc, vh_old, vxc_old, ControlMap);
        delete(RT3);

        size = 1;
        for (i = 0; i < ct.num_blocks; i++) size = std::max(size, ct.block_dim[i] * ct.block_dim[i]);
        size = std::max(size, ct.num_states * (ct.state_end - ct.state_begin));

        work_matrix = new double[size];



        if (pct.imgpe == 0)
            rmg_printf ("init_soft is done\n");

        if (ct.runflag == 200)
        {
            RmgTimer *RT0 = new RmgTimer("2-SCF: orbital_comm");
            OrbitalComm(states);
            delete(RT0);

            RmgTimer *RTk = new RmgTimer("2-SCF: kbpsi");
            KbpsiUpdate(states);
            delete(RTk);

            get_dos(states);
            return;
        }
        if (ct.runflag == 300)
        {
            RmgTimer *RT0 = new RmgTimer("2-SCF: orbital_comm");
            OrbitalComm(states);
            delete(RT0);

            RmgTimer *RTk = new RmgTimer("2-SCF: kbpsi");
            KbpsiUpdate(states);
            delete(RTk);


            get_cond_frommatrix_kyz ();
            /* it will automatically calculate and plot 3Ddos for each peak */
            for (i = 0; i < peakNum; i++)
            {
                get_3Ddos (states, peaks[i]-0.0002 , peaks[i]+0.0002, 3, i);
            }
            return;
        }


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

                    if(ct.runflag == 113)
                    {
                        set_energy_weight_ne (lcr[iprobe].lcr_ne[j].ene_ne, lcr[iprobe].lcr_ne[j].weight_ne, 
                                lcr[idx_delta].bias + lcr[idx_delta].EF_new, 
                                lcr[iprobe].bias + lcr[iprobe].EF_new, &lcr[iprobe].lcr_ne[j].nenergy_ne);
                    }
                    else
                    {
                        lcr[iprobe].lcr_ne[j].nenergy_ne = 0;
                    }


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


        if(pct.imgpe == 0) printf (" apply_potential_drop is done :-) \n");


        /*--------------------------------*/

        /* Dispatch to the correct driver routine */

        RmgTimer *RT4 = new RmgTimer("1-TOTAL: Quench");
        switch (ct.forceflag)
        {

            case MD_QUENCH:            /* Quench the electrons */
                if (pct.imgpe == 0)
                    printf ("\n quench start...\n");
                QuenchNegf (states, states1, vxc, vh, vnuc, vext, vh_old, vxc_old, rho, rhoc, rhocore, rho_tf, vbias);
                if (pct.imgpe == 0)
                    printf ("\n quench done...\n");

                break;

            default:

                exit(0);

        }                           /* end switch */

        delete(RT4);




        /* Save data to output file */
        RmgTimer *RT5 = new RmgTimer("1-TOTAL: Write_data");

        if (ct.runflag == 111)
            write_data_lead (ct.outfile, vh, vxc, vh_old, vxc_old, rho);


        MPI_Barrier(pct.img_comm);
        writeout_matrix_p ();
        delete(RT5);

        MPI_Barrier(pct.img_comm);
        if (pct.imgpe == 0)
            printf ("\n Run done...\n");
        fflush (NULL);

        MPI_Barrier(pct.img_comm);


    }
}
