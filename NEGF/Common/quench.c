/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* QMD-MGDFT/quench.c *****
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
 *   void quench(STATE *states, REAL *vxc, REAL *vh, REAL *vnuc, 
 *               REAL *rho, REAL *rhocore, REAL *rhoc)
 *   For a fixed atomic configuration, quench the electrons to find 
 *   the minimum total energy 
 * INPUTS
 *   states: point to orbital structure (see main.h)
 *   vxc:    exchange correlation potential
 *   vh:     Hartree potential
 *   vnuc:   Pseudopotential 
 *   rho:    total valence charge density
 *   rhocore: core chare density only for non-linear core correction
 *   rhoc:   Gaussian compensating charge density
 * OUTPUT
 *   states, vxc, vh, rho are updated
 * PARENTS
 *   cdfastrlx.c fastrlx.c md.c
 * CHILDREN
 *   scf.c force.c get_te.c subdiag.c get_ke.c
 * SOURCE
 */



#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "pmo.h"


void quench (STATE * states, STATE * states1, STATE *states_distribute, REAL * vxc, REAL * vh, REAL * vnuc, REAL * vext,
             REAL * vh_old, REAL * vxc_old, REAL * rho, REAL * rhoc, REAL * rhocore, REAL * vbias)
{

    int outcount = 0;
    static int CONVERGENCE = FALSE;

    int st1, st2, st11, st22;
    int idx, idx1, nL, iprobe, jprobe;
    int ictxt, mb, nprow, npcol, myrow, mycol;
    int j, k, jj, kk, idx_delta, idx_C;
    int j1, k1, jdiff, kdiff;
    int idx2, FPYZ0_GRID;
    int i;

    double time1, time2;

    time1 = my_crtc ();


    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        vxc_old[idx] = vxc[idx];
        vh_old[idx] = vh[idx];
    }

    if (ct.runflag == 111 | ct.runflag == 112 |ct.runflag == 1121)   /* check */
    {
        for (iprobe = 2; iprobe <= cei.num_probe; iprobe++)
        {
            lcr[iprobe].nenergy = 0;
        }
    }



    idx1 = 0;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        idx = (lcr[iprobe].nenergy + pmo.npe_energy - 1) / pmo.npe_energy;
 
        for (idx_delta = 1; idx_delta < cei.num_probe; idx_delta++)               /* check */
        {
            idx += (lcr[iprobe].lcr_ne[idx_delta - 1].nenergy_ne + pmo.npe_energy - 1) / pmo.npe_energy;
        }
        for (jprobe = 1; jprobe <= cei.num_probe; jprobe++)                       /* check */
        {
        
        idx_C = cei.probe_in_block[jprobe-1];
            nL = pmo.mxllda_cond[idx_C] * pmo.mxlocc_cond[idx_C];
            idx1 += idx * nL;
        }	
    }


    my_malloc_init( sigma_all, idx1, complex double );

    if (ct.runflag != 111)
        sigma_all_energy_point (sigma_all);


    get_all_kbpsi (states, states);
    /* get lcr[0].S00 part */
    get_matB_soft (states, states1, work_matrix);

/*******************************************/



    whole_to_tri_p (lcr[0].Stri, work_matrix, ct.num_blocks, ct.block_dim);

/* ========= interaction between leads is zero ========== */
    zero_lead_image(lcr[0].Stri);





    ictxt = pmo.ictxt[pmo.myblacs];
    mb = pmo.mblock;


    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);


    if (ct.runflag == 111 && cei.num_probe == 2)
    {

        for(j=0; j < pmo.mxllda_lead[1]; j++)
        {
            for(k =0; k < pmo.mxlocc_lead[1]; k++)
            {

                jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
                kk = (k/mb ) * npcol * mb + mycol * mb + k - k/mb * mb;

                st11 = jj + lcr[1].num_states;
                st22 = kk + lcr[1].num_states;
                lcr[1].S00[j + k * pmo.mxllda_lead[1]] =
                    work_matrix[st11 * ct.num_states + st22];
                lcr[2].S00[j + k * pmo.mxllda_lead[1]] =
                    work_matrix[st11 * ct.num_states + st22];

                st11 = jj + lcr[1].num_states + lcr[0].num_states;
                lcr[1].S01[j + k * pmo.mxllda_lead[1]] =
                    work_matrix[st11 * ct.num_states + st22];
                lcr[1].SCL[j + k * pmo.mxllda_lead[1]] =
                    work_matrix[st11 * ct.num_states + st22];

                st11 = jj + lcr[1].num_states;
                st22 = kk + lcr[1].num_states + lcr[0].num_states;
                lcr[2].S01[j + k * pmo.mxllda_lead[1]] =
                    work_matrix[st11 * ct.num_states + st22];
                lcr[2].SCL[j + k * pmo.mxllda_lead[1]] =
                    work_matrix[st11 * ct.num_states + st22];

            }
        }

    }


    /* the corner elements of the matrix should be unchanged */
    setback_corner_matrix_S();  


    for (idx = 0; idx < FP0_BASIS; idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx] + vext[idx];


/*  apply_potential_drop( vtot ); */

    FPYZ0_GRID = FPY0_GRID * FPZ0_GRID;

    for (i = 0; i < FPX0_GRID; i++)
    {
        for (j = 0; j < FPY0_GRID; j++)
        {
            idx = j + i * FPY0_GRID;

            for (k = 0; k < FPZ0_GRID; k++)
            {
                idx2 = k + j * FPZ0_GRID + i * FPYZ0_GRID;
                vtot[idx2] += vbias[idx];
            }
        }
    }



    get_vtot_psi(vtot_c, vtot, FG_NX);

    get_ddd (vtot);


    /* get lcr[0].H00 part */
    get_Hij (states, states1, vtot_c, work_matrix);


    whole_to_tri_p (lcr[0].Htri, work_matrix, ct.num_blocks, ct.block_dim);


/* ========= interaction between L3-L4 is zero ========== */

      zero_lead_image(lcr[0].Htri);  

/* corner elements keep unchanged */
    setback_corner_matrix_H();  


 

    for (ct.steps = 0; ct.steps < ct.max_scf_steps; ct.steps++)
    {

        if (pct.gridpe == 0)
            printf ("\n\n\n ITERATION     %d\n", ct.steps);
        /* Perform a single self-consistent step */
        if (!CONVERGENCE)
        {

            scf (sigma_all, states, states_distribute, vxc, vh, vnuc, vext, rho, rhoc,
                    rhocore, vxc_old, vh_old, vbias, &CONVERGENCE);



        }

        if (!ct.steps)
            CONVERGENCE = FALSE;

        if (CONVERGENCE)
        {

            if (pct.gridpe == 0)
                printf ("\n\n     convergence has been achieved. stopping ...\n");


            break;

        }                       /* end if */

        /* Check if we need to output intermediate results */
        if (outcount >= ct.outcount)
        {

            outcount = 0;

        }                       /* end if */

        outcount++;

    }                           /* end for */


   /* added by shuchun for force calculation */
    if (ct.forceflag !=0 )
    {
        /* Calculate the force */
        force (rho, rho, rhoc, vh, vxc, vnuc, states);
        /* write out the force */
        if (pct.gridpe == 0)
            write_force ();
    }
    /* end of addition */




    if(sigma_all != NULL) my_free(sigma_all);

    if (pct.gridpe == 0)
        printf ("\n Quench is done \n");


    time2 = my_crtc ();
    rmg_timings (QUENCH_TIME, time2 - time1);

}                               /* end quench */

/******/


