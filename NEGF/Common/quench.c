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
 *   void quench(STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc, 
 *               rmg_double_t *rho, rmg_double_t *rhocore, rmg_double_t *rhoc)
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
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"


void quench (STATE * states, STATE * states1, STATE *states_distribute, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc, rmg_double_t * vext,
             rmg_double_t * vh_old, rmg_double_t * vxc_old, rmg_double_t * rho, rmg_double_t * rhoc, rmg_double_t * rhocore, rmg_double_t * vbias)
{

    int outcount = 0;
    static int CONVERGENCE = FALSE;

    int st1, st2, st11, st22;
    int idx, idx1, nL, iprobe, jprobe;
    int j, k, jj, kk, idx_delta, idx_C;
    int j1, k1, jdiff, kdiff;
    int idx2, FPYZ0_GRID;
    int i;


    void *RT = BeginRmgTimer("2-Quench");


    for (idx = 0; idx < get_FP0_BASIS(); idx++)
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
    iprobe = cei.probe_noneq;
    for (iprobe = 1; iprobe <= cei.num_probe; iprobe++)
    {
        if(cei.probe_noneq > 0) iprobe = cei.probe_noneq;
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

        if(cei.probe_noneq > 0) break;
    }




    my_malloc_init( sigma_all, idx1, complex double );

    void *RT1 = BeginRmgTimer("2-Quench: sigma_all");
    if (ct.runflag != 111)
        sigma_all_energy_point (sigma_all, ct.kp[pct.kstart].kpt[1], ct.kp[pct.kstart].kpt[2]);
    my_barrier();
    EndRmgTimer(RT1);
    if(pct.gridpe==0) dprintf("\n sigma_all done");


    void *RT2 = BeginRmgTimer("2-Quench: kbpsi");
    get_all_kbpsi (states, states, ion_orbit_overlap_region_nl, projectors, kbpsi);
    EndRmgTimer(RT2);
    /* get lcr[0].S00 part */

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vh[idx] + vxc[idx] + vnuc[idx] + vext[idx];


    /*  apply_potential_drop( vtot ); */

    FPYZ0_GRID = get_FPY0_GRID() * get_FPZ0_GRID();

    for (i = 0; i < get_FPX0_GRID(); i++)
    {
        for (j = 0; j < get_FPY0_GRID(); j++)
        {
            idx = j + i * get_FPY0_GRID();

            for (k = 0; k < get_FPZ0_GRID(); k++)
            {
                idx2 = k + j * get_FPZ0_GRID() + i * FPYZ0_GRID;
                vtot[idx2] += vbias[idx];
            }
        }
    }



    get_vtot_psi(vtot_c, vtot, get_FG_RATIO());

    get_ddd (vtot);


    void *RT3 = BeginRmgTimer("2-Quench: set H and S first");
    /* get lcr[0].H00 part */
    get_HS(states, states1, vtot_c, Hij_00, Bij_00);


    for(i = 0; i < pmo.ntot;i++) lcr[0].Htri[i] = 0.0;
    for(i = 0; i < pmo.ntot;i++) lcr[0].Stri[i] = 0.0;
    row_to_tri_p (lcr[0].Htri, Hij_00, ct.num_blocks, ct.block_dim);


    /* ========= interaction between L3-L4 is zero ========== */

    zero_lead_image(lcr[0].Htri);  

    /* corner elements keep unchanged */
    setback_corner_matrix_H();  

    /*******************************************/


    row_to_tri_p (lcr[0].Stri, Bij_00, ct.num_blocks, ct.block_dim);

    /* ========= interaction between leads is zero ========== */
    zero_lead_image(lcr[0].Stri);


    if (ct.runflag == 111 && cei.num_probe == 2)
    {
        int size_of_matrix = pmo.mxllda_lead[1] * pmo.mxlocc_lead[1];

        int numst = lcr[1].num_states, ione=1, *desca;
        double one = 1.0, zero = 0.0;
        desca = &pmo.desc_lead[0];

        dcopy (&size_of_matrix, &lcr[0].Stri[2*size_of_matrix], &ione, lcr[1].S00, &ione);
        dcopy (&size_of_matrix, &lcr[0].Stri[2*size_of_matrix], &ione, lcr[2].S00, &ione);


        PDTRAN(&numst, &numst, &one, &lcr[0].Stri[size_of_matrix], &ione, &ione, desca,
                                &zero, lcr[1].S01, &ione, &ione, desca);
        dcopy (&size_of_matrix, lcr[1].S01, &ione, lcr[1].SCL, &ione);


        dcopy (&size_of_matrix, &lcr[0].Stri[3*size_of_matrix], &ione, lcr[2].S01, &ione);
        dcopy (&size_of_matrix, lcr[2].S01, &ione, lcr[2].SCL, &ione);

    }


    /* the corner elements of the matrix should be unchanged */
    setback_corner_matrix_S();  



    EndRmgTimer(RT3);


    for (ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {

        if (pct.gridpe == 0)
            printf ("\n\n\n ITERATION     %d\n", ct.scf_steps);
        /* Perform a single self-consistent step */
        if (!CONVERGENCE)
        {

            void *RT4 = BeginRmgTimer("2-Quench: SCF");
            scf (sigma_all, states, states_distribute, vxc, vh, vnuc, vext, rho, rhoc,
                    rhocore, vxc_old, vh_old, vbias, &CONVERGENCE);

            EndRmgTimer(RT4);


        }

        if (!ct.scf_steps)
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

    EndRmgTimer(RT);


}                               /* end quench */

/******/


