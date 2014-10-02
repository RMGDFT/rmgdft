/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/bandstructure.c *****
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
 *   void bandstructure(STATE *states, rmg_double_t *vxc, rmg_double_t *vh, rmg_double_t *vnuc)
 *
 *   After we got converged rho, vxc ... by specicial kpoints, 
 *   calculate the band structure for some specified k-point
 * INPUTS
 *   states: points to orbital structure
 *   vxc: exchange correlation potential
 *   vh:  Hartree potential
 *   vnuc: pseudopotential
 * OUTPUT
 *   states are updated
 * PARENTS
 *   main.c
 * CHILDREN
 *   mg_eig_state.c pe0_write_eigenvalues.c
 * SOURCE
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


static void write_bs_eigenvalues (char *fn, STATE * states, int ik);



void bandstructure (STATE * states, rmg_double_t * vxc, rmg_double_t * vh, rmg_double_t * vnuc)
{

// Fix up later
#if 0
    int ik, st1, idx;
    rmg_double_t *vtot, *vtot_psi;
    rmg_double_t time1, time2;
    char newname[MAX_PATH + 20];
    char name[] = "wavefunc.dat";

    /* Make the new output file name */
    sprintf (newname, "%s%d", name, pct.gridpe);

    /* Open output files for wfs and energies */
    /*my_open( wave_f, newname, O_CREAT | O_TRUNC | O_RDWR, S_IREAD | S_IWRITE ); */


    my_malloc (vtot, get_FP0_BASIS(), rmg_double_t);
    my_malloc (vtot_psi, get_P0_BASIS(), rmg_double_t);

    /*  get total potential  */
    for (idx = 0; idx < get_FP0_BASIS(); idx++)
        vtot[idx] = vxc[idx] + vh[idx] + vnuc[idx];

    get_ddd (vtot);
    get_vtot_psi (vtot_psi, vtot, get_FG_RATIO());

    /*  loop for k point */
    for (ik = 0; ik < ct.num_kpts; ik++)
    {

        betaxpsi1 (states, ik);

        /* Reinitialize some stuff in states object */
        for (st1 = 0; st1 < ct.num_states; st1++)
        {
            states[st1].kidx = ik;
            states[st1].istate = st1;
            states[st1].vel = get_vel();
            states[st1].eig[0] = 0.00;
            states[st1].res = 0.00;
        }


        /* Do iterations */
        for (ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
        {
            time1 = my_crtc ();

            /* Update the wavefunctions */
	    /* Maybe important: mg_eig_state calls app_nl_state, which by default uses
	     * mixed <beta|psi> while this function requires non-mixed dot products. This works, since
	     * bandstructure calculation is run as a standalone process with a restart. The behavior might break
	     * if this function is directly called at the end without restart. A possible fix might be to set up
	     * projector mixing to 1.0 */
            for (st1 = 0; st1 < ct.num_states; st1++)
                mg_eig_state (&states[st1], 0, vtot_psi);

            time2 = my_crtc ();

            betaxpsi1 (states, ik);

            ortho (states, ik);
            
            betaxpsi1 (states, ik);

            if (ct.sortflag)
                sortpsi (ct.kp[0].kstate);

            /* Output the energies */
            output_eigenvalues (states, ik, ct.scf_steps);
            printf ("\nTotal charge in supercell = %16.8f\n", ct.tcharge);

            wvfn_residual (states);

        }


        /* Dump converged wavefunctions */
        /*output_wave( states, ik, wave_f ); */

        /* Dump eigenvalues for this k-point */
        if (pct.gridpe == 0)
            write_bs_eigenvalues ("eigenvalues.dat", states, ik);


    }                           /* end for ik */


    /* Close output files */
    /*close (wave_f); */

    /* Force change mode of output file */
    /*
       chmod (newname, S_IREAD | S_IWRITE );
     */

    /* release memory */
    my_free (vtot);
    my_free (vtot_psi);
#endif

}                               /* end bandstructure  */

/********************************************/





static void write_bs_eigenvalues (char *fn, STATE * states, int ik)
{
    int is;
    FILE *bs_f;
    if (ik == 0)
        my_fopen (bs_f, fn, "w");
    else
        my_fopen (bs_f, fn, "a");


    fprintf (bs_f, "%4d  %8.4f  %8.4f  %8.4f",
             ik, ct.kp[ik].kpt[0], ct.kp[ik].kpt[1], ct.kp[ik].kpt[2]);

    for (is = 0; is < ct.num_states; is++)
        fprintf (bs_f, "  %8.4f", states[is].eig[0] * Ha_eV);

    fprintf (bs_f, "\n");

    fclose (bs_f);
}
