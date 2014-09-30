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



void output_band_plot()
{
    int is;
    FILE *bs_f;
    my_fopen (bs_f, "band.dat", "w");

    int ik;

    for (is = 0; is < ct.num_states; is++)
    {
        for(ik = 0; ik < ct.num_kpts; ik++)
        {

            fprintf (bs_f, "\n %4d  %16.8f ", ik, ct.kp[ik].kstate[is].eig[0] * Ha_eV);


        }
        
        fprintf (bs_f, "\n &&");
    }

    fclose (bs_f);
}
