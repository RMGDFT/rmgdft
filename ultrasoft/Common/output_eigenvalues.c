/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/****f* QMD-MGDFT/output.c *****
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
 *   void output_eigenvalues(STATE *states, int ikbs, int iscf)
 *   driver routine to write eigenvalues to stdout
 * INPUTS
 *   states: point to orbital structure (see main.h)
 *   ikbs: the index of the k-point in a bandstructure calculation
 *   iscf: the index of the scf iteration
 * OUTPUT
 *   no explicit output
 * PARENTS
 *   main.c quench.c bandstructure.c
 * CHILDREN
  * SOURCE
 */


#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"


void output_eigenvalues (STATE * states, int ikbs, int iscf)
{
    int ik, jk, nk, is, il;

    int bs = verify ("calculation_mode", "Band Structure Only");

    KPOINT *kpt;
    STATE *st;


    nk = (bs) ? 1 : ct.num_kpts;


    for (ik = 0; ik < nk; ik++)
    {

        if (bs)
        {
            kpt = &ct.kp[states[0].kidx];
            st = ct.kp[0].kstate;
            jk = ikbs;
        }
        else
        {
            kpt = &ct.kp[ik];
            st = ct.kp[ik].kstate;
            jk = ik;
        }


        printf ("\n\nKOHN SHAM EIGENVALUES [eV] AT K-POINT [%3d]:   %12.6f  %12.6f  %12.6f\n\n",
                jk, kpt->kpt[0], kpt->kpt[1], kpt->kpt[2]);

	if (pct.spin_flag)
        {
		/* Spin up*/
		printf("Spin up eigenvalues:\n");
        il = 0;
        for (is = 0; is < ct.num_states; is++)
        {
            if (is % 4 == 0)
                printf ("[kpt %3d %3d %3d]", jk, iscf, il++);

            printf ("   %9.4f [%5.3f]%s",
                    st[is].eig * Ha_eV, st[is].occupation, ((is % 4 == 3) ? "\n" : ""));
        }
        printf ("\n");

		/* Spin down*/
		printf("Spin down eigenvalues:\n");
		il = 0;
        	for (is = 0; is < ct.num_states_oppo; is++)
        	{
            		if (is % 4 == 0)
                		printf ("[kpt %3d %3d %3d]", jk, iscf, il++);

            		printf ("   %9.4f [%5.3f]%s",
                    		st[is].eig_oppo * Ha_eV, st[is].occupation_oppo, ((is % 4 == 3) ? "\n" : ""));
        	}
        	printf ("\n");

	}
	else
	{
		il = 0;
        	for (is = 0; is < ct.num_states; is++)
        	{
            		if (is % 4 == 0)
               			 printf ("[kpt %3d %3d %3d]", jk, iscf, il++);

            		printf ("   %9.4f [%5.3f]%s",
                    		st[is].eig * Ha_eV, st[is].occupation, ((is % 4 == 3) ? "\n" : ""));
        	}
        	printf ("\n");
	}
    }
}

/******/
