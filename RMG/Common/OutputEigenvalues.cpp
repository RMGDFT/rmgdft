/************************** SVN Revision Information **************************
 **    $Id: output_eigenvalues.c 2402 2014-07-24 00:37:44Z ebriggs $    **
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


#include <complex>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "Kpoint.h"
#include "transition.h"

template void OutputEigenvalues<double> (Kpoint<double> **, int, int);
template void OutputEigenvalues<std::complex<double> >(Kpoint<std::complex<double> > **, int, int);

template <typename KpointType>
void OutputEigenvalues (Kpoint<KpointType> **Kptr, int ikbs, int iscf)
{
    int ik, jk, nk, is, il, idx, nspin = (ct.spin_flag + 1);

    int bs = verify ("calculation_mode", "Band Structure Only");

    Kpoint<KpointType> *kptr;
    nk = (bs) ? 1 : ct.num_kpts;

    for (ik = 0; ik < nk; ik++)
    {

        if (bs)
        {
            kptr = Kptr[0];
            jk = ikbs;
        }
        else
        {
            kptr = Kptr[ik];
            jk = ik;
        }


        rmg_printf ("\n\nKOHN SHAM EIGENVALUES [eV] AT K-POINT [%3d]:   %12.6f  %12.6f  %12.6f\n\n",
                jk, kptr->kpt[0], kptr->kpt[1], kptr->kpt[2]);

        for (idx = 0; idx < nspin; idx++)
        {
            if ( (nspin == 2) && (idx == 0))	
                rmg_printf("\n------------- SPIN UP ---------------\n\n");
            else if ( (nspin == 2) && (idx == 1))	
                rmg_printf("\n------------ SPIN DOWN --------------\n\n"); 
            il = 0;
            for (is = 0; is < ct.num_states; is++)
            {
                if (is % 4 == 0)
                    rmg_printf ("[kpt %3d %3d %3d]", jk, iscf, il++);

                rmg_printf ("   %8.4f [%5.3f]%s",
                        kptr->Kstates[is].eig[idx] * Ha_eV, kptr->Kstates[is].occupation[idx], ((is % 4 == 3) ? "\n" : ""));
            }
            rmg_printf ("\n");
        }
    }
}

