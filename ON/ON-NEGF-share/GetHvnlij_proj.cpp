/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*
get_Hvnlij:

Get the elements of the Hamiltonian matrix due to the non-local
potential, and add them into Aij.


 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"
#include "blas.h"
#include "Kbpsi.h"



void GetHvnlij_proj(double *Aij, double *Bij, double *Kbpsi_mat1, double *Kbpsi_mat2, 
        int num_orb1, int num_orb2, int num_proj, bool flag_overlap)
{
    double zero = 0.0, one = 1.0;
    double *dnmI, *qnmI;
    double *temA;

    if(num_orb1 < 1 || num_orb2 < 1 || num_proj < 1) return;
    temA = new double[num_orb1 * ct.max_nl];

    int proj_count = 0;
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        ION *iptr = &Atoms[ion];
        SPECIES *sp = &Species[iptr->species];

        int nh = sp->num_projectors;

        if(nh == 0) continue;
        dnmI = pct.dnmI[ion];
        qnmI = pct.qqq[ion];




        dgemm ("T", "N", &num_orb1, &nh, &nh, &one, &Kbpsi_mat1[proj_count], &num_proj, dnmI, &nh, &zero, temA, &num_orb1);
        dgemm ("N", "N", &num_orb1, &num_orb2, &nh, &one, temA, &num_orb1, &Kbpsi_mat2[proj_count], &num_proj, &one, Aij, &num_orb1);


        if(!ct.norm_conserving_pp && flag_overlap)
        {
            dgemm ("T", "N", &num_orb1, &nh, &nh, &one, &Kbpsi_mat1[proj_count], &num_proj, qnmI, &nh, &zero, temA, &num_orb1);
            dgemm ("N", "N", &num_orb1, &num_orb2, &nh, &one, temA, &num_orb1, &Kbpsi_mat2[proj_count], &num_proj, &one, Aij, &num_orb1);
        }

        proj_count += nh;
    }


    delete [] temA;

}
