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



void GetHvnlij_dis(double *Aij, double *Bij, double *Kbpsi_mat, int num_orb, int num_proj)
{
    int nh;
    double zero = 0.0, one = 1.0;
    double *dnmI, *qnmI;
    double *temA;
    double *dnm, *qnm;

    dnm = new double[num_proj * num_proj];
    qnm = new double[num_proj * num_proj];
    temA = new double[num_orb * num_proj];

    for(int idx = 0; idx < num_proj * num_proj; idx++) 
    {
        dnm[idx] = 0.0;
        qnm[idx] = 0.0;
    }

    int proj_count = 0;
    for (int ion = 0; ion < ct.num_ions; ion++)
    {
        ION *iptr = &Atoms[ion];
        SPECIES *sp = &ct.sp[iptr->species];

        int nh = sp->num_projectors;

        if(nh == 0) continue;
        dnmI = pct.dnmI[ion];
        qnmI = pct.qqq[ion];

        for(int i = 0; i < nh; i++)
        for(int j = 0; j < nh; j++)
        {
            int ii = i + proj_count;
            int jj = j + proj_count;
            dnm[ii *num_proj + jj] = dnmI[i * nh + j];
            qnm[ii *num_proj + jj] = qnmI[i * nh + j];
        }

        proj_count += nh;
    }
           
    

    dgemm ("T", "N", &num_orb, &num_proj, &num_proj, &one, Kbpsi_mat, &num_proj, dnm, &num_proj, &zero, temA, &num_orb);
    dgemm ("N", "N", &num_orb, &num_orb, &num_proj, &one, temA, &num_orb, Kbpsi_mat, &num_proj, &one, Aij, &num_orb);


    dgemm ("T", "N", &num_orb, &num_proj, &num_proj, &one, Kbpsi_mat, &num_proj, qnm, &num_proj, &zero, temA, &num_orb);
    dgemm ("N", "N", &num_orb, &num_orb, &num_proj, &one, temA, &num_orb, Kbpsi_mat, &num_proj, &one, Bij, &num_orb);



    delete [] temA;
    delete [] qnm;
    delete [] dnm;
    
}
