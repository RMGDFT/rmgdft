/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/


#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


template void GetOppositeEigvals(Kpoint<double> **);
template void GetOppositeEigvals(Kpoint<std::complex<double> > **);

template <typename KpointType>
void GetOppositeEigvals (Kpoint<KpointType> ** Kptr)
{

    MPI_Status status; 
    double *eigval_sd, *eigval_rv;   
    int st, st1, kpt;


    /* allocate memory for eigenvalue send array and receive array */
    eigval_sd = new double[2 * ct.num_kpts_pe * ct.num_states];
    eigval_rv = eigval_sd + ct.num_kpts_pe * ct.num_states;


    /*Prepare the sending buffer of eigenvalues */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            eigval_sd[st] = Kptr[kpt]->Kstates[st1].eig[0];
            st += 1;
        }	
    }


    /*Communicate for spin up and spin down energy eigenvalues*/    
    MPI_Sendrecv(eigval_sd, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
            eigval_rv, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);


    /* Unpack the received eigenvalue to state structure */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            Kptr[kpt]->Kstates[st1].eig[1] = eigval_rv[st];
            st += 1;

        }	
    } 

    /*Prepare the sending buffer for the feig values (used for the scf energy correction term */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            eigval_sd[st] = Kptr[kpt]->Kstates[st1].feig[0];
            st += 1;
        }	
    }


    /*Communicate for spin up and spin down energy eigenvalues*/    
    MPI_Sendrecv(eigval_sd, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
            eigval_rv, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);


    /* Unpack the received eigenvalue to state structure */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            Kptr[kpt]->Kstates[st1].feig[1] = eigval_rv[st];
            st += 1;

        }	
    } 


    delete [] eigval_sd;
}                               /* end scf */


