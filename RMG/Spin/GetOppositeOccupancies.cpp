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


template void GetOppositeOccupancies(Kpoint<double> **);
template void GetOppositeOccupancies(Kpoint<std::complex<double> > **);

template <typename KpointType>
void GetOppositeOccupancies (Kpoint<KpointType> ** Kptr)
{

    MPI_Status status; 
    double *occ_sd, *occ_rv;   
    int st, st1, kpt;


    /* allocate memory for occupation send array and receive array */
    occ_sd = new double[2 * ct.num_kpts_pe * ct.num_states];
    occ_rv = occ_sd + ct.num_kpts_pe * ct.num_states;


    /*Prepare the sending buffer of occupancies */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            occ_sd[st] = Kptr[kpt]->Kstates[st1].occupation[0];
            st += 1;
        }	
    }


    /*Communicate for spin up and spin down occupations*/    
    MPI_Sendrecv(occ_sd, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe,
            occ_rv, st, MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);


    /* Unpack the received occupations to state structure */
    st = 0;
    for (kpt =0; kpt < ct.num_kpts_pe; kpt++)
    {
        for (st1 = 0; st1 < ct.num_states; st1++)
        {	
            Kptr[kpt]->Kstates[st1].occupation[1] = occ_rv[st];
            st += 1;

        }	
    } 


    delete [] occ_sd;
}                               /* end scf */


/******/
