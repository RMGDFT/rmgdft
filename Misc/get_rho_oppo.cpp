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



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "common_prototypes.h"
#include "transition.h"




void  get_rho_oppo (double * rho, double * rho_oppo)
{
    if(ct.AFM)
    {
         Rmg_Symm->symmetrize_rho_AFM(rho, &rho[get_FP0_BASIS()]);
    }
    else
    {
        MPI_Status status;
        int tag = pct.kpsub_rank * pct.grid_npes + pct.gridpe;
        int send_tag = pct.spinpe * pct.pe_kpoint * pct.grid_npes + tag;
        int recv_tag = (pct.spinpe+1)%2 * pct.pe_kpoint * pct.grid_npes + tag;
        MPI_Sendrecv(rho,(int) get_FP0_BASIS(), MPI_DOUBLE, (pct.spinpe+1)%2, send_tag, 
                rho_oppo,(int) get_FP0_BASIS(), MPI_DOUBLE, (pct.spinpe+1)%2, recv_tag, pct.spin_comm, &status);
        MPI_Barrier(pct.img_comm);
    }
}                               /* end scf */


/******/
