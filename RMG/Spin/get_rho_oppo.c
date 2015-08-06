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




void  get_rho_oppo (double * rho, double * rho_oppo)
{
    MPI_Status status;

    MPI_Sendrecv(rho,(int) get_FP0_BASIS(), MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, 
	    rho_oppo,(int) get_FP0_BASIS(), MPI_DOUBLE, (pct.spinpe+1)%2, pct.gridpe, pct.spin_comm, &status);
}                               /* end scf */


/******/
