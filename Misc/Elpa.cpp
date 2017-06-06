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

#include "portability.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>


#if USE_ELPA
#include "RmgException.h"
#include "Elpa.h"
#include "blacs.h"


void set_up_blacsgrid_from_fortran(int mpi_comm_world, int* my_blacs_ctxt,
                                    int *np_rows, int *np_cols, int *nprow, int *npcol,
                                    int *my_prow, int *my_pcol);


// Clean up. Elpa inherits scalapack so scalapack destructor is called automatically on exit
Elpa::~Elpa(void)
{
}

// Gets row and column communicators needed by elpa routines
void Elpa::GetCommunicators(void)
{
    int mpierr;
    int np_rows, np_cols;
    int np_row, np_col;
    int my_prow, my_pcol;
    int my_blacs_ctxt;

//    set_up_blacsgrid_from_fortran(MPI_Comm_c2f(this->comm), &my_blacs_ctxt, &this->group_rows, &this->group_cols, &np_row, &np_col, &this->my_row, &this->my_col);
    mpierr = elpa_get_communicators(MPI_Comm_c2f(this->root_comm), this->my_row, this->my_col, &this->elpa_comm_rows, &this->elpa_comm_cols);
    printf("MPIERR = %d  %d  %d  %d  %d\n",mpierr, this->my_row, this->my_col, this->elpa_comm_rows, this->elpa_comm_cols);
}

int Elpa::GetElpaCommRows(void)
{
    return this->elpa_comm_rows;
}

int Elpa::GetElpaCommCols(void)
{
    return this->elpa_comm_cols;
}
#endif
