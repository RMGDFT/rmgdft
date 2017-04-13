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


#include "RmgException.h"
#include "Elpa.h"
#include "blacs.h"




// Clean up
Elpa::~Elpa(void)
{
    Cblacs_gridexit(this->context);
    MPI_Comm_free(&this->comm);
    delete [] this->local_desca;
    delete [] this->dist_desca;
}

// Gets row and column communicators needed by elpa routines
void Elpa::GetCommunicators(void)
{
    elpa_get_communicators(MPI_Comm_c2f(this->comm), this->my_row, this->my_col, &this->elpa_comm_rows, &this->elpa_comm_cols);
}

int Elpa::GetElpaCommRows(void)
{
    return this->elpa_comm_rows;
}

int Elpa::GetElpaCommCols(void)
{
    return this->elpa_comm_cols;
}
