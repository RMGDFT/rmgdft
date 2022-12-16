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

#include <vector>
#include <algorithm>
#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>


#if USE_ELPA
#include "RmgException.h"
#include <elpa/elpa.h>
#include "Elpa.h"
#include "blacs.h"

static int once;
void Elpa::Init(void)
{
    int error;
    int api_version = 20221109;
    if(!once)
    {
        elpa_init(api_version);
        once = 1;
    }
    this->handle = (elpa_t *)elpa_allocate(&error);
    if(error != ELPA_OK)
    {
        throw RmgFatalException() << "Error initializeing elpa in "  << __FILE__ << " at line " << __LINE__ << ".\n";

    }

    elpa_set((elpa_t)this->handle, "na", this->N, &error);
    elpa_set((elpa_t)this->handle, "nev", this->N, &error);
    elpa_set((elpa_t)this->handle, "local_nrows", this->m_dist, &error);
    elpa_set((elpa_t)this->handle, "local_ncols", this->n_dist, &error);
    elpa_set((elpa_t)this->handle, "nblk", this->NB, &error);
    elpa_set((elpa_t)this->handle, "mpi_comm_parent", MPI_Comm_c2f(this->comm), &error);
    elpa_set((elpa_t)this->handle, "process_row", this->my_row, &error);
    elpa_set((elpa_t)this->handle, "process_col", this->my_col, &error);
    elpa_set((elpa_t)this->handle, "num_process_rows", this->group_rows, &error);
    elpa_set((elpa_t)this->handle, "num_process_cols", this->group_cols, &error);
    elpa_set((elpa_t)this->handle, "blacs_context", this->context, &error);
    elpa_setup((elpa_t)this->handle);
}

// Clean up. Elpa inherits scalapack so scalapack destructor is called automatically on exit
Elpa::~Elpa(void)
{
    int error;
    elpa_deallocate((elpa_t)this->handle, &error);
    //elpa_uninit(&error);
}

void Elpa::generalized_eigenvectors(double *a, double *b, double *ev, double *q, int is_already_decomposed, int *error)
{
//    elpa_set((elpa_t)this->handle, "solver", ELPA_SOLVER_2STAGE, error);
//    elpa_set((elpa_t)this->handle, "real_kernel", ELPA_2STAGE_REAL_AMD_GPU, error);
    elpa_generalized_eigenvectors((elpa_t)this->handle, a, b, ev, q, is_already_decomposed, error);
}

void Elpa::generalized_eigenvectors(std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, int is_already_decomposed, int *error)
{
    elpa_generalized_eigenvectors((elpa_t)this->handle, a, b, ev, q, is_already_decomposed, error);
}

#if 0
int Elpa::GetElpaCommRows(void)
{
    return this->elpa_comm_rows;
}

int Elpa::GetElpaCommCols(void)
{
    return this->elpa_comm_cols;
}
#endif
#endif
