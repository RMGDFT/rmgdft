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
#include <math.h>
#include <complex>

#include "RmgException.h"
#include "Scalapack.h"
#include "blacs.h"

Scalapack::Scalapack(int ngroups, int thisimg, int images_per_node, int M, int N, int MB, int NB, MPI_Comm rootcomm)
{

    this->ngroups = ngroups;
    this->msize = msize;
    this->M = M;
    this->N = N;
    this->MB = MB;
    this->NB = NB;
    this->root_comm = rootcomm;

    MPI_Comm_size(rootcomm, &this->npes);
    MPI_Comm_rank(rootcomm, &this->root_rank);

    MPI_Group grp_world, grp_this;


    // Make sure we have enough pes for the number of groups requested. We'll reset group_pes later
    this->group_pes = this->npes / this->ngroups / images_per_node;
    if(this->group_pes < 1) 
        throw RmgFatalException() << "Too many Scalapack groups requested in " << __FILE__ << " at line " << __LINE__ << ".\n";

    this->group_sizes = new int[this->ngroups]();
    this->group_starts = new int[this->ngroups]();
    int j = 0;
    for(int i = 0;i < this->npes;i++) {
        this->group_sizes[j]++;
        j++;
        j = j % this->ngroups; 
    }

    this->group_index = 0;
    for(int i = 1;i < this->ngroups;i++) {
        this->group_starts[i] = this->group_starts[i-1] + this->group_sizes[i-1]; 
        if(this->root_rank >= this->group_starts[i]) this->group_index++;
    }
    this->group_pes = this->group_sizes[this->group_index];

    // Get 2-d grid dimensions for this group
    int sqrtnpe = (int) (sqrt (this->group_pes)) + 1;
    for (int i = 1; i <= sqrtnpe; i++) {
        if (this->group_pes % i == 0) this->group_rows = i;
    }
    this->group_cols = this->group_pes / this->group_rows;

    int *pmap, *tgmap;

    pmap = new int[this->npes];
    tgmap = new int[this->npes];



    Cblacs_get (0, 0, &this->context);

    // We need to setup the MPI world rank range in this group to be mapped to blacs
    for (int i = 0; i < this->npes; i++) tgmap[i] = i;

    MPI_Comm_split(rootcomm, this->group_index + 1, 0, &this->comm);

    MPI_Comm_group (MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group (this->comm, &grp_this);
    MPI_Group_translate_ranks (grp_this, this->group_pes, tgmap, grp_world, pmap);
    MPI_Comm_rank(this->comm, &this->comm_rank);


    /* Assign nprow*npcol processes to blacs for calculations */
    int item;
    item = thisimg % images_per_node;

    Cblacs_gridmap (&this->context, &pmap[item * this->group_rows * this->group_cols], this->group_rows, this->group_rows, this->group_cols);

    // Figures out blacs information, used so that we can find 
    // which processors are participating in scalapack operations
    Cblacs_gridinfo (this->context, &this->group_rows, &this->group_cols, &this->my_row, &this->my_col);

//    printf("\n  myrow, mycol nprow npcol %d %d %d %d", this->my_row, this->my_col, this->group_rows, this->group_cols);

    // Variable participates will tell use whether the PE participates in scalapack calculations
    this->participates = (this->my_row >= 0);
//std::cout << "JJJ " << this->group_index << "  " <<this->my_row << "  " << this->my_col << "  " << this->participates << std::endl;
    #if 0

        if(pct.scalapack_pe) {
            pct.scalapack_mpi_rank[myrow*npcol + mycol] = pct.gridpe;
        }
        MPI_Allreduce(MPI_IN_PLACE, pct.scalapack_mpi_rank, npes, MPI_INT, MPI_SUM, pct.grid_comm);

        pct.scalapack_max_dist_size = numroc_ (&ct.num_states, &NB, &myrow, &izero, &nprow) *
                                      numroc_ (&ct.num_states, &NB, &mycol, &izero, &npcol);

        MPI_Allreduce(MPI_IN_PLACE, &pct.scalapack_max_dist_size, 1, MPI_INT, MPI_MAX, pct.grid_comm);
    #endif


    // Set up descriptors for a local matrix (one block on (0,0)
    int izero = 0, info = 0;
    int lld = std::max( numroc_( &this->N, &this->N, &this->my_row, &izero, &this->group_rows ), 1 );
    this->local_desca = new int[this->npes*DLEN];
    descinit_( this->local_desca, &this->M, &this->N, &this->M, &this->N, &izero, &izero, &this->context, &lld, &info );

    // Get dimensions of the local part of the distributed matrix
    this->m_dist = numroc_( &this->M, &this->MB, &this->my_row, &izero, &this->group_rows );
    this->n_dist = numroc_( &this->M, &this->NB, &this->my_col, &izero, &this->group_cols );
    
    // descriptors for the distributed matrix 
    int lld_distr = std::max( this->m_dist, 1 );
    this->dist_desca = new int[this->npes*DLEN];
    descinit_( this->dist_desca, &this->M, &this->N, &this->MB, &this->NB, &izero, &izero, &this->context, &lld_distr, &info );



    delete [] tgmap;
    delete [] pmap;

}

int Scalapack::GetRootRank(void)
{
    return this->root_rank;
}

int Scalapack::GetCommRank(void)
{
    return this->comm_rank;
}

int Scalapack::GetDistMdim(void)
{
    return this->m_dist;
}

int Scalapack::GetDistNdim(void)
{
    return this->n_dist;
}

int *Scalapack::GetDistDesca(void)
{
    return this->dist_desca;
}

// Returns ipiv size required by PDGESV
int Scalapack::GetIpivSize(void)
{
    return numroc_ (&this->dist_desca[2], &this->dist_desca[4], &this->my_row, &this->dist_desca[6],
                                &this->group_rows) + this->dist_desca[4];

}

void Scalapack::DistributeMatrix(double *A, double *A_dist, int m, int n)
{
    // Check that this scalapack instance matches the specified matrix sizes
    if((m != this->M) || (n != this->N)) {
        throw RmgFatalException() << "Error: Scalapack instance was set up with (M,N)=(" << this->M << "," << this->N << ") but request was for (M,N)=(" << m << "," << n << "). Terminating.\n";
    }

    // Call pdgeadd_ to distribute matrix (i.e. copy A into A_dist)
    int ione = 1;
    double rone = 1.0, rzero = 0.0;
    pdgeadd_( "N", &this->M, &this->N, &rone, A, &ione, &ione, this->local_desca, &rzero, A_dist, &ione, &ione, this->dist_desca );
}


void Scalapack::GatherMatrix(double *A, double *A_dist, int m, int n)
{

    // Check that this scalapack instance matches the specified matrix sizes
    if((m != this->M) || (n != this->N)) {
        throw RmgFatalException() << "Error: Scalapack instance was set up with (M,N)=(" << this->M << "," << this->N << ") but request was for (M,N)=(" << m << "," << n << "). Terminating.\n";
    }

    // Call pdgeadd_ to gather matrix (i.e. copy A_dist into A)
    int ione = 1;
    double rone = 1.0, rzero = 0.0;
    pdgeadd_( "N", &this->M, &this->N, &rone, A_dist, &ione, &ione, this->dist_desca, &rzero, A, &ione, &ione, this->local_desca );

}


void Scalapack::DistributeMatrix(std::complex<double> *A, std::complex<double> *A_dist, int m, int n)
{
    // Check that this scalapack instance matches the specified matrix sizes
    if((m != this->M) || (n != this->N)) {
        throw RmgFatalException() << "Error: Scalapack instance was set up with (M,N)=(" << this->M << "," << this->N << ") but request was for (M,N)=(" << m << "," << n << "). Terminating.\n";
    }

    // Call pdgeadd_ to distribute matrix (i.e. copy A into A_dist)
    int ione = 1;
    double rone = 1.0, rzero = 0.0;
    pzgeadd_( "N", &this->M, &this->N, &rone, (double *)A, &ione, &ione, this->local_desca, &rzero, (double *)A_dist, &ione, &ione, this->dist_desca );
}


void Scalapack::GatherMatrix(std::complex<double> *A, std::complex<double> *A_dist, int m, int n)
{

    // Check that this scalapack instance matches the specified matrix sizes
    if((m != this->M) || (n != this->N)) {
        throw RmgFatalException() << "Error: Scalapack instance was set up with (M,N)=(" << this->M << "," << this->N << ") but request was for (M,N)=(" << m << "," << n << "). Terminating.\n";
    }

    // Call pdgeadd_ to gather matrix (i.e. copy A_dist into A)
    int ione = 1;
    double rone = 1.0, rzero = 0.0;
    pzgeadd_( "N", &this->M, &this->N, &rone, (double *)A_dist, &ione, &ione, this->dist_desca, &rzero, (double *)A, &ione, &ione, this->local_desca );

}



// Clean up
Scalapack::~Scalapack(void)
{
    Cblacs_gridexit(this->context);
    MPI_Comm_free(&this->comm);
    delete [] this->group_sizes;
    delete [] this->group_starts;
    delete [] this->local_desca;
    delete [] this->dist_desca;
}
