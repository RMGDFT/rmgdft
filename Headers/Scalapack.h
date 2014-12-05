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


#ifndef RMG_Scalapack_H
#define RMG_Scalapack_H 1

#include <mpi.h>


/* Blacs dimension */
#define DLEN    9

// npes is the total number of pes, ngroups is the number of
// scalapack groups that should be created out of npes where
// npes <= group_pes * ngroups * images_per_node

class Scalapack {

public:

    Scalapack(int ngroups, int thisimg, int images_per_node, int block_factor, MPI_Comm rootcomm);

private:

    int context;        // blacs context of this group of pes
    int npes;           // total number of pes
    int ngroups;        // total number of groups
    int group_pes;      // number of pes in this group
    int group_index;    // index of this group
    int group_rows;     // rows in this group
    int group_cols;     // cols in this group
    int my_row;         // blacs row of this PE
    int my_col;         // blacs col of this PE
    int root_rank;      // rank in rootcomm
    int comm_rank;      // rank in this comm
    int *group_sizes;   // number of pes in this group
    int *group_starts;  // starting rank in rootrank assigned to this group
    int *desca;         // scalapack desca for this group
    bool participates;  // whether or not this PE participates in scalapack calculations
    MPI_Comm comm;      // communicator for this object

};



#endif
