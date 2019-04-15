/*
 *
 * Copyright (c) 2019, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#ifndef RMG_Projector_H
#define RMG_Projector_H 1

#include <boost/multi_array.hpp>
#include <cassert>

// This class is used to manage projector objects where a projector is defined by
// 
//     Pij = <P_i | psi_j>
//
//     where P_i is the ith projector and psi_j is the jth wavefunction
//
// Projector objects may include beta functions, atomic orbitals and spherical harmonics.
//


template <typename KpointType> class Projector {

public:
    typedef boost::multi_array<int, 2> int_2d_array;
    Projector(Kpoint<KpointType> *K, int projector_type, int num_pes, int num_ions);
    ~Projector(void);
    void set_storage(KpointType *storage);
    void project(KpointType *p, int offset, int n);
    int get_num_nonloc_ions(void);
    int get_num_owned_ions(void);
    int *get_owned_ions_list(void);
    int *get_owned_pe_list(void);
    int *get_nonloc_ions_list(void);

    // kpoint this Projector is associated with
    Kpoint<KpointType> *kptr;

    // Number of projectors
    int num_tot_proj;

    // Storage area for the projector
    KpointType *weight;
 
    // basis size for an individual projector
    int pbasis;

    // Type LOCALIZED or DELOCALIZED
    int type;

    int *idxptrlen;
    double xcstart;
    double ycstart;
    double zcstart;


private:
    int num_owned_ions;
    int owned_ions_list[MAX_NONLOC_IONS];

    int num_nonloc_ions;
    int *nonloc_ions_list;
    bool *nonloc_ion_ownflag;
    int num_nonloc_pes;

    /*For ions owned by current PE */
    /* Number of cores to cummunicate with about owned ions*/
    int num_owned_pe;
    /*List of ranks of cores to comunicate with about owned ions*/
    int *owned_pe_list;
    /* Number of owned ions to communicate about for cores from owned_pe_list */
    int *num_owned_ions_per_pe;
    /*List of ion indices to communicate about for core from owned_pe_list 
     * These are indices within nonloc ions, not absolute ones*/
    //int list_owned_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS];
    int_2d_array list_owned_ions_per_pe;

    /*For ions NOT owned by current PE*/
    /* Number of cores to cummunicate with about non-owned ions*/
    int num_owners;
    /*Indices of cores to cummunicate with about non-owned ions*/
    int *owners_list;
    /* Number of non-owned ions to communicate about for cores from owners_list  */
    int *num_nonowned_ions_per_pe;
    /*List of ion indices to communicate about for cores from owners_list
     * These are indices within nonloc ions, not absolute ones*/
    //int list_ions_per_owner[MAX_NONLOC_PROCS][MAX_NONLOC_IONS];
    int_2d_array list_ions_per_owner;

    int num_loc_ions;
    int loc_ions_list[MAX_NONLOC_IONS];

    void betaxpsi_calculate (Kpoint<KpointType> * kptr, KpointType * sint_ptr, KpointType * psi, int num_states);
    void betaxpsi_receive (KpointType * recv_buff, int num_pes,
                               int *pe_list, int *num_ions_per_pe,
                               MPI_Request * req_recv, int num_states);
    void betaxpsi_send (KpointType * send_buff, int num_pes,
                            int *pe_list, int *num_ions_per_pe,
                            MPI_Request * req_send, int num_states);
    void betaxpsi_pack (KpointType * sint, KpointType * fill_buff,
                            int num_pes, int *num_ions_per_pe,
                            int_2d_array &list_ions_per_pe, int num_states);
    void betaxpsi_sum_owned (KpointType * recv_buff, KpointType * sint,
                                 int num_pes, int *num_ions_per_pe,
                                 int_2d_array &list_ions_per_pe, int num_states);
    void betaxpsi_write_non_owned (KpointType * sint, KpointType * recv_buff,
                                       int num_pes,
                                       int *num_ions_per_pe,
                                       int_2d_array &list_ions_per_pe, int num_states);

};

#endif
