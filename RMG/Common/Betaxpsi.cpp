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


#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "GpuAlloc.h"
#include "RmgGemm.h"
#include "blas.h"
#include "../Headers/prototypes.h"


#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template <typename KpointType>
void betaxpsi_calculate (Kpoint<KpointType> * sint_ptr, KpointType * psi, int num_states, KpointType *);

template <typename KpointType>
void betaxpsi_receive (KpointType * recv_buff, int num_pes,
                               int *pe_list, int *num_ions_per_pe,
                               MPI_Request * req_recv, int num_states);

template <typename KpointType>
void betaxpsi_send (KpointType * send_buff, int num_pes,
                            int *pe_list, int *num_ions_per_pe,
                            MPI_Request * req_send, int num_states);

template <typename KpointType>
void betaxpsi_pack (KpointType * sint, KpointType * fill_buff,
                            int num_pes, int *num_ions_per_pe,
                            int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS], int num_states);

template <typename KpointType>
void betaxpsi_sum_owned (KpointType * recv_buff, KpointType * sint,
                                 int num_pes, int *num_ions_per_pe,
                                 int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS], int num_states);

template <typename KpointType>
void betaxpsi_write_non_owned (KpointType * sint, KpointType * recv_buff,
                                       int num_pes,
                                       int *num_ions_per_pe,
                                       int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS], int num_states);

template void Betaxpsi<double>(Kpoint<double> *, int, int, double *, double *);
template void Betaxpsi<std::complex<double> >(Kpoint<std::complex<double>> *, int, int, std::complex<double> *, std::complex<double> *);


template <typename KpointType>
void Betaxpsi (Kpoint<KpointType> *kptr, int first_state, int nstates, KpointType *sint_local, KpointType *nlweight)
{


    if((ct.localize_projectors == false) || (Rmg_G->get_NPES() == 1))
    {

        /*Loop over ions and calculate local projection between beta functions and wave functions */
        int factor = 1;
        KpointType rzero(0.0);
        KpointType alpha(get_vel());
        KpointType *NULLptr = NULL;

        if(typeid(KpointType) == typeid(std::complex<double>)) factor = 2;
        char *transt = "t", *transn = "n", *transc = "c";
        char *transa;
        transa = transc;
        if(ct.is_gamma) transa = transt;

        int length = factor * ct.num_ions * nstates * ct.max_nl;
        RmgGemm (transa, transn, pct.num_tot_proj, nstates, kptr->pbasis, alpha,
            nlweight, kptr->pbasis, &kptr->orbital_storage[first_state*kptr->pbasis], kptr->pbasis,
            rzero, sint_local, pct.num_tot_proj, NULLptr, NULLptr, NULLptr, false, false, false, true);


        if(Rmg_G->get_NPES() != 1)
            MPI_Allreduce(MPI_IN_PLACE, (double *)sint_local, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        return;
    }


    KpointType *own_buff = NULL, *nown_buff = NULL;
    KpointType *send_buff, *recv_buff;
    MPI_Request *req_send, *req_recv;
    MPI_Request *req_own = NULL, *req_nown = NULL;
    int size_own, size_nown, pe;

    /*Allocate memory for communication */
    size_own = 0;
    for (pe = 0; pe < pct.num_owned_pe; pe++)
        size_own += pct.num_owned_ions_per_pe[pe];
    size_own *= nstates * ct.max_nl;

    size_nown = 0;
    for (pe = 0; pe < pct.num_owners; pe++)
        size_nown += pct.num_nonowned_ions_per_pe[pe];
    size_nown *= nstates * ct.max_nl;


    if (size_own) {
        own_buff = new KpointType[size_own]();
    }
    if (size_nown) {
        nown_buff = new KpointType[size_nown]();
    }

    if (pct.num_owned_pe) {
        req_own = new MPI_Request[pct.num_owned_pe];
    }
    if (pct.num_owners) {
        req_nown = new MPI_Request[pct.num_owners];
    }

    /*First owning cores will receive data from non-owners */
    send_buff = nown_buff;
    recv_buff = own_buff;

    req_send = req_nown;
    req_recv = req_own;

    /*First post non-blocking receives for data about owned ions from cores who do not own the ions */
    betaxpsi_receive (recv_buff, pct.num_owned_pe, pct.owned_pe_list,
                       pct.num_owned_ions_per_pe, req_recv, nstates);


    // Set sint array
    //sint = &kptr->newsint_local[first_state*pct.num_nonloc_ions*ct.max_nl];
    KpointType *sint = new KpointType[pct.num_nonloc_ions * nstates * ct.max_nl]();


    /*Loop over ions and calculate local projection between beta functions and wave functions */
    betaxpsi_calculate (kptr, sint, &kptr->orbital_storage[first_state*kptr->pbasis], nstates, nlweight);


    /*Pack data for sending */
    betaxpsi_pack (sint, send_buff, pct.num_owners,
                    pct.num_nonowned_ions_per_pe, pct.list_ions_per_owner, nstates);

    /*Send <beta|psi> contributions  to the owning PE */
    betaxpsi_send (send_buff, pct.num_owners, pct.owners_list,
                    pct.num_nonowned_ions_per_pe, req_send, nstates);

    /*Wait until all data is received */
    if(pct.num_owned_pe)
        MPI_Waitall (pct.num_owned_pe, req_recv, MPI_STATUSES_IGNORE);
    if(pct.num_owners)
        MPI_Waitall (pct.num_owners, req_send, MPI_STATUSES_IGNORE);


    /*Unpack received data and sum contributions from all pes for owned ions */
    betaxpsi_sum_owned (recv_buff, sint, pct.num_owned_pe,
                         pct.num_owned_ions_per_pe, pct.list_owned_ions_per_pe, nstates);

    /*In the second stage, owning cores will send summed data to non-owners */
    send_buff = own_buff;
    recv_buff = nown_buff;

    req_send = req_own;
    req_recv = req_nown;
    
    
    /*Receive summed data for non-owned ions from owners */
    betaxpsi_receive (recv_buff, pct.num_owners, pct.owners_list,
                       pct.num_nonowned_ions_per_pe, req_recv, nstates);

    /*Pack summed data for owned ions to send to non-owners */
    betaxpsi_pack (sint, send_buff, pct.num_owned_pe,
                    pct.num_owned_ions_per_pe, pct.list_owned_ions_per_pe, nstates);

    /*Send packed data for owned ions to non-owners */
    betaxpsi_send (send_buff, pct.num_owned_pe, pct.owned_pe_list,
                    pct.num_owned_ions_per_pe, req_send, nstates);

    /*Wait until all data is received */
    if(pct.num_owned_pe)
        MPI_Waitall (pct.num_owned_pe, req_send, MPI_STATUSES_IGNORE);
    if(pct.num_owners)
        MPI_Waitall (pct.num_owners, req_recv, MPI_STATUSES_IGNORE);

    /*Finaly, write received data about non-owned ions into sint array */
    betaxpsi_write_non_owned (sint, recv_buff, pct.num_owners,
                               pct.num_nonowned_ions_per_pe, pct.list_ions_per_owner, nstates);


    if (pct.num_owned_pe)
        delete [] req_own;
    if (pct.num_owners)
        delete [] req_nown;
    if (size_nown)
        delete [] nown_buff;
    if (size_own)
        delete [] own_buff;


    int idx= 0 ;
//    KpointType *tsint = &sint_local[first_state*pct.num_nonloc_ions*ct.max_nl];
    KpointType *tsint = sint_local;
    for(int st = 0;st < nstates;st++) {
        for(int ion = 0;ion < pct.num_nonloc_ions;ion++) {
            for(int ip = 0;ip < ct.max_nl;ip++) {
                tsint[idx] = sint[ion*nstates*ct.max_nl + st*ct.max_nl + ip];
                idx++;
            }
        }
    }
    delete [] sint;
}



template <typename KpointType>
void betaxpsi_calculate (Kpoint<KpointType> *kptr, KpointType * sint_ptr, KpointType *psi, int num_states, KpointType *nlweight)
{
    KpointType rzero(0.0);
    KpointType *NULLptr = NULL;
    char *transt = "t", *transn = "n", *transc = "c";
    char *transa;
   
    transa = transc;
    if(ct.is_gamma) transa = transt;

    KpointType alpha(get_vel());

    if(pct.num_tot_proj == 0) return;
    int pbasis = kptr->pbasis;

#if GPU_ENABLED
    KpointType *nlarray = (KpointType *)GpuMallocHost(sizeof(KpointType) * pct.num_tot_proj * num_states);
#else
    KpointType *nlarray = new KpointType[pct.num_tot_proj * num_states]();
#endif
    RmgGemm (transa, transn, pct.num_tot_proj, num_states, pbasis, alpha, 
            nlweight, pbasis, psi, pbasis, 
            rzero, nlarray, pct.num_tot_proj, NULLptr, NULLptr, NULLptr, false, false, false, true);

    for (int nion = 0; nion < pct.num_nonloc_ions; nion++)
    {
        for(int istate = 0; istate < num_states; istate++)
        {
            for(int ip = 0; ip < ct.max_nl; ip++)
            {
                int proj_index = nion * ct.max_nl + ip;
                int idx1 = nion * num_states * ct.max_nl + istate * ct.max_nl + ip;
                //idx2 = proj_index * num_states + istate;
                int idx2 = istate * pct.num_tot_proj + proj_index;
                sint_ptr[idx1] = nlarray[idx2];
            }
        }
    }

#if GPU_ENABLED
    GpuFreeHost(nlarray);
#else
    delete [] nlarray;
#endif

}



/*This receives data from other PEs for ions owned by current PE*/
template <typename KpointType>
void betaxpsi_receive (KpointType * recv_buff, int num_pes,
        int *pe_list, int *num_ions_per_pe,
        MPI_Request * req_recv, int num_states)
{
    KpointType *tpr;
    int tag, pe, source, size;

    tpr = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        source = pe_list[pe];
        /*Tag is sender_rank * NPES + receiver_rank */
        tag = 111;
        size = num_ions_per_pe[pe] * num_states * ct.max_nl;
        int transfer_size = size;
        if(!ct.is_gamma) transfer_size *= 2;

        MPI_Irecv (tpr, transfer_size, MPI_DOUBLE, source, tag, pct.grid_comm, &req_recv[pe]);

        tpr += size;
    }

}

template <typename KpointType>
void betaxpsi_send (KpointType * send_buff, int num_pes,
        int *pe_list, int *num_ions_per_pe,
        MPI_Request * req_send, int num_states)
{
    KpointType *tpr;
    int target, num_ions, size, tag, pe;

    tpr = send_buff;

    for (pe = 0; pe < num_pes; pe++)
    {

        target = pe_list[pe];
        /*Tag is sender_rank * NPES + receiver_rank */
        tag = 111;
        num_ions = num_ions_per_pe[pe];
        size = num_ions * num_states * ct.max_nl;
        int transfer_size = size;
        if(!ct.is_gamma) transfer_size *= 2;

        MPI_Isend (tpr, transfer_size, MPI_DOUBLE, target, tag, pct.grid_comm, &req_send[pe]);

        tpr += size;
    }

}


template <typename KpointType>
void betaxpsi_pack (KpointType * sint, KpointType * fill_buff, 
        int num_pes, int *num_ions_per_pe,
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS], int num_states)
{
    KpointType *tpr_buff, *sint_tpr;
    int size, num_ions, ion, nlion, pe;

    tpr_buff = fill_buff;

    size = num_states * ct.max_nl;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        /*Loop over ions that need to be sent to pe */
        for (ion = 0; ion < num_ions; ion++)
        {
            nlion = list_ions_per_pe[pe][ion];

            sint_tpr = &sint[nlion * num_states * ct.max_nl];

            /*Pack data into send array */
            for(int idx=0;idx < size;idx++) tpr_buff[idx] = sint_tpr[idx];
            tpr_buff += size;

        }

    }

}


template <typename KpointType>
void betaxpsi_sum_owned (KpointType * recv_buff, KpointType * sint,
        int num_pes, int *num_ions_per_pe,
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS], int num_states)
{
    KpointType *tpr1, *tpr2;
    int size, num_ions, ion_index, pe, ion;

    size = num_states * ct.max_nl;
    tpr1 = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sint[ion_index * num_states * ct.max_nl];
            for(int idx=0;idx < size;idx++) tpr2[idx] += tpr1[idx];
            //my_axpy (1.0, tpr1, tpr2, size);

            tpr1 += size;
        }
    }


}



template <typename KpointType>
void betaxpsi_write_non_owned (KpointType * sint, KpointType * recv_buff, int num_pes,
        int *num_ions_per_pe,
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS], int num_states)
{
    KpointType *tpr1, *tpr2;
    int size, num_ions, ion_index, pe, ion;

    size = num_states * ct.max_nl;
    tpr1 = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sint[ion_index * num_states * ct.max_nl];
            for(int idx=0;idx < size;idx++) tpr2[idx] = tpr1[idx];

            tpr1 += size;
        }
    }


}

