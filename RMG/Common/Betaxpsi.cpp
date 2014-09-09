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
#include "RmgGemm.h"
#include "blas.h"
#include "../Headers/prototypes.h"

// Turning this off for now while we are in transition
#undef GPU_ENABLED
#define GPU_ENABLED 0

#if GPU_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#endif


template <typename KpointType>
void betaxpsi_calculate (Kpoint<KpointType> * sint_ptr, KpointType * psi);

template <typename KpointType>
void betaxpsi_receive (KpointType * recv_buff, int num_pes,
                               int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
                               MPI_Request * req_recv);

template <typename KpointType>
void betaxpsi_send (KpointType * send_buff, int num_pes,
                            int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
                            MPI_Request * req_send);

template <typename KpointType>
void betaxpsi_pack (KpointType * sint, KpointType * fill_buff,
                            int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
                            int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS]);

template <typename KpointType>
void betaxpsi_sum_onwed (KpointType * recv_buff, KpointType * sint,
                                 int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
                                 int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS]);

template <typename KpointType>
void betaxpsi_write_non_owned (KpointType * sint, KpointType * recv_buff,
                                       int num_pes,
                                       int num_ions_per_pe[MAX_NONLOC_PROCS],
                                       int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS]);

template void Betaxpsi<double>(Kpoint<double> *);
template void Betaxpsi<std::complex<double> >(Kpoint<std::complex<double>> *);

template <typename KpointType>
void Betaxpsi (Kpoint<KpointType> *kptr)
{
    KpointType ZERO_t(0.0);

    KpointType *own_buff = NULL, *nown_buff = NULL, *sint = NULL;
    KpointType *send_buff, *recv_buff;
    MPI_Request *req_send, *req_recv;
    MPI_Request *req_own = NULL, *req_nown = NULL;
    int size_own, size_nown, pe, i;

    /*Allocate memory for communication */
    size_own = 0;
    for (pe = 0; pe < pct.num_owned_pe; pe++)
        size_own += pct.num_owned_ions_per_pe[pe];
    size_own *= ct.num_states * ct.max_nl;

    size_nown = 0;
    for (pe = 0; pe < pct.num_owners; pe++)
        size_nown += pct.num_nonowned_ions_per_pe[pe];
    size_nown *= ct.num_states * ct.max_nl;


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
                       pct.num_owned_ions_per_pe, req_recv);


    // Set sint array
    sint = kptr->newsint_local;

    for (i = 0; i < pct.num_nonloc_ions * ct.num_states * ct.max_nl; i++)
    {
        sint[i] = ZERO_t;
    }

    /*Loop over ions and calculate local projection between beta functions and wave functions */
    betaxpsi_calculate (kptr, sint, kptr->orbital_storage);


    /*Pack data for sending */
    betaxpsi_pack (sint, send_buff, pct.num_owners,
                    pct.num_nonowned_ions_per_pe, pct.list_ions_per_owner);

    /*Send <beta|psi> contributions  to the owning PE */
    betaxpsi_send (send_buff, pct.num_owners, pct.owners_list,
                    pct.num_nonowned_ions_per_pe, req_send);

    /*Wait until all data is received */
    if(pct.num_owned_pe)
        MPI_Waitall (pct.num_owned_pe, req_recv, MPI_STATUSES_IGNORE);
    if(pct.num_owners)
        MPI_Waitall (pct.num_owners, req_send, MPI_STATUSES_IGNORE);


    /*Unpack received data and sum contributions from all pes for owned ions */
    betaxpsi_sum_onwed (recv_buff, sint, pct.num_owned_pe,
                         pct.num_owned_ions_per_pe, pct.list_owned_ions_per_pe);

    /*In the second stage, owning cores will send summed data to non-owners */
    send_buff = own_buff;
    recv_buff = nown_buff;

    req_send = req_own;
    req_recv = req_nown;
    
    
    /*Receive summed data for non-owned ions from owners */
    betaxpsi_receive (recv_buff, pct.num_owners, pct.owners_list,
                       pct.num_nonowned_ions_per_pe, req_recv);

    /*Pack summed data for owned ions to send to non-owners */
    betaxpsi_pack (sint, send_buff, pct.num_owned_pe,
                    pct.num_owned_ions_per_pe, pct.list_owned_ions_per_pe);

    /*Send packed data for owned ions to non-owners */
    betaxpsi_send (send_buff, pct.num_owned_pe, pct.owned_pe_list,
                    pct.num_owned_ions_per_pe, req_send);

    /*Wait until all data is received */
    if(pct.num_owned_pe)
        MPI_Waitall (pct.num_owned_pe, req_send, MPI_STATUSES_IGNORE);
    if(pct.num_owners)
        MPI_Waitall (pct.num_owners, req_recv, MPI_STATUSES_IGNORE);

    /*Finaly, write received data about non-owned ions into sint array */
    betaxpsi_write_non_owned (sint, recv_buff, pct.num_owners,
                               pct.num_nonowned_ions_per_pe, pct.list_ions_per_owner);


    if (pct.num_owned_pe)
        delete [] req_own;
    if (pct.num_owners)
        delete [] req_nown;
    if (size_nown)
        delete [] nown_buff;
    if (size_own)
        delete [] own_buff;

}



template <typename KpointType>
void betaxpsi_calculate (Kpoint<KpointType> *kptr, KpointType * sint_ptr, KpointType *psi)
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

    KpointType *nlarray = new KpointType[pct.num_tot_proj * ct.num_states];
    
    RmgGemm (transa, transn, pct.num_tot_proj, kptr->nstates, pbasis, alpha, 
            kptr->nl_weight, pbasis, psi, pbasis, 
            rzero, nlarray, pct.num_tot_proj, kptr->nl_weight_gpu, NULLptr, NULLptr, false, false, false, true);

    for (int nion = 0; nion < pct.num_nonloc_ions; nion++)
    {
        for(int istate = 0; istate < ct.num_states; istate++)
        {
            for(int ip = 0; ip < ct.max_nl; ip++)
            {
                int proj_index = nion * ct.max_nl + ip;
                int idx1 = nion * ct.num_states * ct.max_nl + istate * ct.max_nl + ip;
                //idx2 = proj_index * ct.num_states + istate;
                int idx2 = istate * pct.num_tot_proj + proj_index;
                sint_ptr[idx1] = nlarray[idx2];
            }
        }
    }

    delete [] nlarray;
}



/*This receives data from other PEs for ions owned by current PE*/
template <typename KpointType>
void betaxpsi_receive (KpointType * recv_buff, int num_pes,
        int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
        MPI_Request * req_recv)
{
    KpointType *tpr;
    int tag, pe, source, size;

    tpr = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        source = pe_list[pe];
        /*Tag is sender_rank * NPES + receiver_rank */
        tag = 111;
        size = num_ions_per_pe[pe] * ct.num_states * ct.max_nl;
        int transfer_size = size;
        if(!ct.is_gamma) transfer_size *= 2;

        MPI_Irecv (tpr, transfer_size, MPI_DOUBLE, source, tag, pct.grid_comm, &req_recv[pe]);

        tpr += size;
    }

}

template <typename KpointType>
void betaxpsi_send (KpointType * send_buff, int num_pes,
        int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
        MPI_Request * req_send)
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
        size = num_ions * ct.num_states * ct.max_nl;
        int transfer_size = size;
        if(!ct.is_gamma) transfer_size *= 2;

        MPI_Isend (tpr, transfer_size, MPI_DOUBLE, target, tag, pct.grid_comm, &req_send[pe]);

        tpr += size;
    }

}


template <typename KpointType>
void betaxpsi_pack (KpointType * sint, KpointType * fill_buff, 
        int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS])
{
    KpointType *tpr_buff, *sint_tpr;
    int size, num_ions, ion, nlion, pe;

    tpr_buff = fill_buff;

    size = ct.num_states * ct.max_nl;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        /*Loop over ions that need to be sent to pe */
        for (ion = 0; ion < num_ions; ion++)
        {
            nlion = list_ions_per_pe[pe][ion];

            sint_tpr = &sint[nlion * ct.num_states * ct.max_nl];

            /*Pack data into send array */
            for(int idx=0;idx < size;idx++) tpr_buff[idx] = sint_tpr[idx];
            tpr_buff += size;

        }

    }

}


template <typename KpointType>
void betaxpsi_sum_onwed (KpointType * recv_buff, KpointType * sint,
        int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS])
{
    KpointType *tpr1, *tpr2;
    int size, num_ions, ion_index, pe, ion;

    size = ct.num_states * ct.max_nl;
    tpr1 = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sint[ion_index * ct.num_states * ct.max_nl];
            for(int idx=0;idx < size;idx++) tpr2[idx] += tpr1[idx];
            //my_axpy (1.0, tpr1, tpr2, size);

            tpr1 += size;
        }
    }


}



template <typename KpointType>
void betaxpsi_write_non_owned (KpointType * sint, KpointType * recv_buff, int num_pes,
        int num_ions_per_pe[MAX_NONLOC_PROCS],
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS])
{
    KpointType *tpr1, *tpr2;
    int size, num_ions, ion_index, pe, ion;

    size = ct.num_states * ct.max_nl;
    tpr1 = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sint[ion_index * ct.num_states * ct.max_nl];
            for(int idx=0;idx < size;idx++) tpr2[idx] = tpr1[idx];

            tpr1 += size;
        }
    }


}



#if 0
void betaxpsi_calculate_one(STATE *st, int ion, int nion, double *sintR, int kpt, double *weiptr_base) {

    int idx, stop, alloc, ip, incx=1, ipindex, istate, ist, st_stop, P0_BASIS;
    ION *iptr;
    SPECIES *sp;
    double *nlarrayR, *nlarrayI, *psiR, *psiI, *weiptr;
    double *pR, *pI;

    istate = st->istate;

    P0_BASIS = get_P0_BASIS();
    alloc = P0_BASIS;
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;

    my_calloc (nlarrayR, 2 * alloc, double);
    nlarrayI = nlarrayR + alloc;

    iptr = &ct.ions[ion];
    sp = &ct.sp[iptr->species];

    stop = P0_BASIS;
#if 0
    st_stop = ct.num_states / ct.THREADS_PER_NODE;
    st_stop = st_stop * ct.THREADS_PER_NODE;
#endif
    st_stop = ct.num_states;

//    for(ist = istate;ist < st_stop;ist+=ct.THREADS_PER_NODE) {
    for(ist = istate;ist < st_stop;ist++) {

        psiR = st->psiR;
        if(!ct.is_gamma) {
            psiI = st->psiI;
            pR = pct.phaseptr[ion];
            pR += 2 * kpt * stop;
            pI = pR + stop;
        }

        if(ct.is_gamma) {
            /* Copy wavefunction into temporary array */
            for (idx = 0; idx < stop; idx++)
                nlarrayR[idx] = psiR[idx];
        }
        else { 
            for (idx = 0; idx < stop; idx++)
                nlarrayR[idx] = psiR[idx] * pR[idx] - psiI[idx] * pI[idx];

            for (idx = 0; idx < stop; idx++)
                nlarrayI[idx] = psiI[idx] * pR[idx] + psiR[idx] * pI[idx];
        }

        /* <Beta|psi>                                       */

        weiptr = weiptr_base;
        ipindex = st->istate * ct.max_nl;

        for (ip = 0; ip < sp->nh; ip++)
        {

            sintR[ipindex] = get_vel() * QMD_ddot (stop, nlarrayR, incx, weiptr, incx);
            if(!ct.is_gamma) {
                sintI[ipindex] = get_vel() * QMD_ddot (stop, nlarrayI, incx, weiptr, incx);
            }
            weiptr += P0_BASIS;
            ipindex++;

        }

//        st += ct.THREADS_PER_NODE;
        st++;
    }
    my_free (nlarrayR);

}
#endif
