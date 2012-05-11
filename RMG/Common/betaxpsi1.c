/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

#if HYBRID_MODEL
#include "hybrid.h"
#endif

static void betaxpsi1_calculate (REAL * sintR_ptr, REAL * sintI_ptr, STATE * states, int kpt);

static void betaxpsi1_receive (REAL * recv_buff, REAL * recv_buffI, int num_pes,
                               int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
                               MPI_Request * req_recv, MPI_Request * req_recvI);

static void betaxpsi1_send (REAL * send_buff, REAL * send_buffI, int num_pes,
                            int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
                            MPI_Request * req_send, MPI_Request * req_sendI);

static void betaxpsi1_pack (REAL * sintR, REAL * sintI, REAL * fill_buff, REAL * fill_buffI,
                            int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
                            int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS]);

static void betaxpsi1_sum_onwed (REAL * recv_buff, REAL * recv_buffI, REAL * sintR, REAL * sintI,
                                 int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
                                 int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS]);

static void betaxpsi1_write_non_owned (REAL * sintR, REAL * sintI, REAL * recv_buff,
                                       REAL * recv_buffI, int num_pes,
                                       int num_ions_per_pe[MAX_NONLOC_PROCS],
                                       int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS]);




void betaxpsi1 (STATE * states, int kpt)
{
    REAL *own_buff = NULL, *own_buffI = NULL, *nown_buff = NULL, *nown_buffI = NULL, *sintR = NULL, *sintI = NULL;
    REAL *send_buff, *send_buffI, *recv_buff, *recv_buffI;
    MPI_Request *req_send, *req_sendI, *req_recv, *req_recvI;
    MPI_Request *req_own = NULL, *req_ownI = NULL, *req_nown = NULL, *req_nownI = NULL;
    int size_own, size_nown, pe, koffset, i, nlion, idx;

    /*Allocate memory for communication */
    size_own = 0;
    for (pe = 0; pe < pct.num_owned_pe; pe++)
        size_own += pct.num_owned_ions_per_pe[pe];
    size_own *= ct.num_states * ct.max_nl;

    size_nown = 0;
    for (pe = 0; pe < pct.num_owners; pe++)
        size_nown += pct.num_nonowned_ions_per_pe[pe];
    size_nown *= ct.num_states * ct.max_nl;


    if (size_own)
        my_calloc (own_buff, size_own, REAL);
    if (size_nown)
        my_calloc (nown_buff, size_nown, REAL);
#if !GAMMA_PT
    if (size_own)
        my_calloc (own_buffI, size_own, REAL);
    if (size_nown)
        my_calloc (nown_buffI, size_nown, REAL);
#endif

    if (pct.num_owned_pe)
        my_malloc (req_own, pct.num_owned_pe, MPI_Request);
    if (pct.num_owners)
        my_malloc (req_nown, pct.num_owners, MPI_Request);
#if !GAMMA_PT
    if (pct.num_owned_pe)
        my_malloc (req_ownI, pct.num_owned_pe, MPI_Request);
    if (pct.num_owners)
        my_malloc (req_nownI, pct.num_owners, MPI_Request);
#endif

    /*First owning cores will receive data from non-owners */
    send_buff = nown_buff;
    recv_buff = own_buff;
    send_buffI = nown_buffI;
    recv_buffI = own_buffI;

    req_send = req_nown;
    req_recv = req_own;
    req_sendI = req_nownI;
    req_recvI = req_ownI;

    /*First post non-blocking receives for data about owned ions from cores who do not own the ions */
    betaxpsi1_receive (recv_buff, recv_buffI, pct.num_owned_pe, pct.owned_pe_list,
                       pct.num_owned_ions_per_pe, req_recv, req_recvI);


    /*Offset due to a kpoint */
    koffset = kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl;
    sintR = &pct.newsintR_local[koffset];
#if !GAMMA_PT
    sintI = &pct.newsintI_local[koffset];
#endif

    for (i = 0; i < pct.num_nonloc_ions * ct.num_states * ct.max_nl; i++)
    {
        sintR[i] = 0.0;
#if !GAMMA_PT
        sintI[i] = 0.0;
#endif
    }

    /*Loop over ions and calculate local projection between beta functions and wave functions */
    betaxpsi1_calculate (sintR, sintI, states, kpt);


    /*Pack data for sending */
    betaxpsi1_pack (sintR, sintI, send_buff, send_buffI, pct.num_owners,
                    pct.num_nonowned_ions_per_pe, pct.list_ions_per_owner);

    /*Send <beta|psi> contributions  to the owning PE */
    betaxpsi1_send (send_buff, send_buffI, pct.num_owners, pct.owners_list,
                    pct.num_nonowned_ions_per_pe, req_send, req_sendI);

    /*Wait until all data is received */
    MPI_Waitall (pct.num_owned_pe, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall (pct.num_owners, req_send, MPI_STATUSES_IGNORE);

#if !GAMMA_PT
    MPI_Waitall (pct.num_owned_pe, req_recvI, MPI_STATUSES_IGNORE);
    MPI_Waitall (pct.num_owners, req_sendI, MPI_STATUSES_IGNORE);
#endif

    /*Unpack received data and sum contributions from all pes for owned ions */
    betaxpsi1_sum_onwed (recv_buff, recv_buffI, sintR, sintI, pct.num_owned_pe,
                         pct.num_owned_ions_per_pe, pct.list_owned_ions_per_pe);

    /*In the second stage, owning cores will send summed data to non-owners */
    send_buff = own_buff;
    recv_buff = nown_buff;
    send_buffI = own_buffI;
    recv_buffI = nown_buffI;

    req_send = req_own;
    req_recv = req_nown;
    req_sendI = req_ownI;
    req_recvI = req_nownI;
    
    
    /*Receive summed data for non-owned ions from owners */
    betaxpsi1_receive (recv_buff, recv_buffI, pct.num_owners, pct.owners_list,
                       pct.num_nonowned_ions_per_pe, req_recv, req_recvI);

    /*Pack summed data for owned ions to send to non-owners */
    betaxpsi1_pack (sintR, sintI, send_buff, send_buffI, pct.num_owned_pe,
                    pct.num_owned_ions_per_pe, pct.list_owned_ions_per_pe);

    /*Send packed data for owned ions to non-owners */
    betaxpsi1_send (send_buff, send_buffI, pct.num_owned_pe, pct.owned_pe_list,
                    pct.num_owned_ions_per_pe, req_send, req_sendI);

    /*Wait until all data is received */
    MPI_Waitall (pct.num_owned_pe, req_send, MPI_STATUSES_IGNORE);
    MPI_Waitall (pct.num_owners, req_recv, MPI_STATUSES_IGNORE);

#if !GAMMA_PT
    MPI_Waitall (pct.num_owned_pe, req_sendI, MPI_STATUSES_IGNORE);
    MPI_Waitall (pct.num_owners, req_recvI, MPI_STATUSES_IGNORE);
#endif

    /*Finaly, write received data about non-owned ions into sintR array */
    betaxpsi1_write_non_owned (sintR, sintI, recv_buff, recv_buffI, pct.num_owners,
                               pct.num_nonowned_ions_per_pe, pct.list_ions_per_owner);


    if (size_own)
        my_free (own_buff);
    if (size_nown)
        my_free (nown_buff);
    if (pct.num_owned_pe)
        my_free (req_own);
    if (pct.num_owners)
        my_free (req_nown);

#if !GAMMA_PT
    if (size_own)
        my_free (own_buffI);
    if (size_nown)
        my_free (nown_buffI);
    if (pct.num_owned_pe)
        my_free (req_ownI);
    if (pct.num_owners)
        my_free (req_nownI);
#endif


}



static void betaxpsi1_calculate (REAL * sintR_ptr, REAL * sintI_ptr, STATE * states, int kpt)
{
    int alloc, nion, ion, *pidx, istate, idx, ipindex, stop, ip, incx = 1, start_state, istop, ist;
    REAL *nlarrayR, *nlarrayI, *sintR, *sintI, *pR, *pI;
    REAL *weiptr, *psiR, *psiI;
    ION *iptr;
    SPECIES *sp;
    STATE *st;


    alloc = P0_BASIS;
    if (alloc < ct.max_nlpoints)
        alloc = ct.max_nlpoints;

#if HYBRID_MODEL
    my_calloc (nlarrayR, 2 * alloc * THREADS_PER_NODE, REAL);
#else
    my_calloc (nlarrayR, 2 * alloc, REAL);
#endif
    nlarrayI = nlarrayR + alloc;



    /* Loop over ions on this processor */
    for (nion = 0; nion < pct.num_nonloc_ions; nion++)
    {

        /*Actual index of the ion under consideration */
        ion = pct.nonloc_ions_list[nion];


        iptr = &ct.ions[ion];
        sp = &ct.sp[iptr->species];
        stop = pct.idxptrlen[ion];


        if (stop)
        {

            pidx = pct.nlindex[ion];
#if !GAMMA_PT
            pR = pct.phaseptr[ion];
            pR += 2 * kpt * stop;
            pI = pR + stop;
#endif

            sintR = &sintR_ptr[nion * ct.num_states * ct.max_nl];
#if !GAMMA_PT
            sintI = &sintI_ptr[nion * ct.num_states * ct.max_nl];
#endif

// Parallelized over threads here.
            start_state = 0;
#if HYBRID_MODEL
            enter_threaded_region();
            scf_barrier_init(THREADS_PER_NODE);
            istop = ct.num_states / THREADS_PER_NODE;
            istop = istop * THREADS_PER_NODE;

            for(ist = 0;ist < THREADS_PER_NODE;ist++) {
                  thread_control[ist].job = HYBRID_BETAX_PSI1_CALCULATE;
                  thread_control[ist].sp = &states[ist];
                  thread_control[ist].ion = ion;
                  thread_control[ist].nion = nion;
                  thread_control[ist].sintR = sintR;
                  thread_control[ist].sintI = sintI;
            }

            // Thread tasks are set up so wake them
            wake_threads(THREADS_PER_NODE);

            // Then wait for them to finish this task
            wait_for_threads(THREADS_PER_NODE);

            scf_barrier_destroy();
            leave_threaded_region();
            start_state = istop;

#endif

            // Handle the remainder of the states in serial fashion. If
            // not running in hybrid mode then start_state is 0.
            for (istate = start_state; istate < ct.num_states; istate++)
            {

                st = &states[istate];
                psiR = st->psiR;
#if !GAMMA_PT
                psiI = st->psiI;
#endif


#if GAMMA_PT
                /* Copy wavefunction into temporary array */
                for (idx = 0; idx < stop; idx++)
                    nlarrayR[idx] = psiR[pidx[idx]];
#else
                for (idx = 0; idx < stop; idx++)
                    nlarrayR[idx] = psiR[pidx[idx]] * pR[idx] - psiI[pidx[idx]] * pI[idx];

                for (idx = 0; idx < stop; idx++)
                    nlarrayI[idx] = psiI[pidx[idx]] * pR[idx] + psiR[pidx[idx]] * pI[idx];
#endif

                /* <Beta|psi>                                       */

                weiptr = pct.weight[ion];
                ipindex = istate * ct.max_nl;

                for (ip = 0; ip < sp->nh; ip++)
                {

                    sintR[ipindex] = ct.vel * QMD_sdot (stop, nlarrayR, incx, weiptr, incx);
#if !GAMMA_PT
                    sintI[ipindex] = ct.vel * QMD_sdot (stop, nlarrayI, incx, weiptr, incx);
#endif

                    weiptr += pct.idxptrlen[ion];
                    ipindex++;

                }

            }

        }


    }
    my_free (nlarrayR);
}

void betaxpsi1_calculate_one(STATE *st, int ion, int nion, REAL *sintR, REAL *sintI, int kpt) {

  int idx, stop, alloc, ip, incx=1, ipindex, istate, *pidx, ist, st_stop;
  ION *iptr;
  SPECIES *sp;
  REAL *nlarrayR, *nlarrayI, *psiR, *psiI, *weiptr;
    REAL *pR, *pI;

  istate = st->istate;

  alloc = P0_BASIS;
  if (alloc < ct.max_nlpoints)
      alloc = ct.max_nlpoints;

  my_malloc (nlarrayR, 2 * alloc, REAL);
  nlarrayI = nlarrayR + alloc;

  iptr = &ct.ions[ion];
  sp = &ct.sp[iptr->species];
  stop = pct.idxptrlen[ion];
  pidx = pct.nlindex[ion];

  st_stop = ct.num_states / THREADS_PER_NODE;
  st_stop = st_stop * THREADS_PER_NODE;

  for(ist = istate;ist < st_stop;ist+=THREADS_PER_NODE) {
      
      psiR = st->psiR;
#if !GAMMA_PT
      psiI = st->psiI;
      pR = pct.phaseptr[ion];
      pR += 2 * kpt * stop;
      pI = pR + stop;
#endif

#if GAMMA_PT
      /* Copy wavefunction into temporary array */
      for (idx = 0; idx < stop; idx++)
          nlarrayR[idx] = psiR[pidx[idx]];
#else
      for (idx = 0; idx < stop; idx++)
          nlarrayR[idx] = psiR[pidx[idx]] * pR[idx] - psiI[pidx[idx]] * pI[idx];

      for (idx = 0; idx < stop; idx++)
          nlarrayI[idx] = psiI[pidx[idx]] * pR[idx] + psiR[pidx[idx]] * pI[idx];
#endif

      /* <Beta|psi>                                       */

      weiptr = pct.weight[ion];
      ipindex = st->istate * ct.max_nl;

      for (ip = 0; ip < sp->nh; ip++)
      {

          sintR[ipindex] = ct.vel * QMD_sdot (stop, nlarrayR, incx, weiptr, incx);
#if !GAMMA_PT
          sintI[ipindex] = ct.vel * QMD_sdot (stop, nlarrayI, incx, weiptr, incx);
#endif

          weiptr += pct.idxptrlen[ion];
          ipindex++;

      }

      st += THREADS_PER_NODE;
  }
  my_free (nlarrayR);

}


/*This receives data from other PEs for ions owned by current PE*/
static void betaxpsi1_receive (REAL * recv_buff, REAL * recv_buffI, int num_pes,
        int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
        MPI_Request * req_recv, MPI_Request * req_recvI)
{
    REAL *tpr, *tprI;
    int tag, pe, source, size;

    tpr = recv_buff;
    tprI = recv_buffI;

    for (pe = 0; pe < num_pes; pe++)
    {
        source = pe_list[pe];
        /*Tag is sender_rank * NPES + receiver_rank */
        tag = 111;
        size = num_ions_per_pe[pe] * ct.num_states * ct.max_nl;

        MPI_Irecv (tpr, size, MPI_DOUBLE, source, tag, pct.grid_comm, &req_recv[pe]);
	Dprintf("Posting nonblock receive from PE %d tag is %d size is %d", source, tag, size);
#if !GAMMA_PT
        tag *= 2;
        MPI_Irecv (tprI, size, MPI_DOUBLE, source, tag, pct.grid_comm, &req_recvI[pe]);
	Dprintf("Posting nonblock receive for imaginary part from PE %d tag is %d size is %d", source, tag, size);
#endif

        tpr += size;
        tprI += size;
    }

}

static void betaxpsi1_send (REAL * send_buff, REAL * send_buffI, int num_pes,
        int pe_list[MAX_NONLOC_PROCS], int num_ions_per_pe[MAX_NONLOC_PROCS],
        MPI_Request * req_send, MPI_Request * req_sendI)
{
    REAL *tpr, *tprI;
    int target, num_ions, size, tag, pe;

    tpr = send_buff;
    tprI = send_buffI;

    for (pe = 0; pe < num_pes; pe++)
    {

        target = pe_list[pe];
        num_ions = num_ions_per_pe[pe];
        size = num_ions * ct.num_states * ct.max_nl;

        /*Tag is sender_rank * NPES + receiver_rank */
        tag = 111;

        MPI_Isend (tpr, size, MPI_DOUBLE, target, tag, pct.grid_comm, &req_send[pe]);
	Dprintf("Sending data to PE %d, tag is %d and size is %d", target, tag, size);

#if !GAMMA_PT
        tag *= 2;
        MPI_Isend (tprI, size, MPI_DOUBLE, target, tag, pct.grid_comm, &req_sendI[pe]);
	Dprintf("Sending imaginary part of data to PE %d, tag is %d and size is %d", target, tag, size);
#endif

        tpr += size;
        tprI += size;
    }

}


static void betaxpsi1_pack (REAL * sintR, REAL * sintI, REAL * fill_buff, REAL * fill_buffI,
        int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS])
{
    REAL *tpr_buff, *tpr_buffI, *sintR_tpr, *sintI_tpr;
    int size, num_ions, ion, nlion, pe;

    tpr_buff = fill_buff;
    tpr_buffI = fill_buffI;

    size = ct.num_states * ct.max_nl;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        /*Loop over ions that need to be sent to pe */
        for (ion = 0; ion < num_ions; ion++)
        {
            nlion = list_ions_per_pe[pe][ion];

            sintR_tpr = &sintR[nlion * ct.num_states * ct.max_nl];

            /*Pack data into send array */
            my_copy (sintR_tpr, tpr_buff, size);
            tpr_buff += size;

#if !GAMMA_PT
            sintI_tpr = &sintI[nlion * ct.num_states * ct.max_nl];
            my_copy (sintI_tpr, tpr_buffI, &size);
            tpr_buffI += size;
#endif
        }

    }

}


static void betaxpsi1_sum_onwed (REAL * recv_buff, REAL * recv_buffI, REAL * sintR, REAL * sintI,
        int num_pes, int num_ions_per_pe[MAX_NONLOC_PROCS],
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS])
{
    REAL *tpr1, *tpr1I, *tpr2, *tpr2I;
    int size, num_ions, ion_index, pe, ion;


    size = ct.num_states * ct.max_nl;
    tpr1 = recv_buff;
    tpr1I = recv_buffI;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sintR[ion_index * ct.num_states * ct.max_nl];
            my_axpy (1.0, tpr1, tpr2, size);

#if !GAMMA_PT
            tpr2I = &sintI[ion_index * ct.num_states * ct.max_nl];
            my_axpy (1.0, tpr1I, tpr2I, size);
#endif
            tpr1 += size;
            tpr1I += size;
        }
    }


}



static void betaxpsi1_write_non_owned (REAL * sintR, REAL * sintI, REAL * recv_buff,
        REAL * recv_buffI, int num_pes,
        int num_ions_per_pe[MAX_NONLOC_PROCS],
        int list_ions_per_pe[MAX_NONLOC_PROCS][MAX_NONLOC_IONS])
{
    REAL *tpr1, *tpr1I, *tpr2, *tpr2I;
    int size, num_ions, ion_index, pe, ion;


    size = ct.num_states * ct.max_nl;
    tpr1 = recv_buff;
    tpr1I = recv_buffI;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sintR[ion_index * ct.num_states * ct.max_nl];
            my_copy (tpr1, tpr2, size);

#if !GAMMA_PT
            tpr2I = &sintI[ion_index * ct.num_states * ct.max_nl];
            my_copy (tpr1I, tpr2I, size);
#endif
            tpr1 += size;
            tpr1I += size;
        }
    }


}
