template <typename DataType>
void ScatterGrid(rmg::grid *G, int n, DataType *A, DataType *B, int factor)
{
    int chunksize = n / factor;
    int my_pe_x, pe_y, pe_z;
    G->pe2xyz(pct.gridpe, &my_pe_x, &pe_y, &pe_z);
    int pe_offset = my_pe_x % factor;
    CopyAndConvert(chunksize, &A[pe_offset*chunksize], &B[pe_offset*chunksize]);
    if(!ct.coalesce_states || (factor == 1)) return;

    BaseThread *T = BaseThread::getBaseThread(0);
    int tid = T->get_thread_tid();
    SCF_THREAD_CONTROL *s = (SCF_THREAD_CONTROL *)T->get_pptr(tid);
    int base_istate = 0;

    //DataType *rbuf = new CalcType[n];
    DataType *sbuf = (CalcType *)bufptr_s[tid];
    DataType *rbuf = (CalcType *)bufptr_r[tid];

    mpi_queue_item_t qitems_r[MAX_CFACTOR];
    mpi_queue_item_t qitems_s[MAX_CFACTOR];

    // States are coalesced so we have to get the remote parts of istate
    for(int i=0;i < factor;i++)
    {
        // Queue receives
        int remote_istate = base_istate + i;
        if(i != pe_offset)
        {
            qitems_r[i].comm = T->get_unique_coalesced_comm(remote_istate);
            qitems_r[i].mpi_tag = (remote_istate<<5);
            qitems_r[i].target = Rmg_T->target_node[i - pe_offset + MAX_CFACTOR][1][1];
            qitems_r[i].type = RMG_MPI_IRECV;
//            qitems_r[i].buf = (void *)&B[remote_istate*chunksize];
            qitems_r[i].buf = (void *)&rbuf[i*chunksize];
            qitems_r[i].buflen = sizeof(DataType)*chunksize;

            MPI_Irecv(qitems_r[i].buf, qitems_r[i].buflen, MPI_BYTE, qitems_r[i].target, qitems_r[i].mpi_tag, qitems_r[i].comm, &qitems_r[i].req);

        }
    }

    // Next we send the parts of states that other MPI procs require
    for(int i=0;i < factor;i++)
    {
        // Queue sends
        if(i != pe_offset)
        {
            qitems_s[i].comm = T->get_unique_coalesced_comm(pe_offset);
            qitems_s[i].mpi_tag = (pe_offset<<5);
            qitems_s[i].target = Rmg_T->target_node[i - pe_offset + MAX_CFACTOR][1][1];
            qitems_s[i].type = RMG_MPI_ISEND;
//            qitems_s[i].buf = (void *)&A[i*chunksize];
memcpy(&sbuf[i*chunksize], &A[i*chunksize], chunksize * sizeof(DataType));
            qitems_s[i].buf = (void *)&sbuf[i*chunksize];
            qitems_s[i].buflen = sizeof(DataType)*chunksize;

            MPI_Isend(qitems_s[i].buf, qitems_s[i].buflen, MPI_BYTE, qitems_s[i].target, qitems_s[i].mpi_tag, qitems_s[i].comm, &qitems_s[i].req);

        }
    }

    for(int i=0;i < factor;i++)
    {
        if(i != pe_offset)
        {
            MPI_Wait(&qitems_s[i].req, MPI_STATUS_IGNORE);
        }
        int remote_istate = base_istate + i;
        if(pe_offset != remote_istate)
        {
            MPI_Wait(&qitems_r[i].req, MPI_STATUS_IGNORE);
            CopyAndConvert(chunksize, &rbuf[i*chunksize], &B[remote_istate*chunksize]);
        }
    }

}
