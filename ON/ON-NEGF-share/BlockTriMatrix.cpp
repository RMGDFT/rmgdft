#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "main.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "blas.h"
#include "RmgSumAll.h"

#include "prototypes_on.h"
#include "init_var.h"
#include "PulayMixing.h"
#include "LocalObject.h"
#include "Kbpsi.h"
#include "LdaU_on.h"
#include "GpuAlloc.h"
#include "Exx_on.h"
#include "BlockTriMatrix.h"

#include "BlockTriMatrix.h"
#include "LocalObject.h"
#include "blas.h"
#include "RmgException.h"
#include "transition.h"
#include "prototypes_on.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "GlobalSums.h"
#include "blacs.h"
template BlockTriMatrix<double>::BlockTriMatrix(
        int num_block_in, int mtot_in, int ntot_in, 
        std::vector<int> &m_block_dim_in, std::vector<int> &n_block_dim_in, bool dn_block);
template BlockTriMatrix<double>::~BlockTriMatrix(void);
template <class T> BlockTriMatrix<T>::~BlockTriMatrix()
{
    delete [] storage;
}

template <class T> BlockTriMatrix<T>::BlockTriMatrix(
        int num_block_in, int mtot_in, int ntot_in, 
        std::vector<int> &m_block_dim_in, std::vector<int> &n_block_dim_in, bool dn_block_in) : num_blocks(num_block_in),
    mtot(mtot_in), ntot(ntot_in), m_block_dim(m_block_dim_in), n_block_dim(n_block_dim_in), dn_block(dn_block_in)
{

    diag_blocks.resize(num_blocks);
    up_offdiag_blocks.resize(num_blocks);
    dn_offdiag_blocks.resize(num_blocks);
    n_block_idx0.resize(num_blocks);
    m_block_idx0.resize(num_blocks);
    size_diag=0;
    size_up_offdiag=0;
    size_dn_offdiag=0;

    n_block_idx0[0] = 0;
    m_block_idx0[0] = 0;
    for(int ib = 1; ib < num_blocks; ib++)
    {
        n_block_idx0[ib] += n_block_dim[ib-1];
        m_block_idx0[ib] += m_block_dim[ib-1];
    }

    for(int i = 0; i < num_blocks; i++)
    {
        size_diag += n_block_dim[i] * m_block_dim[i];
        int j = (i+1)%num_blocks;
        size_up_offdiag += m_block_dim[i] * n_block_dim[j];
        size_dn_offdiag += m_block_dim[j] * n_block_dim[i];
    }

    if(!dn_block) 
    {
        size_dn_offdiag = 0;
    }
    storage = new T[size_diag+size_up_offdiag + size_dn_offdiag];

    T *diag_ptr = storage;
    T *up_offdiag_ptr = storage+ size_diag;
    T *dn_offdiag_ptr = storage+ size_diag + size_up_offdiag;

    for(int i = 0; i < num_blocks; i++)
    {
        diag_blocks[i] = diag_ptr;
        up_offdiag_blocks[i] = up_offdiag_ptr;
        if(dn_block) {
            dn_offdiag_blocks[i] = dn_offdiag_ptr;
        }
        if (i == num_blocks -1) break;
        diag_ptr = diag_ptr +  m_block_dim[i] * n_block_dim[i];
        up_offdiag_ptr = up_offdiag_ptr +  m_block_dim[i] * n_block_dim[i+1];
        if(dn_block) {
            dn_offdiag_ptr = dn_offdiag_ptr +  m_block_dim[i+1] * n_block_dim[i];
        }
    }

}


template void BlockTriMatrix<double>::Local2BlockTri(double *mat_local,  LocalObject<double> &A, LocalObject<double> &B);

template <class T> void BlockTriMatrix<T>::Local2BlockTri(T *mat_local,   LocalObject<T> &A, LocalObject<T> &B)
{
    if(A.density != B.density)
        throw RmgFatalException() << "density is different "<< " at file " << __FILE__ << "\n";

    RmgTimer *RT= new RmgTimer("BlockTr");
    for(int ib = 0; ib < num_blocks; ib++) {
        int jb = (ib + 1)%num_blocks;
        for(int st1 = 0; st1 < m_block_dim[ib]; st1++) {
            int st1_glob = m_block_idx0[ib] + st1;
            int st1_local = A.index_global_to_proj[st1_glob];
            if (st1_local < 0 ) continue;

            for(int st2 = 0; st2 < n_block_dim[ib]; st2++) {
                int st2_glob = n_block_idx0[ib] + st2;
                int st2_local = B.index_global_to_proj[st2_glob];
                if (st2_local < 0 ) continue;

                diag_blocks[ib][st2 * m_block_dim[ib] + st1] += mat_local[st2_local * A.num_thispe + st1_local];
            }

            for(int st2 = 0; st2 < n_block_dim[jb]; st2++) {
                int st2_glob = n_block_idx0[jb] + st2;
                int st2_local = B.index_global_to_proj[st2_glob];
                if (st2_local < 0 ) continue;

                up_offdiag_blocks[ib][st2 * m_block_dim[ib] + st1] += mat_local[st2_local * A.num_thispe + st1_local];
            }
        }

        if(dn_block)
        {
            for(int st1 = 0; st1 < m_block_dim[jb]; st1++) {
                int st1_glob = m_block_idx0[jb] + st1;
                int st1_local = A.index_global_to_proj[st1_glob];
                if (st1_local < 0 ) continue;
                for(int st2 = 0; st2 < n_block_dim[ib]; st2++) {
                    int st2_glob = n_block_idx0[ib] + st2;
                    int st2_local = B.index_global_to_proj[st2_glob];
                    if (st2_local < 0 ) continue;
                    dn_offdiag_blocks[ib][st2 * m_block_dim[jb] + st1] += mat_local[st2_local * A.num_thispe + st1_local];
                }
            }
        }
    }
    delete RT;


    MPI_Barrier(MPI_COMM_WORLD);
    RT= new RmgTimer("BlockTrReduc");

    int tot_size = size_diag + size_up_offdiag + size_dn_offdiag;
    //printf("\n tot_sie %d", tot_size);
    //MPI_Allreduce(MPI_IN_PLACE, (double *)storage, tot_size, MPI_DOUBLE, MPI_SUM, A.comm);
    //   GlobalSums(storage_up_offdiag, size_up_offdiag, A.comm);
    MPI_Barrier(MPI_COMM_WORLD);
    delete RT;

}

template void BlockTriMatrix<double>::BlockTri2GlobDist(double *mat_dist, int *desca); 

template <class T> void BlockTriMatrix<T>::BlockTri2GlobDist(T *a_dist, int *desca) 
{

    int mycol, myrow, nprow, npcol;
    int ictxt=desca[1], m = desca[2], n = desca[3], mb=desca[4], nb=desca[5], mxllda = desca[8];

    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    int izero = 0;
    int mxli = numroc(&m, &mb, &myrow, &izero, &nprow);
    int mxlocc = numroc(&n, &nb, &mycol, &izero, &npcol);

    for(int j =0; j < mxli; j++)
    {
        for(int k=0; k < mxlocc; k++)
        {
            a_dist[k * mxllda + j] = 0.0;
        }
    }


    for(int j =0; j < mxli; j++)
    {
        for(int k=0; k < mxlocc; k++)
        {

            /*  distributed local index (j,k) is (jj, kk) in global matrix
             */

            int jj = (j/mb ) * nprow * mb + myrow * mb + j - j/mb * mb;
            int kk = (k/nb ) * npcol * nb + mycol * nb + k - k/nb * nb;

            int jblock=0, kblock=0;
            for(int ib = 0; ib < num_blocks; ib++) {
                if(m_block_idx0[ib] >= jj) {
                    jblock = ib;
                    break;
                } 
            }
            for(int ib = 0; ib < num_blocks; ib++) {
                if(n_block_idx0[ib] >= kk) {
                    kblock = ib;
                    break;
                } 
            }

            int jj_tri = jj - m_block_idx0[jblock];
            int kk_tri = kk - n_block_idx0[kblock];

            if(jblock == kblock) {
                a_dist[k * mxllda + j] = diag_blocks[jblock][kk_tri * m_block_dim[jblock] + jj_tri];
            }
            else if( (jblock + 1 == kblock) || (jblock == num_blocks -1 && kblock == 0) ) {
                a_dist[k * mxllda + j] = up_offdiag_blocks[jblock][kk_tri * m_block_dim[jblock] + jj_tri];
            }
            else if( (jblock  == kblock + 1) || (kblock == num_blocks -1 && jblock == 0) ) {
                if(dn_block) {
                    a_dist[k * mxllda + j] = dn_offdiag_blocks[jblock][kk_tri * m_block_dim[kblock] + jj_tri];
                }
                else {
                    a_dist[k * mxllda + j] = up_offdiag_blocks[jblock][jj_tri * m_block_dim[kblock] + kk_tri];
                }
            }

        }
    }

}
