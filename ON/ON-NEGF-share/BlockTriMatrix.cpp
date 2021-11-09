#include "BlockTriMatrix.h"
#include "LocalObject.h"
#include "blas.h"
#include "RmgException.h"
#include "transition.h"
#include "prototypes_on.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "GlobalSums.h"
template BlockTriMatrix<double>::BlockTriMatrix(
        int num_block_in, int ntot_in, int mtot_in, 
        std::vector<int> n_block_dim_in, std::vector<int> m_block_dim_in, bool dn_block);
template BlockTriMatrix<double>::~BlockTriMatrix(void);
template <class T> BlockTriMatrix<T>::~BlockTriMatrix()
{
    delete [] storage_diag;
    delete [] storage_up_offdiag;
    if(dn_block) {
        delete [] storage_dn_offdiag;
    }
}

template <class T> BlockTriMatrix<T>::BlockTriMatrix(
        int num_block_in, int ntot_in, int mtot_in, 
        std::vector<int> n_block_dim_in, std::vector<int> m_block_dim_in, bool dn_block_in) : num_blocks(num_block_in),
    ntot(ntot_in), mtot(mtot_in), n_block_dim(n_block_dim_in), m_block_dim(m_block_dim_in), dn_block(dn_block_in)
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
        size_up_offdiag += n_block_dim[i] * m_block_dim[j];
        size_dn_offdiag += n_block_dim[j] * m_block_dim[i];
    }

    storage_diag = new T[size_diag];
    storage_up_offdiag = new T[size_up_offdiag];
    if(dn_block) {
        storage_dn_offdiag = new T[size_dn_offdiag];
    }
    else {
        size_dn_offdiag = 0;
    }

    T *diag_ptr = storage_diag;
    T *up_offdiag_ptr = storage_up_offdiag;
    T *dn_offdiag_ptr = storage_dn_offdiag;

    for(int i = 0; i < num_blocks; i++)
    {
        diag_blocks[i] = diag_ptr;
        up_offdiag_blocks[i] = up_offdiag_ptr;
        if(dn_block) {
            dn_offdiag_blocks[i] = dn_offdiag_ptr;
        }
        if (i == num_blocks -1) break;
        diag_ptr = diag_ptr +  n_block_dim[i] * m_block_dim[i];
        up_offdiag_ptr = up_offdiag_ptr +  n_block_dim[i] * m_block_dim[i+1];
        if(dn_block) {
            dn_offdiag_ptr = dn_offdiag_ptr +  n_block_dim[i+1] * m_block_dim[i];
        }
    }

}


template void BlockTriMatrix<double>::Local2BlockTri(double *mat_local,  LocalObject<double> &A, LocalObject<double> &B);

template <class T> void BlockTriMatrix<T>::Local2BlockTri(T *mat_local,   LocalObject<T> &A, LocalObject<T> &B)
{
    if(A.density != B.density)
        throw RmgFatalException() << "density is different "<< " at file " << __FILE__ << "\n";

    for(int ib = 0; ib < num_blocks; ib++) {
        int jb = ib + 1;
        for(int st1 = 0; st1 < n_block_dim[ib]; st1++) {
            int st1_glob = n_block_idx0[ib] + st1;
            int st1_local = A.index_global_to_proj[st1_glob];
            if (st1_local < 0 ) continue;

            for(int st2 = 0; st2 < m_block_dim[ib]; st2++) {
                int st2_glob = m_block_idx0[ib] + st2;
                int st2_local = B.index_global_to_proj[st2_glob];
                if (st2_local < 0 ) continue;

                diag_blocks[ib][st2 * n_block_dim[ib] + st1] += mat_local[st2_local * A.num_thispe + st1_local];
            }

            for(int st2 = 0; st2 < m_block_dim[jb]; st2++) {
                int st2_glob = m_block_idx0[jb] + st2;
                int st2_local = B.index_global_to_proj[st2_glob];
                if (st2_local < 0 ) continue;

               up_offdiag_blocks[ib][st2 * n_block_dim[ib] + st1] += mat_local[st2_local * A.num_thispe + st1_local];
            }
        }

        if(dn_block)
        {
            for(int st1 = 0; st1 < n_block_dim[jb]; st1++) {
                int st1_glob = n_block_idx0[jb] + st1;
                int st1_local = A.index_global_to_proj[st1_glob];
                if (st1_local < 0 ) continue;
                for(int st2 = 0; st2 < m_block_dim[ib]; st2++) {
                    int st2_glob = m_block_idx0[ib] + st2;
                    int st2_local = B.index_global_to_proj[st2_glob];
                    if (st2_local < 0 ) continue;
                    dn_offdiag_blocks[ib][st2 * n_block_dim[jb] + st1] += mat_local[st2_local * A.num_thispe + st1_local];
                }
            }
        }
    }


    GlobalSums(storage_diag, size_diag, A.comm);
    GlobalSums(storage_up_offdiag, size_up_offdiag, A.comm);
    if(dn_block)
        GlobalSums(storage_dn_offdiag, size_dn_offdiag, A.comm);

}
