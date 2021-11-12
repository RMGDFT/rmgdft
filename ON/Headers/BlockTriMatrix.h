
#ifndef BlockTriMatrix_H
#define BlockTriMatrix_H 1


#include <string>
#include <vector>
#include <complex>
#include "LocalObject.h"

/*
| A00, A01, 0, 0, 0, ...         0, A(0,nb-1)    |
| A10, A11, A12,0, .......          0            |
| 0    A21, A22, A23, 0,...         0            |
| ...                                            |
| ...                                            |
|A(nb-1,0),0 ...   0, A(nb-1,nb-2), A(nb-1,nb-1) |
*/

template <typename T> class BlockTriMatrix{

private:
    int num_blocks, mtot, ntot;
    std::vector<int> m_block_dim, n_block_dim;
    std::vector<int> n_block_idx0, m_block_idx0;
    bool dn_block;

public:
    size_t size_diag, size_up_offdiag, size_dn_offdiag;
    std::vector<T *> diag_blocks; 
    std::vector<T *> up_offdiag_blocks; 
    std::vector<T *> dn_offdiag_blocks; 
    BlockTriMatrix(int num_block_in, int m_tot_in, int n_tot_in, 
            std::vector<int> &m_block_dim_in, std::vector<int> &n_block_dim_in, bool dn_block);
    ~BlockTriMatrix(void);
    T *storage;

    void Local2BlockTri(T *mat_local,   LocalObject<T> &A, LocalObject<T> &B);
    void BlockTri2GlobDist(T *mat_dist, int *desca);


};

#endif


