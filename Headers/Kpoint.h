/*
 *
 * Copyright (c) 2014, Emil Briggs
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

#ifndef RMG_Kpoint_H
#define RMG_Kpoint_H 1

#include "BaseGrid.h"
#include "Lattice.h"
#include "TradeImages.h"
#include "State.h"
#include "InputKey.h"
#include "Projector.h"
#include "RmgParallelFft.h"
#include "LdaU.h"
#include <mpi.h>
#include <boost/pool/pool.hpp>

//char *rmg_pool_malloc(const size_t bytes);
// Used for kalloc
char inline *rmg_pool_malloc(const size_t bytes)
    {
    #if HIP_ENABLED
        void *ptr;
        hipHostMalloc(&ptr, bytes, hipHostMallocNumaUser);
        return reinterpret_cast<char *>(ptr);
    #else
        return reinterpret_cast<char *>(std::malloc(bytes));
    #endif
    }


struct rmg_user_allocator
{
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    static char * malloc(const size_t bytes)
    {return reinterpret_cast<char *>(rmg_pool_malloc(bytes)); }
    static void free(char * const block)
#if HIP_ENABLED
    { hipHostFree(block); }
#else
    { std::free(block); }
#endif
};


template <typename KpointType> class Kpoint {

public:

    Kpoint(KSTRUCT &kpin, int index, MPI_Comm newcomm, BaseGrid *newG, TradeImages *newT, Lattice *newL, std::unordered_map<std::string, InputKey *>& ControlMap);

    void set_pool(KpointType *pool);
    void random_init(void);
    int get_nstates(void);
    int get_index(void);
    void orthogonalize(double *storage);
    void orthogonalize(std::complex<double> *storage);
    void init_states(void);
    void write_occ(void);
    void get_nlop(int type);
    void get_ldaUop(int type);
    void get_orbitals(KpointType *orbitals);
    void get_ion_orbitals(ION *iptr, KpointType *orbitals);
    void reset_beta_arrays(void);
    void reset_orbital_arrays(void);
    void Subdiag (double *vtot_eig, double *vxc_psi, int subdiag_driver);
    void ComputeHcore (double *vtot_eig, double *vxc_psi, KpointType *Hcore, KpointType *Hcore_kin);
    void MgridSubspace (double *vtot_psi, double *vxc_psi);
    void Davidson(double *vtot, double *vxc_psi, int &notconv);
    void GetLocalizedWeight (void);
    void GetDelocalizedWeight (void);
    void GetDelocalizedOrbital (void);
    void LcaoGetPsi (void);
    void DeleteNvmeArrays(void);
    void ClearPotentialAcceleration(void);


    // Minimal kpoint structure
    KSTRUCT &kp;

    // Input file internal map
    std::unordered_map<std::string, InputKey *>& ControlMap;

    // BaseGrid class
    BaseGrid *G;

    // TradeImages object to use
    TradeImages *T;

    // Lattice object
    Lattice *L;

    // Beta function projectors
    Projector<KpointType> *BetaProjector;

    // Atomic orbital projectors
    Projector<KpointType> *OrbitalProjector;

    // ldaU object
    LdaU<KpointType> *ldaU;

    // MPI communicator to use for trade images and reduction operations
    MPI_Comm grid_comm;

    // The index of the k-point for backreferencing
    int kidx;

    // Number of orbitals
    int nstates;

    // Block of contiguous storage for the orbitals
    KpointType *orbital_storage;

    // Block of contiguous storage for the orbitals from the previous step which is needed in some cases
    KpointType *prev_orbitals;

    // The orbital structure for this k-point
    State<KpointType> *Kstates;

    // Pointer to sint arrays (Betaxpsi)
    KpointType *newsint_local;

    // Pointer to orbital sint arrays (Orbital projectors)
    KpointType *orbitalsint_local;

    // Size of the sint arrays
    //size_t sint_size;

    // Pointers to nv, ns. Each of these arrays is dimensioned (NL_BLOCK_SIZE, P0_BASIS).
    // since we apply the non-local operators in blocks for efficiency and to save memory
    KpointType *nv;
    KpointType *ns;
    int nl_first_state;  // first state in the buffer

    // Pointers to weight and Bweight
    KpointType *nl_weight;
#if HIP_ENABLED || CUDA_ENABLED
    KpointType *nl_weight_gpu;
#endif
    size_t nl_weight_size;

    //Pointer to vexx
    KpointType *vexx;

    // Orbital weights
    KpointType *orbital_weight;
    size_t orbital_weight_size;

    // Pointer to potential acceleration arrays and size
    double *dvh;
    size_t dvh_size;

    // Number of potential acceleration arrays and the skip factor
    int ndvh;
    int dvh_skip;

    // Number of points in orbital basis
    int pbasis;
    int pbasis_noncoll;


    // Highest occupied orbital
    int highest_occupied;

    // We make this static since only one kpoint runs at a time and a work memory pool
    // can be shared between kpoints.
    static std::vector<boost::pool<rmg_user_allocator> *> kalloc;


private:

    // Mean min, and max wavefunction residuals for occupied space
    double mean_occ_res;
    double min_occ_res;
    double max_occ_res;

    // Mean min, and max wavefunction residuals for unoccupied space
    double mean_unocc_res;
    double min_unocc_res;
    double max_unocc_res;

    // Index of the highest orbital included in the calculation of mean/min/max
    int max_unocc_res_index;

    int nvme_weight_fd;
    int nvme_Bweight_fd;
    int nvme_ldaU_fd;
    std::string nvme_weight_path;
    std::string nvme_Bweight_path;
    std::string nvme_ldaU_path;

};

template <typename T>
using Kpoints = std::vector< Kpoint<T> >;

#include "../../RMG/Headers/prototypes_rmg.h"
#endif
