#pragma once

#include "cartesian_transposer.hpp"
#include "custom_stencils.hpp"
#include "multigrid_preconditioner.hpp"
#include "chebyshev_rmmdiis.hpp"
#include "rayleigh_diis_helpers.hpp"
#include "defect_correction_solver.hpp"

#include <mpi.h>
#include <hip/hip_runtime.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <map>
#include <memory>
#include <vector>

template < typename Canonical64T >
struct DualPrecisionTypes;

template <>
struct DualPrecisionTypes< double > {
    using Wave32T = float;
    using Wave64T = double;
    using Real64T = double;
    using LowT    = float;
    using LowWtT  = float;
};

template <>
struct DualPrecisionTypes< std::complex< double > > {
    using Wave32T = std::complex< float >;
    using Wave64T = std::complex< double >;
    using Real64T = double;
    using LowT    = std::complex< float >;
    using LowWtT  = std::complex< float >;
};

template < typename Canonical64T >
class DynamicPrecisionGpuSolver {
public:
    using Types   = DualPrecisionTypes< Canonical64T >;
    using Wave32T = typename Types::Wave32T;
    using Wave64T = typename Types::Wave64T;
    using Real64T = typename Types::Real64T;
    using LowT    = typename Types::LowT;
    using LowWtT  = typename Types::LowWtT;

    static constexpr int DOT_BLOCKS = 256;
    static constexpr int MAX_DIIS   = 3;

    struct BlockDiisState {
        int hist_len = 0;
        std::array< void*, MAX_DIIS > raw_hist_psi = {{nullptr, nullptr, nullptr}};
        std::array< LowT*, MAX_DIIS > d_hist_r     = {{nullptr, nullptr, nullptr}};
        void** raw_hist_psi_ptrs = nullptr;
        LowT** d_hist_r_ptrs     = nullptr;
    };

    DynamicPrecisionGpuSolver(
        MPI_Comm parent_comm,
        int max_block_states,
        int Nx, int Ny, int Nz,
        SubBox3D my_box,
        int mg_num_levels,
        const std::vector< StencilTap< std::complex< double > > >& fine_taps_8th,
        const std::vector< std::vector< StencilTap< std::complex< double > > > >& coarse_taps_2nd,
        double fine_lambda_min,
        double fine_lambda_max,
        double helmholtz_shift,
        double vel = 1.0,
        bool periodic = true,
        int desired_mpi_batch = 1024,
        int pre_smooth_deg = 2,
        int post_smooth_deg = 2,
        int coarsest_solve_deg = 6)
        : max_block_states_(max_block_states),
          Nx_(Nx), Ny_(Ny), Nz_(Nz),
          vel_(vel),
          periodic_(periodic),
          fine_taps_(fine_taps_8th) {
        MPI_Comm_dup(parent_comm, &comm_);
        MPI_Comm_size(comm_, &n_ranks_);
        MPI_Comm_rank(comm_, &rank_);

        Halo3D halo{4, 4, 4};

        transposer_32_ = std::make_unique< CartesianWavefunctionTransposer< Wave32T > >(
            comm_, max_block_states_, Nx_, Ny_, Nz_, my_box, halo, desired_mpi_batch, true);
        transposer_64_ = std::make_unique< CartesianWavefunctionTransposer< Wave64T > >(
            comm_, max_block_states_, Nx_, Ny_, Nz_, my_box, halo, desired_mpi_batch, true);

        my_max_states_ = transposer_64_->my_band_states();
        Gx_ = transposer_64_->Gx();
        Gy_ = transposer_64_->Gy();
        Gz_ = transposer_64_->Gz();

        mg_ = std::make_unique< MultigridPreconditioner< LowT > >(
            my_max_states_, Nx_, Ny_, Nz_, mg_num_levels,
            fine_taps_8th, coarse_taps_2nd,
            fine_lambda_min, fine_lambda_max, helmholtz_shift,
            pre_smooth_deg, post_smooth_deg, coarsest_solve_deg);

        size_t safe_states = static_cast< size_t >(std::max(1, my_max_states_));
        size_t grid_elems  = static_cast< size_t >(Gx_) * static_cast< size_t >(Gy_) * static_cast< size_t >(Gz_);
        size_t band_elems  = safe_states * grid_elems;

        (void)hipMalloc(&raw_psi_band_,    band_elems * sizeof(Wave64T));
        (void)hipMalloc(&raw_vnl_band_,    band_elems * sizeof(Wave64T));
        (void)hipMalloc(&raw_r_wave_,      band_elems * sizeof(Wave64T));
        (void)hipMalloc(&raw_p_wave_,      band_elems * sizeof(Wave64T));
        (void)hipMalloc(&raw_vloc_wave_,   grid_elems * sizeof(double));
        (void)hipMalloc(&d_vloc_32_,       grid_elems * sizeof(float));
        (void)hipMalloc(&raw_eig_wave_,    safe_states * sizeof(double));
        (void)hipMalloc(&d_eig_32_,        safe_states * sizeof(float));

        (void)hipMalloc(&d_r_32_,          band_elems * sizeof(LowT));
        (void)hipMalloc(&d_e_32_,          band_elems * sizeof(LowT));

        (void)hipMalloc(&d_diis_coeffs_,     safe_states * MAX_DIIS * sizeof(LowT));
        (void)hipMalloc(&d_dot_partials_,    safe_states * DOT_BLOCKS * 6 * sizeof(std::complex< double >));
        (void)hipMalloc(&d_norm_partials_,   safe_states * DOT_BLOCKS * 2 * sizeof(double));
        (void)hipMalloc(&raw_inv_norms_,     safe_states * sizeof(double));

        (void)hipMemset(raw_psi_band_,  0, band_elems * sizeof(Wave64T));
        (void)hipMemset(raw_vnl_band_,  0, band_elems * sizeof(Wave64T));
        (void)hipMemset(raw_r_wave_,    0, band_elems * sizeof(Wave64T));
        (void)hipMemset(raw_p_wave_,    0, band_elems * sizeof(Wave64T));
        (void)hipMemset(d_r_32_,        0, band_elems * sizeof(LowT));
        (void)hipMemset(d_e_32_,        0, band_elems * sizeof(LowT));
    }

    ~DynamicPrecisionGpuSolver() {
        (void)hipDeviceSynchronize();
        mg_.reset();
        transposer_32_.reset();
        transposer_64_.reset();

        for (auto& kv : block_diis_) {
            free_block_diis(kv.second);
        }
        block_diis_.clear();

        if (d_diis_coeffs_     != nullptr) (void)hipFree(d_diis_coeffs_);
        if (d_dot_partials_    != nullptr) (void)hipFree(d_dot_partials_);
        if (d_norm_partials_   != nullptr) (void)hipFree(d_norm_partials_);
        if (raw_inv_norms_     != nullptr) (void)hipFree(raw_inv_norms_);

        if (raw_psi_band_    != nullptr) (void)hipFree(raw_psi_band_);
        if (raw_vnl_band_    != nullptr) (void)hipFree(raw_vnl_band_);
        if (raw_r_wave_      != nullptr) (void)hipFree(raw_r_wave_);
        if (raw_p_wave_      != nullptr) (void)hipFree(raw_p_wave_);
        if (raw_vloc_wave_   != nullptr) (void)hipFree(raw_vloc_wave_);
        if (d_vloc_32_       != nullptr) (void)hipFree(d_vloc_32_);
        if (raw_eig_wave_    != nullptr) (void)hipFree(raw_eig_wave_);
        if (d_eig_32_        != nullptr) (void)hipFree(d_eig_32_);
        if (d_r_32_          != nullptr) (void)hipFree(d_r_32_);
        if (d_e_32_          != nullptr) (void)hipFree(d_e_32_);

        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized && comm_ != MPI_COMM_NULL) {
            MPI_Comm_free(&comm_);
        }
    }

    DynamicPrecisionGpuSolver(const DynamicPrecisionGpuSolver&) = delete;
    DynamicPrecisionGpuSolver& operator=(const DynamicPrecisionGpuSolver&) = delete;

    void reset_diis_history() {
        for (auto& kv : block_diis_) {
            kv.second.hist_len = 0;
        }
    }

    void update_stencil_taps(
        const std::vector< StencilTap< std::complex< double > > >& fine_taps_8th,
        const std::vector< std::vector< StencilTap< std::complex< double > > > >& coarse_taps_2nd) {
        fine_taps_ = fine_taps_8th;
        mg_->update_stencils(fine_taps_8th, coarse_taps_2nd);
        reset_diis_history();
    }

    template < typename VlocT >
    void upload_local_potential(const VlocT* h_vloc_padded, hipStream_t stream = nullptr) {
        size_t grid_elems = static_cast< size_t >(Gx_) * static_cast< size_t >(Gy_) * static_cast< size_t >(Gz_);
        (void)hipMemcpyAsync(raw_vloc_wave_, h_vloc_padded, grid_elems * sizeof(VlocT), hipMemcpyHostToDevice, stream);

        if constexpr (std::is_same_v< VlocT, float >) {
            (void)hipMemcpyAsync(d_vloc_32_, h_vloc_padded, grid_elems * sizeof(float), hipMemcpyHostToDevice, stream);
        } else {
            std::vector< float > h_v32(grid_elems);
            for (size_t i = 0; i < grid_elems; ++i) {
                h_v32.at(i) = static_cast< float >(*(h_vloc_padded + i));
            }
            (void)hipMemcpyAsync(d_vloc_32_, h_v32.data(), grid_elems * sizeof(float), hipMemcpyHostToDevice, stream);
            (void)hipStreamSynchronize(stream);
        }
        (void)hipStreamSynchronize(stream);
        reset_diis_history();
    }

    template < typename WaveT, typename EigHostT >
    void solve_wavefunction_block(
        WaveT* h_psi_block,
        const WaveT* h_vnl_block,
        std::vector< EigHostT >& global_eigvals,
        int block_state_offset,
        int cur_block_states,
        double eig_mix_fraction = 0.7,
        double eig_max_step = 1.0e6,
        double fg_step = 1.0,
        hipStream_t stream = nullptr) {
        static_assert(std::is_same_v< WaveT, Wave32T > || std::is_same_v< WaveT, Wave64T >,
                      "WaveT must match either 32-bit or 64-bit solver precision");
        using RealT = typename PrecisionTraits< WaveT >::RealT;

        assert(cur_block_states > 0 && cur_block_states <= max_block_states_);

        BlockDiisState& bstate = get_or_create_block_diis(block_state_offset);
        if (bstate.hist_len < MAX_DIIS) {
            bstate.hist_len += 1;
        } else {
            void* tmp_psi = bstate.raw_hist_psi.at(0);
            LowT* tmp_r   = bstate.d_hist_r.at(0);
            for (int m = 0; m < MAX_DIIS - 1; ++m) {
                size_t um = static_cast< size_t >(m);
                bstate.raw_hist_psi.at(um) = bstate.raw_hist_psi.at(um + 1);
                bstate.d_hist_r.at(um)     = bstate.d_hist_r.at(um + 1);
            }
            bstate.raw_hist_psi.at(MAX_DIIS - 1) = tmp_psi;
            bstate.d_hist_r.at(MAX_DIIS - 1)     = tmp_r;
            (void)hipMemcpyAsync(bstate.raw_hist_psi_ptrs, bstate.raw_hist_psi.data(),
                                 MAX_DIIS * sizeof(void*), hipMemcpyHostToDevice, stream);
(void)hipMemcpyAsync(bstate.d_hist_r_ptrs, bstate.d_hist_r.data(),
MAX_DIIS * sizeof(LowT*), hipMemcpyHostToDevice, stream);
}
// Partition cur_block_states across ranks identically to CartesianWavefunctionTransposer
int base = cur_block_states / n_ranks_;
int rem = cur_block_states % n_ranks_;
int my_slot_offset = rank_ * base + std::min(rank_, rem);
int my_active_states = base + (rank_ < rem ? 1 : 0);
std::vector< RealT > h_eig_wave(static_cast< size_t >(std::max(1, my_active_states)), RealT{0});
std::vector< float > h_eig_32(static_cast< size_t >(std::max(1, my_active_states)), 0.0f);
for (int b = 0; b < my_active_states; ++b) {
size_t ub = static_cast< size_t >(b);
double ev = static_cast< double >(global_eigvals.at(static_cast< size_t >(block_state_offset + my_slot_offset + b)));
h_eig_wave.at(ub) = static_cast< RealT >(ev);
h_eig_32.at(ub) = static_cast< float >(ev);
}
RealT* d_eig_wave = reinterpret_cast< RealT* >(raw_eig_wave_);
if (my_active_states > 0) {
(void)hipMemcpyAsync(d_eig_wave, h_eig_wave.data(),
static_cast< size_t >(my_active_states) * sizeof(RealT), hipMemcpyHostToDevice, stream);
(void)hipMemcpyAsync(d_eig_32_, h_eig_32.data(),
static_cast< size_t >(my_active_states) * sizeof(float), hipMemcpyHostToDevice, stream);
}
WaveT* d_psi = reinterpret_cast< WaveT* >(raw_psi_band_);
WaveT* d_vnl = reinterpret_cast< WaveT* >(raw_vnl_band_);
if constexpr (std::is_same_v< WaveT, Wave32T >) {
transposer_32_->forward_pair(h_psi_block, h_vnl_block, d_psi, d_vnl);
} else {
transposer_64_->forward_pair(h_psi_block, h_vnl_block, d_psi, d_vnl);
}
if (my_active_states > 0) {
std::array< WaveT*, 3 > hist_psi_typed = {{
reinterpret_cast< WaveT* >(bstate.raw_hist_psi.at(0)),
reinterpret_cast< WaveT* >(bstate.raw_hist_psi.at(1)),
reinterpret_cast< WaveT* >(bstate.raw_hist_psi.at(2))
}};
run_single_defect_correction_gpu< WaveT >(
d_psi, d_vnl,
reinterpret_cast< const RealT* >(raw_vloc_wave_), d_eig_wave,
d_eig_32_, h_eig_wave, h_eig_32,
reinterpret_cast< WaveT* >(raw_r_wave_),
reinterpret_cast< WaveT* >(raw_p_wave_),
d_r_32_, d_e_32_,
hist_psi_typed, bstate.d_hist_r,
reinterpret_cast< WaveT* const* >(bstate.raw_hist_psi_ptrs),
bstate.d_hist_r_ptrs, d_diis_coeffs_,
d_dot_partials_, d_norm_partials_,
reinterpret_cast< RealT* >(raw_inv_norms_),
*mg_, fine_taps_, bstate.hist_len,
my_active_states, Nx_, Ny_, Nz_, Gx_, Gy_, Gz_,
vel_, periodic_, eig_mix_fraction, eig_max_step, fg_step, stream);
}
if constexpr (std::is_same_v< WaveT, Wave32T >) {
transposer_32_->inverse(d_psi, h_psi_block);
} else {
transposer_64_->inverse(d_psi, h_psi_block);
}
std::vector< int > counts(static_cast< size_t >(n_ranks_), 0);
std::vector< int > displs(static_cast< size_t >(n_ranks_), 0);
for (int r = 0; r < n_ranks_; ++r) {
size_t ur = static_cast< size_t >(r);
counts.at(ur) = base + (r < rem ? 1 : 0);
displs.at(ur) = r * base + std::min(r, rem);
}
std::vector< RealT > block_eigs_gathered(static_cast< size_t >(cur_block_states), RealT{0});
MPI_Allgatherv(h_eig_wave.data(), my_active_states, get_mpi_type< RealT >(),
block_eigs_gathered.data(), counts.data(), displs.data(),
get_mpi_type< RealT >(), comm_);
for (int s = 0; s < cur_block_states; ++s) {
global_eigvals.at(static_cast< size_t >(block_state_offset + s)) =
static_cast< EigHostT >(block_eigs_gathered.at(static_cast< size_t >(s)));
}
}
int max_block_states() const { return max_block_states_; }
private:
BlockDiisState& get_or_create_block_diis(int block_state_offset) {
auto it = block_diis_.find(block_state_offset);
if (it != block_diis_.end()) {
return it->second;
}
BlockDiisState st{};
size_t safe_states = static_cast< size_t >(std::max(1, my_max_states_));
size_t grid_elems = static_cast< size_t >(Gx_) * static_cast< size_t >(Gy_) * static_cast< size_t >(Gz_);
size_t band_elems = safe_states * grid_elems;
for (int m = 0; m < MAX_DIIS; ++m) {
size_t um = static_cast< size_t >(m);
(void)hipMalloc(&st.raw_hist_psi.at(um), band_elems * sizeof(Wave64T));
(void)hipMalloc(&st.d_hist_r.at(um), band_elems * sizeof(LowT));
(void)hipMemset(st.raw_hist_psi.at(um), 0, band_elems * sizeof(Wave64T));
(void)hipMemset(st.d_hist_r.at(um), 0, band_elems * sizeof(LowT));
}
(void)hipMalloc(&st.raw_hist_psi_ptrs, MAX_DIIS * sizeof(void*));
(void)hipMalloc(&st.d_hist_r_ptrs, MAX_DIIS * sizeof(LowT*));
(void)hipMemcpy(st.raw_hist_psi_ptrs, st.raw_hist_psi.data(), MAX_DIIS * sizeof(void*), hipMemcpyHostToDevice);
(void)hipMemcpy(st.d_hist_r_ptrs, st.d_hist_r.data(), MAX_DIIS * sizeof(LowT*), hipMemcpyHostToDevice);
auto inserted = block_diis_.emplace(block_state_offset, st);
return inserted.first->second;
}
void free_block_diis(BlockDiisState& st) {
for (int m = 0; m < MAX_DIIS; ++m) {
size_t um = static_cast< size_t >(m);
if (st.raw_hist_psi.at(um) != nullptr) (void)hipFree(st.raw_hist_psi.at(um));
if (st.d_hist_r.at(um) != nullptr) (void)hipFree(st.d_hist_r.at(um));
}
if (st.raw_hist_psi_ptrs != nullptr) (void)hipFree(st.raw_hist_psi_ptrs);
if (st.d_hist_r_ptrs != nullptr) (void)hipFree(st.d_hist_r_ptrs);
}
MPI_Comm comm_ = MPI_COMM_NULL;
int n_ranks_ = 0;
int rank_ = 0;
int max_block_states_ = 0;
int my_max_states_ = 0;
int Nx_ = 0, Ny_ = 0, Nz_ = 0;
int Gx_ = 0, Gy_ = 0, Gz_ = 0;
double vel_ = 1.0;
bool periodic_ = true;
std::map< int, BlockDiisState > block_diis_;
std::vector< StencilTap< std::complex< double > > > fine_taps_;
std::unique_ptr< CartesianWavefunctionTransposer< Wave32T > > transposer_32_;
std::unique_ptr< CartesianWavefunctionTransposer< Wave64T > > transposer_64_;
std::unique_ptr< MultigridPreconditioner< LowT > > mg_;
void* raw_psi_band_ = nullptr;
void* raw_vnl_band_ = nullptr;
void* raw_r_wave_ = nullptr;
void* raw_p_wave_ = nullptr;
void* raw_vloc_wave_ = nullptr;
float* d_vloc_32_ = nullptr;
void* raw_eig_wave_ = nullptr;
float* d_eig_32_ = nullptr;
LowT* d_r_32_ = nullptr;
LowT* d_e_32_ = nullptr;
LowT* d_diis_coeffs_ = nullptr;
std::complex< double >* d_dot_partials_ = nullptr;
double* d_norm_partials_ = nullptr;
void* raw_inv_norms_ = nullptr;
};

