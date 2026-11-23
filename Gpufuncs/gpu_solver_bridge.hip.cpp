#include "gpu_solver_bridge.hpp"
#include "dynamic_gpu_solver.hpp"
#include <cassert>
#include "main.h"
#include "transition.h"

namespace {

inline std::vector< StencilTap< std::complex< double > > > to_internal_taps(
    const GpuStencilTap* taps,
    int num_taps) {
    assert(num_taps >= 0);
    std::vector< StencilTap< std::complex< double > > > out;
    out.reserve(static_cast< size_t >(num_taps));
    for (int i = 0; i < num_taps; ++i) {
        const GpuStencilTap& t = *(taps + i);
        out.push_back(StencilTap< std::complex< double > >{t.dx, t.dy, t.dz, t.weight});
    }
    return out;
}

inline std::vector< std::vector< StencilTap< std::complex< double > > > > to_internal_coarse_taps(
    int mg_num_levels,
    const GpuStencilTap* const* coarse_taps,
    const int* num_coarse_taps_per_level) {
    int num_coarse = (mg_num_levels > 1) ? (mg_num_levels - 1) : 0;
    std::vector< std::vector< StencilTap< std::complex< double > > > > out;
    out.reserve(static_cast< size_t >(num_coarse));
    for (int l = 0; l < num_coarse; ++l) {
        out.push_back(to_internal_taps(*(coarse_taps + l), *(num_coarse_taps_per_level + l)));
    }
    return out;
}

} // namespace

template < typename Canonical64T >
struct GpuSolverBridge< Canonical64T >::Impl {
    int mg_num_levels = 1;
    DynamicPrecisionGpuSolver< Canonical64T > solver;

    Impl(MPI_Comm parent_comm,
         int max_block_states,
         int Nx, int Ny, int Nz,
         GpuSubBox3D my_box,
         int num_levels,
         const GpuStencilTap* fine_taps_8th,
         int num_fine_taps,
         const GpuStencilTap* const* coarse_taps_2nd,
         const int* num_coarse_taps_per_level,
         double fine_lambda_min,
         double fine_lambda_max,
         double helmholtz_shift,
         double vel,
         bool periodic,
         int desired_mpi_batch,
         int pre_smooth_deg,
         int post_smooth_deg,
         int coarsest_solve_deg)
        : mg_num_levels(num_levels),
          solver(parent_comm, max_block_states, Nx, Ny, Nz,
                 SubBox3D{my_box.x0, my_box.y0, my_box.z0, my_box.nx, my_box.ny, my_box.nz},
                 num_levels,
                 to_internal_taps(fine_taps_8th, num_fine_taps),
                 to_internal_coarse_taps(num_levels, coarse_taps_2nd, num_coarse_taps_per_level),
                 fine_lambda_min, fine_lambda_max, helmholtz_shift,
                 vel, periodic, desired_mpi_batch, pre_smooth_deg, post_smooth_deg, coarsest_solve_deg) {}
};

template < typename Canonical64T >
GpuSolverBridge< Canonical64T >::GpuSolverBridge(
    MPI_Comm parent_comm,
    int max_block_states,
    int Nx, int Ny, int Nz,
    GpuSubBox3D my_box,
    int mg_num_levels,
    const GpuStencilTap* fine_taps_8th,
    int num_fine_taps,
    const GpuStencilTap* const* coarse_taps_2nd,
    const int* num_coarse_taps_per_level,
    double fine_lambda_min,
    double fine_lambda_max,
    double helmholtz_shift,
    double vel,
    bool periodic,
    int desired_mpi_batch,
    int pre_smooth_deg,
    int post_smooth_deg,
    int coarsest_solve_deg)
    : impl_(std::make_unique< Impl >(
          parent_comm, max_block_states, Nx, Ny, Nz, my_box, mg_num_levels,
          fine_taps_8th, num_fine_taps, coarse_taps_2nd, num_coarse_taps_per_level,
          fine_lambda_min, fine_lambda_max, helmholtz_shift,
          vel, periodic, desired_mpi_batch, pre_smooth_deg, post_smooth_deg, coarsest_solve_deg)) {
}

template < typename Canonical64T >
GpuSolverBridge< Canonical64T >::~GpuSolverBridge() = default;

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::update_stencil_taps(
    const GpuStencilTap* fine_taps_8th,
    int num_fine_taps,
    const GpuStencilTap* const* coarse_taps_2nd,
    const int* num_coarse_taps_per_level) {
    impl_->solver.update_stencil_taps(
        to_internal_taps(fine_taps_8th, num_fine_taps),
        to_internal_coarse_taps(impl_->mg_num_levels, coarse_taps_2nd, num_coarse_taps_per_level));
}

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::upload_local_potential_32(const float* h_vloc_padded) {
    impl_->solver.upload_local_potential(h_vloc_padded);
}

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::upload_local_potential_64(const double* h_vloc_padded) {
    impl_->solver.upload_local_potential(h_vloc_padded);
}

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::solve_block_32(
    Wave32T* h_psi_block,
    const Wave32T* h_vnl_block,
    double* global_eigvals,
    int total_states,
    int block_state_offset,
    int cur_block_states,
    double eig_mix_fraction,
    double eig_max_step) {
    (void)eig_max_step;
    double eig_mix = (ct.scf_steps && (ct.scf_accuracy < 1.0e-4)) ? 1.0 : eig_mix_fraction;
    std::vector< double > eig_vec(global_eigvals, global_eigvals + total_states);
    impl_->solver.template solve_wavefunction_block< Wave32T >(
        h_psi_block, h_vnl_block, eig_vec,
        block_state_offset, cur_block_states, eig_mix, 1.0e6);
    for (int s = 0; s < cur_block_states; ++s) {
        int idx = block_state_offset + s;
        *(global_eigvals + idx) = eig_vec.at(static_cast< size_t >(idx));
    }
}

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::solve_block_64(
    Wave64T* h_psi_block,
    const Wave64T* h_vnl_block,
    double* global_eigvals,
    int total_states,
    int block_state_offset,
    int cur_block_states,
    double eig_mix_fraction,
    double eig_max_step) {
    (void)eig_max_step;
    std::vector< double > eig_vec(global_eigvals, global_eigvals + total_states);
    double eig_mix = (ct.scf_steps && (ct.scf_accuracy < 1.0e-4)) ? 1.0 : eig_mix_fraction;
    impl_->solver.template solve_wavefunction_block< Wave64T >(
        h_psi_block, h_vnl_block, eig_vec,
        block_state_offset, cur_block_states, eig_mix, 1.0e6);
    for (int s = 0; s < cur_block_states; ++s) {
        int idx = block_state_offset + s;
        *(global_eigvals + idx) = eig_vec.at(static_cast< size_t >(idx));
    }
}

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::pin_host_memory(void* ptr, std::size_t bytes) {
    (void)hipHostRegister(ptr, bytes, hipHostRegisterDefault);
}

template < typename Canonical64T >
void GpuSolverBridge< Canonical64T >::unpin_host_memory(void* ptr) {
    (void)hipHostUnregister(ptr);
}

template class GpuSolverBridge< double >;
template class GpuSolverBridge< std::complex< double > >;

