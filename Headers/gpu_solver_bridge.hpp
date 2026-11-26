#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <mpi.h>
#include <type_traits>
#include <vector>

struct GpuSubBox3D {
    int x0 = 0, y0 = 0, z0 = 0;
    int nx = 0, ny = 0, nz = 0;
};

struct GpuStencilTap {
    int dx = 0;
    int dy = 0;
    int dz = 0;
    std::complex<double> weight{0.0, 0.0};
};

// Pure host C++ Pimpl wrapper (zero HIP headers and zero std::vector across the
// .cpp ABI boundary) Instantiate with Canonical64T = double (Gamma/real) or
// std::complex< double > (k-point/complex)
template <typename Canonical64T> class GpuSolverBridge {
  public:
    using Wave32T =
        typename std::conditional_t<std::is_same_v<Canonical64T, double>, float,
                                    std::complex<float>>;
    using Wave64T = Canonical64T;

    // Primary ABI-safe constructor (raw POD pointers + counts)
    GpuSolverBridge(MPI_Comm parent_comm, int max_block_states, int Nx, int Ny,
                    int Nz, GpuSubBox3D my_box, int mg_num_levels,
                    const GpuStencilTap *fine_taps_8th, int num_fine_taps,
                    const GpuStencilTap *const *coarse_taps_2nd,
                    const int *num_coarse_taps_per_level,
                    double fine_lambda_min, double fine_lambda_max,
                    double helmholtz_shift, double vel = 1.0,
                    bool periodic = true, int desired_mpi_batch = 1024,
                    int pre_smooth_deg = 2, int post_smooth_deg = 2,
                    int coarsest_solve_deg = 6);

    // Header-only convenience constructor forwarding std::vector into the
    // ABI-safe constructor
    GpuSolverBridge(
        MPI_Comm parent_comm, int max_block_states, int Nx, int Ny, int Nz,
        GpuSubBox3D my_box, int mg_num_levels,
        const std::vector<GpuStencilTap> &fine_taps_8th,
        const std::vector<std::vector<GpuStencilTap>> &coarse_taps_2nd,
        double fine_lambda_min, double fine_lambda_max, double helmholtz_shift,
        double vel = 1.0, bool periodic = true, int desired_mpi_batch = 1024,
        int pre_smooth_deg = 2, int post_smooth_deg = 2,
        int coarsest_solve_deg = 6)
        : GpuSolverBridge(
              parent_comm, max_block_states, Nx, Ny, Nz, my_box, mg_num_levels,
              fine_taps_8th.data(), static_cast<int>(fine_taps_8th.size()),
              [&]() {
                  std::vector<const GpuStencilTap *> ptrs;
                  ptrs.reserve(coarse_taps_2nd.size());
                  for (const auto &lvl : coarse_taps_2nd)
                      ptrs.push_back(lvl.data());
                  return ptrs;
              }()
                  .data(),
              [&]() {
                  std::vector<int> cnts;
                  cnts.reserve(coarse_taps_2nd.size());
                  for (const auto &lvl : coarse_taps_2nd)
                      cnts.push_back(static_cast<int>(lvl.size()));
                  return cnts;
              }()
                  .data(),
              fine_lambda_min, fine_lambda_max, helmholtz_shift, vel, periodic,
              desired_mpi_batch, pre_smooth_deg, post_smooth_deg,
              coarsest_solve_deg) {}

    ~GpuSolverBridge();

    GpuSolverBridge(const GpuSolverBridge &) = delete;
    GpuSolverBridge &operator=(const GpuSolverBridge &) = delete;

    // ABI-safe k-point stencil tap update (copies data internally; caller
    // buffers may be reused immediately)
    void update_stencil_taps(const GpuStencilTap *fine_taps_8th,
                             int num_fine_taps,
                             const GpuStencilTap *const *coarse_taps_2nd,
                             const int *num_coarse_taps_per_level);

    // Header-only convenience overload accepting std::vector
    void update_stencil_taps(
        const std::vector<GpuStencilTap> &fine_taps_8th,
        const std::vector<std::vector<GpuStencilTap>> &coarse_taps_2nd) {
        std::vector<const GpuStencilTap *> ptrs;
        std::vector<int> cnts;
        ptrs.reserve(coarse_taps_2nd.size());
        cnts.reserve(coarse_taps_2nd.size());
        for (const auto &lvl : coarse_taps_2nd) {
            ptrs.push_back(lvl.data());
            cnts.push_back(static_cast<int>(lvl.size()));
        }
        update_stencil_taps(fine_taps_8th.data(),
                            static_cast<int>(fine_taps_8th.size()), ptrs.data(),
                            cnts.data());
    }

    // Upload padded local potential V_loc [Nx+8][Ny+8][Nz+8] once per SCF step
    void upload_local_potential_32(const float *h_vloc_padded);
    void upload_local_potential_64(const double *h_vloc_padded);

    template <typename VlocT>
    void upload_local_potential(const VlocT *h_vloc_padded) {
        if constexpr (std::is_same_v<VlocT, float>) {
            upload_local_potential_32(h_vloc_padded);
        } else {
            upload_local_potential_64(h_vloc_padded);
        }
    }

    // ABI-safe blocked solves taking raw double* global_eigvals pointer
    void solve_block_32(Wave32T *h_psi_block, const Wave32T *h_vnl_block,
                        double *global_eigvals, int total_states,
                        int block_state_offset, int cur_block_states,
                        double eig_mix_fraction = 0.7,
                        double eig_max_step = 0.25);

    void solve_block_64(Wave64T *h_psi_block, const Wave64T *h_vnl_block,
                        double *global_eigvals, int total_states,
                        int block_state_offset, int cur_block_states,
                        double eig_mix_fraction = 0.7,
                        double eig_max_step = 0.25);

    // Header-only convenience overloads accepting std::vector< double >&
    // global_eigvals
    void solve_block_32(Wave32T *h_psi_block, const Wave32T *h_vnl_block,
                        std::vector<double> &global_eigvals,
                        int block_state_offset, int cur_block_states,
                        double eig_mix_fraction = 0.7,
                        double eig_max_step = 0.25) {
        solve_block_32(h_psi_block, h_vnl_block, global_eigvals.data(),
                       static_cast<int>(global_eigvals.size()),
                       block_state_offset, cur_block_states, eig_mix_fraction,
                       eig_max_step);
    }

    void solve_block_64(Wave64T *h_psi_block, const Wave64T *h_vnl_block,
                        std::vector<double> &global_eigvals,
                        int block_state_offset, int cur_block_states,
                        double eig_mix_fraction = 0.7,
                        double eig_max_step = 0.25) {
        solve_block_64(h_psi_block, h_vnl_block, global_eigvals.data(),
                       static_cast<int>(global_eigvals.size()),
                       block_state_offset, cur_block_states, eig_mix_fraction,
                       eig_max_step);
    }

    template <typename KpointType>
    void solve_block(KpointType *h_psi_block, const KpointType *h_vnl_block,
                     std::vector<double> &global_eigvals,
                     int block_state_offset, int cur_block_states,
                     double eig_mix_fraction = 0.7,
                     double eig_max_step = 0.25) {
        static_assert(std::is_same_v<KpointType, Wave32T> ||
                          std::is_same_v<KpointType, Wave64T>,
                      "KpointType must match either the 32-bit or 64-bit "
                      "wavefunction type of this bridge");
        if constexpr (std::is_same_v<KpointType, Wave32T>) {
            solve_block_32(h_psi_block, h_vnl_block, global_eigvals,
                           block_state_offset, cur_block_states,
                           eig_mix_fraction, eig_max_step);
        } else {
            solve_block_64(h_psi_block, h_vnl_block, global_eigvals,
                           block_state_offset, cur_block_states,
                           eig_mix_fraction, eig_max_step);
        }
    }

    template <typename KpointType>
    void solve_block(KpointType *h_psi_block, const KpointType *h_vnl_block,
                     double *global_eigvals, int total_states,
                     int block_state_offset, int cur_block_states,
                     double eig_mix_fraction = 0.7,
                     double eig_max_step = 0.25) {
        static_assert(std::is_same_v<KpointType, Wave32T> ||
                          std::is_same_v<KpointType, Wave64T>,
                      "KpointType must match either the 32-bit or 64-bit "
                      "wavefunction type of this bridge");
        if constexpr (std::is_same_v<KpointType, Wave32T>) {
            solve_block_32(h_psi_block, h_vnl_block, global_eigvals,
                           total_states, block_state_offset, cur_block_states,
                           eig_mix_fraction, eig_max_step);
        } else {
            solve_block_64(h_psi_block, h_vnl_block, global_eigvals,
                           total_states, block_state_offset, cur_block_states,
                           eig_mix_fraction, eig_max_step);
        }
    }

    static void pin_host_memory(void *ptr, std::size_t bytes);
    static void unpin_host_memory(void *ptr);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
