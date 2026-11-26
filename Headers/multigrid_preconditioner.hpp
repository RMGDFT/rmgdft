#pragma once

#include "chebyshev_rmmdiis.hpp"
#include "custom_stencils.hpp"
#include "multigrid_transfer.hpp"
#include "periodic_halo.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <hip/hip_runtime.h>
#include <vector>

struct MultigridLevelConfig {
    int Nx = 0, Ny = 0, Nz = 0;
    int Gx = 0, Gy = 0, Gz = 0;
    int halo = 1;
    float lambda_min = 0.1f;
    float lambda_max = 10.0f;
    int pre_smooth_deg = 2;
    int post_smooth_deg = 2;
};

template <typename LowT> class MultigridPreconditioner {
  public:
    using WeightT = typename PrecisionTraits<LowT>::WeightT;

    MultigridPreconditioner(
        int my_band_states, int Nx0, int Ny0, int Nz0, int num_levels,
        const std::vector<StencilTap<std::complex<double>>> &fine_taps_8th,
        const std::vector<std::vector<StencilTap<std::complex<double>>>>
            &coarse_taps_2nd_per_level,
        double fine_lambda_min, double fine_lambda_max, double helmholtz_shift,
        int pre_smooth_deg = 2, int post_smooth_deg = 2,
        int coarsest_solve_deg = 6)
        : my_band_states_(my_band_states), num_levels_(num_levels),
          coarsest_solve_deg_(coarsest_solve_deg),
          helmholtz_shift_(helmholtz_shift) {
        assert(num_levels_ >= 1);
        assert(static_cast<int>(coarse_taps_2nd_per_level.size()) >=
               num_levels_ - 1);

        levels_.resize(static_cast<size_t>(num_levels_));
        d_e_.resize(static_cast<size_t>(num_levels_), nullptr);
        d_rhs_.resize(static_cast<size_t>(num_levels_), nullptr);

        int cur_nx = Nx0, cur_ny = Ny0, cur_nz = Nz0;
        size_t max_state_elems = 0;

        for (int l = 0; l < num_levels_; ++l) {
            size_t ul = static_cast<size_t>(l);
            int h = (l == 0) ? 4 : 1;

            MultigridLevelConfig cfg{};
            cfg.Nx = cur_nx;
            cfg.Ny = cur_ny;
            cfg.Nz = cur_nz;
            cfg.halo = h;
            cfg.Gx = cur_nx + 2 * h;
            cfg.Gy = cur_ny + 2 * h;
            cfg.Gz = cur_nz + 2 * h;
            cfg.pre_smooth_deg =
                (l == num_levels_ - 1) ? coarsest_solve_deg_ : pre_smooth_deg;
            cfg.post_smooth_deg = (l == num_levels_ - 1) ? 0 : post_smooth_deg;
            levels_.at(ul) = cfg;

            size_t grid_elems = static_cast<size_t>(cfg.Gx) *
                                static_cast<size_t>(cfg.Gy) *
                                static_cast<size_t>(cfg.Gz);
            size_t band_elems = std::max<size_t>(
                1, static_cast<size_t>(my_band_states_) * grid_elems);
            if (l == 0)
                max_state_elems = band_elems;

            if (l > 0) {
                (void)hipMalloc(&d_e_.at(ul), band_elems * sizeof(LowT));
                (void)hipMalloc(&d_rhs_.at(ul), band_elems * sizeof(LowT));
                (void)hipMemset(d_e_.at(ul), 0, band_elems * sizeof(LowT));
                (void)hipMemset(d_rhs_.at(ul), 0, band_elems * sizeof(LowT));
            }

            cur_nx = (cur_nx + 1) / 2;
            cur_ny = (cur_ny + 1) / 2;
            cur_nz = (cur_nz + 1) / 2;
        }

        levels_.at(0).lambda_min = static_cast<float>(fine_lambda_min);
        levels_.at(0).lambda_max = static_cast<float>(fine_lambda_max);

        size_t lvl1_grid_elems =
            (num_levels_ > 1) ? static_cast<size_t>(levels_.at(1).Gx) *
                                    static_cast<size_t>(levels_.at(1).Gy) *
                                    static_cast<size_t>(levels_.at(1).Gz)
                              : 1;
        size_t safe_states = static_cast<size_t>(std::max(1, my_band_states_));
        (void)hipMalloc(&d_zero_vloc_, lvl1_grid_elems * sizeof(float));
        (void)hipMalloc(&d_zero_eig_, safe_states * sizeof(float));
        (void)hipMemset(d_zero_vloc_, 0, lvl1_grid_elems * sizeof(float));
        (void)hipMemset(d_zero_eig_, 0, safe_states * sizeof(float));

        (void)hipMalloc(&d_cheb_d_, max_state_elems * sizeof(LowT));
        (void)hipMalloc(&d_cheb_Au_, max_state_elems * sizeof(LowT));
        (void)hipMalloc(&d_res_scratch_, max_state_elems * sizeof(LowT));
        (void)hipMemset(d_cheb_d_, 0, max_state_elems * sizeof(LowT));
        (void)hipMemset(d_cheb_Au_, 0, max_state_elems * sizeof(LowT));
        (void)hipMemset(d_res_scratch_, 0, max_state_elems * sizeof(LowT));

        coarse_stencils_.resize(
            static_cast<size_t>(std::max(0, num_levels_ - 1)));
        update_stencils(fine_taps_8th, coarse_taps_2nd_per_level);
    }

    ~MultigridPreconditioner() {
        (void)hipDeviceSynchronize();
        for (int l = 1; l < num_levels_; ++l) {
            size_t ul = static_cast<size_t>(l);
            if (d_e_.at(ul) != nullptr)
                (void)hipFree(d_e_.at(ul));
            if (d_rhs_.at(ul) != nullptr)
                (void)hipFree(d_rhs_.at(ul));
        }
        if (d_zero_vloc_ != nullptr)
            (void)hipFree(d_zero_vloc_);
        if (d_zero_eig_ != nullptr)
            (void)hipFree(d_zero_eig_);
        if (d_cheb_d_ != nullptr)
            (void)hipFree(d_cheb_d_);
        if (d_cheb_Au_ != nullptr)
            (void)hipFree(d_cheb_Au_);
        if (d_res_scratch_ != nullptr)
            (void)hipFree(d_res_scratch_);
    }

    MultigridPreconditioner(const MultigridPreconditioner &) = delete;
    MultigridPreconditioner &
    operator=(const MultigridPreconditioner &) = delete;

    // Fast k-point stencil update reusing initial fine_lambda_min,
    // fine_lambda_max, and helmholtz_shift
    void update_stencils(
        const std::vector<StencilTap<std::complex<double>>> &fine_taps_8th,
        const std::vector<std::vector<StencilTap<std::complex<double>>>>
            &coarse_taps_2nd_per_level) {
        assert(static_cast<int>(coarse_taps_2nd_per_level.size()) >=
               num_levels_ - 1);

        fine_stencil_ = make_custom_stencil<WeightT, 4, 104>(fine_taps_8th);

        for (int l = 1; l < num_levels_; ++l) {
            size_t ul = static_cast<size_t>(l);
            size_t idx = ul - 1;

            std::vector<StencilTap<std::complex<double>>> helm_taps =
                coarse_taps_2nd_per_level.at(idx);
            bool found_center = false;
            double off_diag_abs_sum = 0.0;
            double diag_real = 0.0;

            for (auto &tap : helm_taps) {
                if (tap.dx == 0 && tap.dy == 0 && tap.dz == 0) {
                    tap.weight += std::complex<double>(helmholtz_shift_, 0.0);
                    diag_real += tap.weight.real();
                    found_center = true;
                } else {
                    off_diag_abs_sum += std::abs(tap.weight);
                }
            }
            if (!found_center) {
                helm_taps.insert(
                    helm_taps.begin(),
                    StencilTap<std::complex<double>>{
                        0, 0, 0, std::complex<double>(helmholtz_shift_, 0.0)});
                diag_real += helmholtz_shift_;
            }

            coarse_stencils_.at(idx) =
                make_custom_stencil<WeightT, 1, 20>(helm_taps);

            double gershgorin_max =
                1.05 * std::max(1e-6, diag_real + off_diag_abs_sum);
#if 0
     // This doesn't work as well as the simple update below
            double gershgorin_min =
                (l == num_levels_ - 1)
                    ? std::max(1e-6, helmholtz_shift_)
                    : std::max(helmholtz_shift_, 0.25 * gershgorin_max);
#endif
            double gershgorin_min = 0.1*gershgorin_max;

            levels_.at(ul).lambda_min = static_cast<float>(gershgorin_min);
            levels_.at(ul).lambda_max = static_cast<float>(gershgorin_max);
        }
    }

    void smooth_fine_level(LowT *d_e0, const LowT *d_r0, const float *d_vloc_0,
                           const float *d_eigvals_fp32, int degree,
                           bool zero_initial_guess, bool periodic,
                           hipStream_t stream) {
        if (my_band_states_ <= 0 || degree <= 0)
            return;
        const auto &fine = levels_.at(0);
        chebyshev_smooth_batched<LowT, 4, 104>(
            d_e0, d_cheb_d_, d_cheb_Au_, d_r0, d_vloc_0, d_eigvals_fp32,
            my_band_states_, fine.Nx, fine.Ny, fine.Nz, fine.Gx, fine.Gy,
            fine.Gz, fine_stencil_, fine.lambda_min, fine.lambda_max, degree,
            zero_initial_guess, periodic, stream);
    }

    void coarse_helmholtz_correction(LowT *d_e0, const LowT *d_r0,
                                     bool periodic, hipStream_t stream) {
        if (my_band_states_ <= 0 || num_levels_ <= 1)
            return;

        // 1. Restrict fine-grid residual d_r0 -> d_rhs_[1]
        {
            const auto &fine = levels_.at(0);
            const auto &crse = levels_.at(1);
            launch_fill_periodic_halo<LowT, 4>(
                const_cast<LowT *>(d_r0), my_band_states_, fine.Nx, fine.Ny,
                fine.Nz, fine.Gx, fine.Gy, fine.Gz, periodic, stream);
            launch_restrict_3d<LowT, 4, 1>(
                d_r0, d_rhs_.at(1), my_band_states_, fine.Nx, fine.Ny, fine.Nz,
                fine.Gx, fine.Gy, fine.Gz, crse.Nx, crse.Ny, crse.Nz, crse.Gx,
                crse.Gy, crse.Gz, stream);
        }

        // 2. Downward sweep across coarse levels 1 .. num_levels_ - 2
        for (int l = 1; l < num_levels_ - 1; ++l) {
            size_t ul = static_cast<size_t>(l);
            size_t un = ul + 1;
            const auto &fine = levels_.at(ul);
            const auto &crse = levels_.at(un);
            const auto &st_c = coarse_stencils_.at(ul - 1);

            chebyshev_smooth_batched<LowT, 1, 20>(
                d_e_.at(ul), d_cheb_d_, d_cheb_Au_, d_rhs_.at(ul), d_zero_vloc_,
                d_zero_eig_, my_band_states_, fine.Nx, fine.Ny, fine.Nz,
                fine.Gx, fine.Gy, fine.Gz, st_c, fine.lambda_min,
                fine.lambda_max, fine.pre_smooth_deg, true, periodic, stream);

            launch_shifted_hamiltonian<LowT, 1, 20>(
                d_e_.at(ul), d_cheb_Au_, d_zero_vloc_, d_zero_eig_,
                my_band_states_, fine.Nx, fine.Ny, fine.Nz, fine.Gx, fine.Gy,
                fine.Gz, st_c, periodic, stream);
            compute_residual_vec<LowT, 1>(d_rhs_.at(ul), d_cheb_Au_,
                                          d_res_scratch_, fine, stream);

            launch_fill_periodic_halo<LowT, 1>(
                d_res_scratch_, my_band_states_, fine.Nx, fine.Ny, fine.Nz,
                fine.Gx, fine.Gy, fine.Gz, periodic, stream);
            launch_restrict_3d<LowT, 1, 1>(
                d_res_scratch_, d_rhs_.at(un), my_band_states_, fine.Nx,
                fine.Ny, fine.Nz, fine.Gx, fine.Gy, fine.Gz, crse.Nx, crse.Ny,
                crse.Nz, crse.Gx, crse.Gy, crse.Gz, stream);
        }
        // 3. Coarsest Grid Helmholtz Solve (Level num_levels_ - 1)
        {
            size_t ulast = static_cast<size_t>(num_levels_ - 1);
            const auto &clvl = levels_.at(ulast);
            const auto &st_last = coarse_stencils_.at(ulast - 1);
            chebyshev_smooth_batched<LowT, 1, 20>(
                d_e_.at(ulast), d_cheb_d_, d_cheb_Au_, d_rhs_.at(ulast),
                d_zero_vloc_, d_zero_eig_, my_band_states_, clvl.Nx, clvl.Ny,
                clvl.Nz, clvl.Gx, clvl.Gy, clvl.Gz, st_last, clvl.lambda_min,
                clvl.lambda_max, clvl.pre_smooth_deg, true, periodic, stream);
        }
        // 4. Upward sweep across coarse levels num_levels_ - 2 .. 1
        for (int l = num_levels_ - 2; l >= 1; --l) {
            size_t ul = static_cast<size_t>(l);
            size_t un = ul + 1;
            const auto &fine = levels_.at(ul);
            const auto &crse = levels_.at(un);
            const auto &st_c = coarse_stencils_.at(ul - 1);
            launch_fill_periodic_halo<LowT, 1>(
                d_e_.at(un), my_band_states_, crse.Nx, crse.Ny, crse.Nz,
                crse.Gx, crse.Gy, crse.Gz, periodic, stream);
            launch_prolongate_add_3d<LowT, 1, 1>(
                d_e_.at(un), d_e_.at(ul), my_band_states_, fine.Nx, fine.Ny,
                fine.Nz, fine.Gx, fine.Gy, fine.Gz, crse.Nx, crse.Ny, crse.Nz,
                crse.Gx, crse.Gy, crse.Gz, stream);
            chebyshev_smooth_batched<LowT, 1, 20>(
                d_e_.at(ul), d_cheb_d_, d_cheb_Au_, d_rhs_.at(ul), d_zero_vloc_,
                d_zero_eig_, my_band_states_, fine.Nx, fine.Ny, fine.Nz,
                fine.Gx, fine.Gy, fine.Gz, st_c, fine.lambda_min,
                fine.lambda_max, fine.post_smooth_deg, false, periodic, stream);
        }
        // 5. Prolongate Level 1 Helmholtz error correction directly into d_e0
        {
            const auto &fine = levels_.at(0);
            const auto &crse = levels_.at(1);
            launch_fill_periodic_halo<LowT, 1>(
                d_e_.at(1), my_band_states_, crse.Nx, crse.Ny, crse.Nz, crse.Gx,
                crse.Gy, crse.Gz, periodic, stream);
            size_t elems0 = static_cast<size_t>(my_band_states_) *
                            static_cast<size_t>(fine.Gx) *
                            static_cast<size_t>(fine.Gy) *
                            static_cast<size_t>(fine.Gz);
            (void)hipMemsetAsync(d_e0, 0, elems0 * sizeof(LowT), stream);
            launch_prolongate_add_3d<LowT, 4, 1>(
                d_e_.at(1), d_e0, my_band_states_, fine.Nx, fine.Ny, fine.Nz,
                fine.Gx, fine.Gy, fine.Gz, crse.Nx, crse.Ny, crse.Nz, crse.Gx,
                crse.Gy, crse.Gz, stream);
        }
    }
    int num_levels() const { return num_levels_; }
    int pre_smooth_deg() const { return levels_.at(0).pre_smooth_deg; }
    int post_smooth_deg() const { return levels_.at(0).post_smooth_deg; }
    float fine_lambda_min() const { return levels_.at(0).lambda_min; }
    float fine_lambda_max() const { return levels_.at(0).lambda_max; }

  private:
    template <typename T, int HaloRadius>
    void compute_residual_vec(const T *d_rhs, const T *d_Au, T *d_res_out,
                              const MultigridLevelConfig &lvl,
                              hipStream_t stream) {
        dim3 threads(64, 4, 1);
        dim3 blocks(static_cast<unsigned int>((lvl.Nz + 63) / 64),
                    static_cast<unsigned int>((lvl.Ny + 3) / 4),
                    static_cast<unsigned int>(my_band_states_ * lvl.Nx));
        chebyshev_zero_guess_init_kernel<T, HaloRadius>
            <<<blocks, threads, 0, stream>>>(d_res_out, d_cheb_d_, d_rhs, 0.0f,
                                             lvl.Nx, lvl.Ny, lvl.Nz, lvl.Gx,
                                             lvl.Gy, lvl.Gz);
        chebyshev_recurrence_update_kernel<T, HaloRadius>
            <<<blocks, threads, 0, stream>>>(d_res_out, d_cheb_d_, d_rhs, d_Au,
                                             1.0f, 0.0f, lvl.Nx, lvl.Ny, lvl.Nz,
                                             lvl.Gx, lvl.Gy, lvl.Gz);
    }
    int my_band_states_ = 0;
    int num_levels_ = 1;
    int coarsest_solve_deg_ = 6;
    double helmholtz_shift_ = 0.01;
    FineStencil8th<WeightT> fine_stencil_{};
    std::vector<CoarseStencil2nd<WeightT>> coarse_stencils_;
    std::vector<MultigridLevelConfig> levels_;
    std::vector<LowT *> d_e_;
    std::vector<LowT *> d_rhs_;
    float *d_zero_vloc_ = nullptr;
    float *d_zero_eig_ = nullptr;
    LowT *d_cheb_d_ = nullptr;
    LowT *d_cheb_Au_ = nullptr;
    LowT *d_res_scratch_ = nullptr;
};
