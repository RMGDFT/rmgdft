#pragma once

#include "chebyshev_rmmdiis.hpp"
#include "custom_stencils.hpp"
#include "mg_smoother_gpu.hpp"
#include "multigrid_preconditioner.hpp"
#include "periodic_halo.hpp"
#include "rayleigh_diis_helpers.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <hip/hip_runtime.h>
#include <vector>

template <typename WaveT>
HIP_GLOBAL void downcast_field_to_fp32_kernel(
    const WaveT *HIP_RESTRICT src,
    typename PrecisionTraits<WaveT>::LowT *HIP_RESTRICT dst, int Nx, int Ny,
    int Nz, int Gx, int Gy, int Gz) {
    constexpr int HALO = 4;
    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int ix = static_cast<int>(blockIdx.z) % Nx;
    int b = static_cast<int>(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny)
        return;

    size_t plane_stride = static_cast<size_t>(Gy) * static_cast<size_t>(Gz);
    size_t state_stride = static_cast<size_t>(Gx) * plane_stride;

    size_t idx = static_cast<size_t>(b) * state_stride +
                 (static_cast<size_t>(ix + HALO) * static_cast<size_t>(Gy) +
                  static_cast<size_t>(iy + HALO)) *
                     static_cast<size_t>(Gz) +
                 static_cast<size_t>(iz + HALO);

    dst[idx] = PrecisionTraits<WaveT>::to_low(src[idx]);
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void project_correction_orthogonal_kernel(
    const WaveT *HIP_RESTRICT psi,
    typename PrecisionTraits<WaveT>::LowT *HIP_RESTRICT corr_32,
    const std::complex<double> *HIP_RESTRICT overlaps, double vel, int Nx,
    int Ny, int Nz, int Gx, int Gy, int Gz) {
    using LowT = typename PrecisionTraits<WaveT>::LowT;

    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int ix = static_cast<int>(blockIdx.z) % Nx;
    int b = static_cast<int>(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny)
        return;

    size_t plane_stride = static_cast<size_t>(Gy) * static_cast<size_t>(Gz);
    size_t state_stride = static_cast<size_t>(Gx) * plane_stride;

    size_t idx =
        static_cast<size_t>(b) * state_stride +
        (static_cast<size_t>(ix + HaloRadius) * static_cast<size_t>(Gy) +
         static_cast<size_t>(iy + HaloRadius)) *
            static_cast<size_t>(Gz) +
        static_cast<size_t>(iz + HaloRadius);

    double scale_vol = (vel > 0.0) ? vel : 1.0;
    std::complex<double> ov = overlaps[b] * scale_vol;
    std::complex<double> pv = rayleigh_to_complex64(psi[idx]);
    std::complex<double> cv = rayleigh_to_complex64(corr_32[idx]);

    corr_32[idx] = rayleigh_cast_from_complex64<LowT>(cv - ov * pv);
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void build_rmmres_signed_kernel(
    const WaveT *HIP_RESTRICT saved_psi, const WaveT *HIP_RESTRICT smoothed_psi,
    const typename PrecisionTraits<WaveT>::LowT *HIP_RESTRICT mg_corr_32,
    bool has_mg_corr, double mg_weight, WaveT *HIP_RESTRICT rmmres_wave,
    typename PrecisionTraits<WaveT>::LowT *HIP_RESTRICT rmmres_32, int Nx,
    int Ny, int Nz, int Gx, int Gy, int Gz) {
    using RealT = typename PrecisionTraits<WaveT>::RealT;

    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int ix = static_cast<int>(blockIdx.z) % Nx;
    int b = static_cast<int>(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny)
        return;

    size_t plane_stride = static_cast<size_t>(Gy) * static_cast<size_t>(Gz);
    size_t state_stride = static_cast<size_t>(Gx) * plane_stride;

    size_t idx =
        static_cast<size_t>(b) * state_stride +
        (static_cast<size_t>(ix + HaloRadius) * static_cast<size_t>(Gy) +
         static_cast<size_t>(iy + HaloRadius)) *
            static_cast<size_t>(Gz) +
        static_cast<size_t>(iz + HaloRadius);

    WaveT diff = saved_psi[idx] - smoothed_psi[idx];
    if (has_mg_corr && mg_corr_32 != nullptr) {
        diff -= static_cast<RealT>(mg_weight) *
                cast_from_low<WaveT>(mg_corr_32[idx]);
    }
    rmmres_wave[idx] = diff;
    rmmres_32[idx] = PrecisionTraits<WaveT>::to_low(diff);
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void rmm_step_axpy_per_band_kernel(
    const WaveT *HIP_RESTRICT base_psi, const WaveT *HIP_RESTRICT rmmres_wave,
    WaveT *HIP_RESTRICT psi_out,
    const typename PrecisionTraits<WaveT>::RealT *HIP_RESTRICT lambdas, int Nx,
    int Ny, int Nz, int Gx, int Gy, int Gz) {
    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int ix = static_cast<int>(blockIdx.z) % Nx;
    int b = static_cast<int>(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny)
        return;

    size_t plane_stride = static_cast<size_t>(Gy) * static_cast<size_t>(Gz);
    size_t state_stride = static_cast<size_t>(Gx) * plane_stride;

    size_t idx =
        static_cast<size_t>(b) * state_stride +
        (static_cast<size_t>(ix + HaloRadius) * static_cast<size_t>(Gy) +
         static_cast<size_t>(iy + HaloRadius)) *
            static_cast<size_t>(Gz) +
        static_cast<size_t>(iz + HaloRadius);

    psi_out[idx] = base_psi[idx] + lambdas[b] * rmmres_wave[idx];
}

template <typename WaveT>
void run_single_defect_correction_gpu(
    WaveT *d_psi_band, const WaveT *d_vnl_band,
    const typename PrecisionTraits<WaveT>::RealT *d_vloc_wave,
    typename PrecisionTraits<WaveT>::RealT *d_eig_wave, float *d_eig_32,
    std::vector<typename PrecisionTraits<WaveT>::RealT> &h_eig_wave,
    std::vector<float> &h_eig_32, WaveT *d_r_wave, WaveT *d_p_wave,
    typename PrecisionTraits<WaveT>::LowT *d_r_32,
    typename PrecisionTraits<WaveT>::LowT *d_e_32,
    std::array<WaveT *, 3> &d_hist_psi,
    std::array<typename PrecisionTraits<WaveT>::LowT *, 3> &d_hist_r,
    WaveT *const *d_hist_psi_ptrs,
    typename PrecisionTraits<WaveT>::LowT *const *d_hist_r_ptrs,
    typename PrecisionTraits<WaveT>::LowT *d_diis_coeffs,
    std::complex<double> *d_dot_partials, double *d_norm_partials,
    typename PrecisionTraits<WaveT>::RealT *d_inv_norms,
    MultigridPreconditioner<typename PrecisionTraits<WaveT>::LowT> &mg,
    const std::vector<StencilTap<std::complex<double>>> &fine_taps,
    int diis_hist_len, int my_states, int Nx, int Ny, int Nz, int Gx, int Gy,
    int Gz, double vel = 1.0, bool periodic = true,
    double eig_mix_fraction = 0.7, double eig_max_step = 1.0e6,
    double fg_step = 1.0, hipStream_t stream = nullptr) {
    if (my_states <= 0)
        return;
    using LowT = typename PrecisionTraits<WaveT>::LowT;
    using RealT = typename PrecisionTraits<WaveT>::RealT;
    constexpr int DOT_BLOCKS = 256;

    dim3 threads(64, 4, 1);
    dim3 blocks(static_cast<unsigned int>((Nz + 63) / 64),
                static_cast<unsigned int>((Ny + 3) / 4),
                static_cast<unsigned int>(my_states * Nx));

    int slot = diis_hist_len - 1;
    WaveT *d_saved_psi = d_hist_psi.at(static_cast<size_t>(slot));

    // 1. Save input wavefunction saved_psi into RMM-DIIS history slot
    copy_interior_field_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
        d_psi_band, d_saved_psi, Nx, Ny, Nz, Gx, Gy, Gz);

    // 2. Fine-grid Chebyshev smoother (RMGDFT mg_smoother.cpp, order = 3)
    int smooth_order = (mg.pre_smooth_deg() > 0) ? mg.pre_smooth_deg() : 3;
    launch_rmg_mgsmoother_gpu<WaveT>(
        d_psi_band, d_vnl_band, d_vloc_wave, d_eig_wave, d_eig_32, h_eig_wave,
        h_eig_32, d_r_wave, d_p_wave, d_r_32, d_norm_partials, d_inv_norms,
        fine_taps, mg.fine_lambda_min(), mg.fine_lambda_max(), smooth_order,
        my_states, Nx, Ny, Nz, Gx, Gy, Gz, vel, periodic, eig_mix_fraction,
        eig_max_step, stream, false);

    // 3. Coarse-grid multigrid residual correction (active when mg_num_levels >
    // 1)
    bool has_mg = (mg.num_levels() > 1);
    double mg_weight = 0.25;
    if (has_mg) {
        launch_fill_periodic_halo<WaveT, 4>(d_psi_band, my_states, Nx, Ny, Nz,
                                            Gx, Gy, Gz, periodic, stream);
        launch_eval_full_residual_wave<WaveT>(
            d_psi_band, d_vnl_band, d_vloc_wave, d_eig_wave, d_r_wave, d_r_32,
            my_states, Nx, Ny, Nz, Gx, Gy, Gz, fine_taps, stream);

        mg.coarse_helmholtz_correction(d_e_32, d_r_32, periodic, stream);

        downcast_field_to_fp32_kernel<WaveT><<<blocks, threads, 0, stream>>>(
            d_psi_band, d_r_32, Nx, Ny, Nz, Gx, Gy, Gz);
        auto ov_gram = compute_batched_gram_m3<LowT, 4>(
            d_r_32, d_e_32, d_e_32, d_dot_partials, my_states, Nx, Ny, Nz, Gx,
            Gy, Gz, DOT_BLOCKS, stream);

        std::vector<std::complex<double>> h_overlaps(
            static_cast<size_t>(my_states));
        for (int b = 0; b < my_states; ++b) {
            h_overlaps.at(static_cast<size_t>(b)) =
                ov_gram.at(static_cast<size_t>(b)).at(1);
        }
        (void)hipMemcpyAsync(d_dot_partials, h_overlaps.data(),
                             static_cast<size_t>(my_states) *
                                 sizeof(std::complex<double>),
                             hipMemcpyHostToDevice, stream);

        project_correction_orthogonal_kernel<WaveT, 4>
            <<<blocks, threads, 0, stream>>>(d_psi_band, d_e_32, d_dot_partials,
                                             vel, Nx, Ny, Nz, Gx, Gy, Gz);
    }

    // 4. Construct total preconditioned residual rmmres = (saved_psi -
    // smoothed_psi) - mg_weight * mg_corr
    LowT *d_rmmres_slot_32 = d_hist_r.at(static_cast<size_t>(slot));
    build_rmmres_signed_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
        d_saved_psi, d_psi_band, d_e_32, has_mg, mg_weight, d_p_wave,
        d_rmmres_slot_32, Nx, Ny, Nz, Gx, Gy, Gz);
    // 5. RMM-DIIS combination + line-minimized lambda step (matches RmmDiis.cpp
    // compute_lambda)
    if (diis_hist_len >= 2) {
        auto grams = compute_diis_gram_matrix<LowT, 4>(
            d_hist_r, diis_hist_len, d_dot_partials, my_states, Nx, Ny, Nz, Gx,
            Gy, Gz, DOT_BLOCKS, stream);
        std::vector<LowT> h_coeffs(
            static_cast<size_t>(my_states * diis_hist_len));
        std::vector<std::complex<double>> c_vec;
        for (int b = 0; b < my_states; ++b) {
            solve_small_diis_system(diis_hist_len,
                                    grams.at(static_cast<size_t>(b)), c_vec);
            for (int m = 0; m < diis_hist_len; ++m) {
                h_coeffs.at(static_cast<size_t>(b * diis_hist_len + m)) =
                    cast_scalar_from_complex64<LowT>(
                        c_vec.at(static_cast<size_t>(m)));
            }
        }
        (void)hipMemcpyAsync(d_diis_coeffs, h_coeffs.data(),
                             static_cast<size_t>(my_states * diis_hist_len) *
                                 sizeof(LowT),
                             hipMemcpyHostToDevice, stream);
        // Line minimization for lambda on vcycle >= 2 (diis_hist_len >= 3,
        // matches MgEigState.cpp line 82)
        std::vector<RealT> h_lambdas(static_cast<size_t>(my_states),
                                     static_cast<RealT>(-fg_step));
        if (diis_hist_len >= 2) {
            // d_r_32 = r_init = eig * u0 - H * u0 (using d_saved_psi,
            // d_vnl_band, and d_r_wave scratch)
            launch_fill_periodic_halo<WaveT, 4>(d_saved_psi, my_states, Nx, Ny,
                                                Nz, Gx, Gy, Gz, periodic,
                                                stream);
            launch_eval_full_residual_wave<WaveT>(
                d_saved_psi, d_vnl_band, d_vloc_wave, d_eig_wave, d_r_wave,
                d_r_32, my_states, Nx, Ny, Nz, Gx, Gy, Gz, fine_taps, stream);
            // Zero d_r_wave to act as zero V_nl for H_loc * rmmres, then
            // evaluate r_pr0 into d_e_32
            size_t band_bytes = static_cast<size_t>(my_states) *
                                static_cast<size_t>(Gx) *
                                static_cast<size_t>(Gy) *
                                static_cast<size_t>(Gz) * sizeof(WaveT);
            (void)hipMemsetAsync(d_r_wave, 0, band_bytes, stream);
            launch_fill_periodic_halo<WaveT, 4>(d_p_wave, my_states, Nx, Ny, Nz,
                                                Gx, Gy, Gz, periodic, stream);
            launch_eval_full_residual_wave<WaveT>(
                d_p_wave, d_r_wave, d_vloc_wave, d_eig_wave, d_r_wave, d_e_32,
                my_states, Nx, Ny, Nz, Gx, Gy, Gz, fine_taps, stream);
            auto lam_gram = compute_batched_gram_m3<LowT, 4>(
                d_e_32, d_r_32, d_e_32, d_dot_partials, my_states, Nx, Ny, Nz,
                Gx, Gy, Gz, DOT_BLOCKS, stream);
            for (int b = 0; b < my_states; ++b) {
                const auto &G = lam_gram.at(static_cast<size_t>(b));
                double denom = G.at(0).real(); // <r_pr0 | r_pr0>
                double numer = G.at(1).real(); // Re<r_pr0 | r_init>
                if (denom > 1e-24) {
                    double step_mag = std::clamp(numer / denom, 0.25 * fg_step,
                                                 2.5 * fg_step);
                    h_lambdas.at(static_cast<size_t>(b)) =
                        static_cast<RealT>(-step_mag);
                }
            }
        }
        // Compute optimal DIIS linear combination d_r_wave = sum_m c_m *
        // saved_psi_m
        combine_diis_history_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
            d_hist_psi_ptrs, d_diis_coeffs, d_r_wave, diis_hist_len, Nx, Ny, Nz,
            Gx, Gy, Gz);
        // Reuse d_inv_norms (already allocated to safe_states * sizeof(RealT))
        // to hold per-band lambdas
        (void)hipMemcpyAsync(d_inv_norms, h_lambdas.data(),
                             static_cast<size_t>(my_states) * sizeof(RealT),
                             hipMemcpyHostToDevice, stream);
        rmm_step_axpy_per_band_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
            d_r_wave, d_p_wave, d_psi_band, d_inv_norms, Nx, Ny, Nz, Gx, Gy,
            Gz);
    } else if (has_mg) {
        launch_apply_fp32_correction<WaveT, 4>(
            d_psi_band, d_e_32, my_states, Nx, Ny, Nz, Gx, Gy, Gz, stream,
            static_cast<float>(mg_weight * fg_step));
    }
    launch_renormalize_wavefunctions<WaveT, 4>(
        d_psi_band, d_norm_partials, d_inv_norms, vel, my_states, Nx, Ny, Nz,
        Gx, Gy, Gz, DOT_BLOCKS, stream);
}
