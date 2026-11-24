#pragma once

#include "chebyshev_rmmdiis.hpp"
#include "custom_stencils.hpp"
#include "periodic_halo.hpp"
#include "rayleigh_diis_helpers.hpp"

#include <cmath>
#include <complex>
#include <hip/hip_runtime.h>
#include <vector>

template <typename WaveT>
inline void launch_eval_full_residual_wave(
    const WaveT *d_u, const WaveT *d_vnl,
    const typename PrecisionTraits<WaveT>::RealT *d_vloc_wave,
    const typename PrecisionTraits<WaveT>::RealT *d_eig_wave, WaveT *d_r_wave,
    typename PrecisionTraits<WaveT>::LowT *d_r_32, int my_states, int Nx,
    int Ny, int Nz, int Gx, int Gy, int Gz,
    const std::vector<StencilTap<std::complex<double>>> &fine_taps,
    hipStream_t stream = nullptr) {
    launch_fine_residual_and_downcast<WaveT>(
        d_u, d_vnl, d_vloc_wave, d_eig_wave, d_r_32, my_states, Nx, Ny, Nz, Gx,
        Gy, Gz, fine_taps, stream, d_r_wave);
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void mgsmoother_init_p_wave_kernel(
    WaveT *HIP_RESTRICT p_wave, const WaveT *HIP_RESTRICT r_wave,
    const typename PrecisionTraits<WaveT>::RealT *HIP_RESTRICT vloc,
    double w0_diag, int Nx, int Ny, int Nz, int Gx, int Gy, int Gz) {
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

    size_t grid_idx =
        (static_cast<size_t>(ix + HaloRadius) * static_cast<size_t>(Gy) +
         static_cast<size_t>(iy + HaloRadius)) *
            static_cast<size_t>(Gz) +
        static_cast<size_t>(iz + HaloRadius);
    size_t idx = static_cast<size_t>(b) * state_stride + grid_idx;

    double denom = fabs(w0_diag - static_cast<double>(vloc[grid_idx]));
    RealT dinv = (denom > 1e-12) ? static_cast<RealT>(1.0 / denom) : RealT{1};
    p_wave[idx] = r_wave[idx] * dinv;
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void mgsmoother_update_p_wave_kernel(
    WaveT *HIP_RESTRICT p_wave, const WaveT *HIP_RESTRICT r_wave,
    const typename PrecisionTraits<WaveT>::RealT *HIP_RESTRICT vloc,
    double w0_diag, double b_coeff, int Nx, int Ny, int Nz, int Gx, int Gy,
    int Gz) {
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

    size_t grid_idx =
        (static_cast<size_t>(ix + HaloRadius) * static_cast<size_t>(Gy) +
         static_cast<size_t>(iy + HaloRadius)) *
            static_cast<size_t>(Gz) +
        static_cast<size_t>(iz + HaloRadius);
    size_t idx = static_cast<size_t>(b) * state_stride + grid_idx;

    double denom = fabs(w0_diag - static_cast<double>(vloc[grid_idx]));
    RealT dinv = (denom > 1e-12) ? static_cast<RealT>(1.0 / denom) : RealT{1};
    RealT bc = static_cast<RealT>(b_coeff);
    p_wave[idx] = r_wave[idx] * dinv + p_wave[idx] * bc;
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void mgsmoother_axpy_u_wave_kernel(WaveT *HIP_RESTRICT u_wave,
                                              const WaveT *HIP_RESTRICT p_wave,
                                              double a_coeff, int Nx, int Ny,
                                              int Nz, int Gx, int Gy, int Gz) {
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

    u_wave[idx] += static_cast<RealT>(a_coeff) * p_wave[idx];
}

template <typename WaveT>
inline void launch_rmg_mgsmoother_gpu(
    WaveT *d_u, const WaveT *d_vnl,
    const typename PrecisionTraits<WaveT>::RealT *d_vloc_wave,
    typename PrecisionTraits<WaveT>::RealT *d_eig_wave, float *d_eig_32,
    std::vector<typename PrecisionTraits<WaveT>::RealT> &h_eig_wave,
    std::vector<float> &h_eig_32, WaveT *d_r_wave, WaveT *d_p_wave,
    typename PrecisionTraits<WaveT>::LowT *d_r_32, double *d_norm_partials,
    typename PrecisionTraits<WaveT>::RealT *d_inv_norms,
    const std::vector<StencilTap<std::complex<double>>> &fine_taps, float lmin,
    float lmax, int order, int my_states, int Nx, int Ny, int Nz, int Gx,
    int Gy, int Gz, double vel, bool periodic, double eig_mix_fraction = 0.7,
    double eig_max_step = 0.25, hipStream_t stream = nullptr,
    bool print_diag = false) {
    (void)print_diag;
    if (my_states <= 0)
        return;
    constexpr int DOT_BLOCKS = 256;

    const double theta =
        0.5 * (static_cast<double>(lmax) + static_cast<double>(lmin));
    const double delta =
        0.5 * (static_cast<double>(lmax) - static_cast<double>(lmin));
    const double sigma = (std::abs(delta) > 1e-12) ? (theta / delta) : 1.0;

    //    double w0_diag = fine_taps.empty() ? 1.0 :
    //    fine_taps.at(0).weight.real();
    double w0_diag = 0.0;
    for (const auto &tap : fine_taps) {
        if (tap.dx == 0 && tap.dy == 0 && tap.dz == 0) {
            w0_diag += tap.weight.real();
        }
    }
    if (std::abs(w0_diag) < 1e-12)
        w0_diag = 1.0;

    dim3 threads(64, 4, 1);
    dim3 blocks(static_cast<unsigned int>((Nz + 63) / 64),
                static_cast<unsigned int>((Ny + 3) / 4),
                static_cast<unsigned int>(my_states * Nx));

    launch_fill_periodic_halo<WaveT, 4>(d_u, my_states, Nx, Ny, Nz, Gx, Gy, Gz,
                                        periodic, stream);

    launch_eval_full_residual_wave<WaveT>(d_u, d_vnl, d_vloc_wave, d_eig_wave,
                                          d_r_wave, d_r_32, my_states, Nx, Ny,
                                          Nz, Gx, Gy, Gz, fine_taps, stream);

    launch_exact_rayleigh_update<WaveT, 4>(
        d_u, d_r_32, d_norm_partials, h_eig_wave, h_eig_32, d_eig_wave,
        d_eig_32, eig_mix_fraction, eig_max_step, vel, my_states, Nx, Ny, Nz,
        Gx, Gy, Gz, DOT_BLOCKS, stream);

    launch_renormalize_wavefunctions<WaveT, 4>(
        d_u, d_norm_partials, d_inv_norms, vel, my_states, Nx, Ny, Nz, Gx, Gy,
        Gz, DOT_BLOCKS, stream);

    launch_fill_periodic_halo<WaveT, 4>(d_u, my_states, Nx, Ny, Nz, Gx, Gy, Gz,
                                        periodic, stream);

    launch_eval_full_residual_wave<WaveT>(d_u, d_vnl, d_vloc_wave, d_eig_wave,
                                          d_r_wave, d_r_32, my_states, Nx, Ny,
                                          Nz, Gx, Gy, Gz, fine_taps, stream);

    mgsmoother_init_p_wave_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
        d_p_wave, d_r_wave, d_vloc_wave, w0_diag, Nx, Ny, Nz, Gx, Gy, Gz);

    double a = 0.0;
    double b = 0.0;

    for (int k = 0; k < order; ++k) {
        if (k == 0) {
            a = 2.0 / theta;
            b = 0.0;
        } else if (k == 1) {
            a = 2.0 / theta;
            b = (a * delta * delta) / 4.0;
        } else {
            const double beta_new = 1.0 / (2.0 * sigma - b);
            b = beta_new;
            a = 2.0 * beta_new / delta;
        }

        mgsmoother_axpy_u_wave_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
            d_u, d_p_wave, a, Nx, Ny, Nz, Gx, Gy, Gz);

        launch_renormalize_wavefunctions<WaveT, 4>(
            d_u, d_norm_partials, d_inv_norms, vel, my_states, Nx, Ny, Nz, Gx,
            Gy, Gz, DOT_BLOCKS, stream);

        launch_fill_periodic_halo<WaveT, 4>(d_u, my_states, Nx, Ny, Nz, Gx, Gy,
                                            Gz, periodic, stream);

        launch_eval_full_residual_wave<WaveT>(
            d_u, d_vnl, d_vloc_wave, d_eig_wave, d_r_wave, d_r_32, my_states,
            Nx, Ny, Nz, Gx, Gy, Gz, fine_taps, stream);

        launch_exact_rayleigh_update<WaveT, 4>(
            d_u, d_r_32, d_norm_partials, h_eig_wave, h_eig_32, d_eig_wave,
            d_eig_32, eig_mix_fraction, eig_max_step, vel, my_states, Nx, Ny,
            Nz, Gx, Gy, Gz, DOT_BLOCKS, stream);

        if (k + 1 < order) {
            launch_eval_full_residual_wave<WaveT>(
                d_u, d_vnl, d_vloc_wave, d_eig_wave, d_r_wave, d_r_32,
                my_states, Nx, Ny, Nz, Gx, Gy, Gz, fine_taps, stream);

            mgsmoother_update_p_wave_kernel<WaveT, 4>
                <<<blocks, threads, 0, stream>>>(d_p_wave, d_r_wave,
                                                 d_vloc_wave, w0_diag, b, Nx,
                                                 Ny, Nz, Gx, Gy, Gz);
        }
    }
}
