#pragma once

#include "custom_stencils.hpp"
#include "chebyshev_rmmdiis.hpp"
#include <hip/hip_runtime.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

#ifndef HIP_DEVICE
#define HIP_DEVICE __attribute__((device))
#endif

HIP_DEVICE inline std::complex< double > rayleigh_to_complex64(float x) {
    return std::complex< double >(static_cast< double >(x), 0.0);
}
HIP_DEVICE inline std::complex< double > rayleigh_to_complex64(double x) {
    return std::complex< double >(x, 0.0);
}
HIP_DEVICE inline std::complex< double > rayleigh_to_complex64(std::complex< float > x) {
    return std::complex< double >(static_cast< double >(x.real()), static_cast< double >(x.imag()));
}
HIP_DEVICE inline std::complex< double > rayleigh_to_complex64(std::complex< double > x) {
    return x;
}

template < typename OutT >
HIP_DEVICE inline OutT rayleigh_cast_from_complex64(std::complex< double > v) {
    if constexpr (std::is_same_v< OutT, float >) {
        return static_cast< float >(v.real());
    } else if constexpr (std::is_same_v< OutT, double >) {
        return v.real();
    } else if constexpr (std::is_same_v< OutT, std::complex< float > >) {
        return std::complex< float >(static_cast< float >(v.real()), static_cast< float >(v.imag()));
    } else {
        return v;
    }
}

// Exact FP64 inner product Re<psi | r> using WaveT psi and FP32 residual r_32
template < typename WaveT, int HaloRadius >
HIP_GLOBAL void batched_psi_res_dot_kernel(
    const WaveT* HIP_RESTRICT psi,
    const typename PrecisionTraits< WaveT >::LowT* HIP_RESTRICT r_32,
    double* HIP_RESTRICT block_partials,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz) {
    int b   = static_cast< int >(blockIdx.y);
    int tid = static_cast< int >(threadIdx.x);

    size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
    size_t state_stride = static_cast< size_t >(Gx) * plane_stride;

    const WaveT* pb = psi  + static_cast< size_t >(b) * state_stride;
    const auto*  rb = r_32 + static_cast< size_t >(b) * state_stride;

    double local_dot = 0.0;
    int total_inner = Nx * Ny * Nz;

    for (int idx = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + tid;
         idx < total_inner;
         idx += static_cast< int >(gridDim.x) * static_cast< int >(blockDim.x)) {
        int iz  = idx % Nz;
        int tmp = idx / Nz;
        int iy  = tmp % Ny;
        int ix  = tmp / Ny;

        size_t g_idx = (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
                     +  static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
                     +  static_cast< size_t >(iz + HaloRadius);

        std::complex< double > pv = rayleigh_to_complex64(pb[g_idx]);
        std::complex< double > rv = rayleigh_to_complex64(rb[g_idx]);
        local_dot += pv.real() * rv.real() + pv.imag() * rv.imag();
    }

    __shared__ double sdata[256];
    sdata[tid] = local_dot;
    __syncthreads();

    for (int s = 128; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }

    if (tid == 0) {
        size_t out_idx = static_cast< size_t >(b) * static_cast< size_t >(gridDim.x) + static_cast< size_t >(blockIdx.x);
        block_partials[out_idx] = sdata[0];
    }
}

// Exact FP64 inner product Re<psi | V_nl * psi> to compute local-only shift eig_loc = eig - vel * Re<psi | V_nl * psi>
template < typename WaveT, int HaloRadius >
HIP_GLOBAL void batched_psi_vnl_dot_kernel(
    const WaveT* HIP_RESTRICT psi,
    const WaveT* HIP_RESTRICT vnl,
    double* HIP_RESTRICT block_partials,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz) {
    int b   = static_cast< int >(blockIdx.y);
    int tid = static_cast< int >(threadIdx.x);

    size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
    size_t state_stride = static_cast< size_t >(Gx) * plane_stride;

    const WaveT* pb  = psi + static_cast< size_t >(b) * state_stride;
    const WaveT* vnb = vnl + static_cast< size_t >(b) * state_stride;

    double local_dot = 0.0;
    int total_inner = Nx * Ny * Nz;

    for (int idx = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + tid;
         idx < total_inner;
         idx += static_cast< int >(gridDim.x) * static_cast< int >(blockDim.x)) {
        int iz  = idx % Nz;
        int tmp = idx / Nz;
        int iy  = tmp % Ny;
        int ix  = tmp / Ny;

        size_t g_idx = (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
                     +  static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
                     +  static_cast< size_t >(iz + HaloRadius);

        std::complex< double > pv  = rayleigh_to_complex64(pb[g_idx]);
        std::complex< double > vnv = rayleigh_to_complex64(vnb[g_idx]);
        local_dot += pv.real() * vnv.real() + pv.imag() * vnv.imag();
    }

    __shared__ double sdata[256];
    sdata[tid] = local_dot;
    __syncthreads();

    for (int s = 128; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }

    if (tid == 0) {
        size_t out_idx = static_cast< size_t >(b) * static_cast< size_t >(gridDim.x) + static_cast< size_t >(blockIdx.x);
        block_partials[out_idx] = sdata[0];
    }
}

// Computes local-only smoother eigenvalues eig_loc[b] = eig[b] - vel * Re<psi[b] | V_nl * psi[b]>
// so the local Chebyshev operator (T_kin + V_loc - eig_loc) has zero expectation on the target state
template < typename WaveT, int HaloRadius >
inline void launch_compute_local_smoother_shifts(
    const WaveT* d_psi,
    const WaveT* d_vnl,
    double* d_dot_partials,
    const std::vector< float >& h_eig_32,
    float* d_eig_loc_32,
    double vel,
    int my_band_states,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz,
    int num_blocks = 256,
    hipStream_t stream = nullptr) {
    if (my_band_states <= 0) return;

    dim3 threads(256, 1, 1);
    dim3 blocks(static_cast< unsigned int >(num_blocks), static_cast< unsigned int >(my_band_states), 1);

    batched_psi_vnl_dot_kernel< WaveT, HaloRadius ><<<blocks, threads, 0, stream>>>(
        d_psi, d_vnl, d_dot_partials, Nx, Ny, Nz, Gx, Gy, Gz);

    size_t total_partials = static_cast< size_t >(my_band_states) * static_cast< size_t >(num_blocks);
    std::vector< double > h_partials(total_partials);
    (void)hipMemcpyAsync(h_partials.data(), d_dot_partials,
                         total_partials * sizeof(double), hipMemcpyDeviceToHost, stream);
    (void)hipStreamSynchronize(stream);

    double scale_vol = (vel > 0.0) ? vel : 1.0;
    std::vector< float > h_eig_loc(static_cast< size_t >(my_band_states), 0.0f);

    for (int b = 0; b < my_band_states; ++b) {
        size_t ub = static_cast< size_t >(b);
        double vnl_dot = 0.0;
        size_t base = ub * static_cast< size_t >(num_blocks);
        for (int blk = 0; blk < num_blocks; ++blk) {
            vnl_dot += h_partials.at(base + static_cast< size_t >(blk));
        }
        double eig_loc = static_cast< double >(h_eig_32.at(ub)) - scale_vol * vnl_dot;
        h_eig_loc.at(ub) = static_cast< float >(eig_loc);
    }

    (void)hipMemcpyAsync(d_eig_loc_32, h_eig_loc.data(),
                         static_cast< size_t >(my_band_states) * sizeof(float), hipMemcpyHostToDevice, stream);
}

// Updates eigenvalues using exact full-Hamiltonian Rayleigh quotient:
//   eig_new = eig_old - beta_mix * clamp(vel * Re<psi | r>, -max_step, max_step)
template < typename WaveT, int HaloRadius >
inline void launch_exact_rayleigh_update(
    const WaveT* d_psi,
    const typename PrecisionTraits< WaveT >::LowT* d_r_32,
    double* d_dot_partials,
    std::vector< typename PrecisionTraits< WaveT >::RealT >& h_eig_wave,
    std::vector< float >& h_eig_32,
    typename PrecisionTraits< WaveT >::RealT* d_eig_wave,
    float* d_eig_32,
    double beta_mix,
    double max_step,
    double vel,
    int my_band_states,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz,
    int num_blocks = 256,
    hipStream_t stream = nullptr) {
    if (my_band_states <= 0 || beta_mix == 0.0) return;
    using RealT = typename PrecisionTraits< WaveT >::RealT;

    dim3 threads(256, 1, 1);
    dim3 blocks(static_cast< unsigned int >(num_blocks), static_cast< unsigned int >(my_band_states), 1);

    batched_psi_res_dot_kernel< WaveT, HaloRadius ><<<blocks, threads, 0, stream>>>(
        d_psi, d_r_32, d_dot_partials, Nx, Ny, Nz, Gx, Gy, Gz);

    size_t total_partials = static_cast< size_t >(my_band_states) * static_cast< size_t >(num_blocks);
    std::vector< double > h_partials(total_partials);
    (void)hipMemcpyAsync(h_partials.data(), d_dot_partials,
                         total_partials * sizeof(double), hipMemcpyDeviceToHost, stream);
    (void)hipStreamSynchronize(stream);

    double scale_vol = (vel > 0.0) ? vel : 1.0;

    for (int b = 0; b < my_band_states; ++b) {
        size_t ub = static_cast< size_t >(b);
        double dot_val = 0.0;
        size_t base = ub * static_cast< size_t >(num_blocks);
        for (int blk = 0; blk < num_blocks; ++blk) {
            dot_val += h_partials.at(base + static_cast< size_t >(blk));
        }
        double delta_eps = std::clamp(-scale_vol * dot_val, -max_step, max_step);
        double updated   = static_cast< double >(h_eig_wave.at(ub)) + beta_mix * delta_eps;
        h_eig_wave.at(ub) = static_cast< RealT >(updated);
        h_eig_32.at(ub)   = static_cast< float >(updated);
    }
    static bool printed_initial_eigs = false;
#if 0
    if (!printed_initial_eigs) {
        printed_initial_eigs = true;
        int world_rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        for (int b = 0; b < my_band_states; ++b) {
            size_t ub = static_cast< size_t >(b);
            double dot_val = 0.0;
            size_t base = ub * static_cast< size_t >(num_blocks);
            for (int blk = 0; blk < num_blocks; ++blk) {
                dot_val += h_partials.at(base + static_cast< size_t >(blk));
            }
            double eig_in     = static_cast< double >(h_eig_wave.at(ub));
            double rayleigh_H = eig_in - scale_vol * dot_val;
            printf("[GPU-EIG-FIRST] rank=%d b=%d: eig_in=%.10f, rayleigh_H=%.10f, diff=%.3e\n",
                         world_rank, b, eig_in, rayleigh_H, rayleigh_H - eig_in);
        }
        fflush(NULL);
    }
#endif

    (void)hipMemcpyAsync(d_eig_32, h_eig_32.data(),
                         static_cast< size_t >(my_band_states) * sizeof(float), hipMemcpyHostToDevice, stream);
    if (d_eig_wave != nullptr && reinterpret_cast< const void* >(d_eig_wave) != reinterpret_cast< const void* >(d_eig_32)) {
        (void)hipMemcpyAsync(d_eig_wave, h_eig_wave.data(),
static_cast< size_t >(my_band_states) * sizeof(RealT), hipMemcpyHostToDevice, stream);
}
}
template < typename FieldT, int HaloRadius >
HIP_GLOBAL void copy_interior_field_kernel(
const FieldT* HIP_RESTRICT src,
FieldT* HIP_RESTRICT dst,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz) {
int iz = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + static_cast< int >(threadIdx.x);
int iy = static_cast< int >(blockIdx.y) * static_cast< int >(blockDim.y) + static_cast< int >(threadIdx.y);
int ix = static_cast< int >(blockIdx.z) % Nx;
int b = static_cast< int >(blockIdx.z) / Nx;
if (iz >= Nz || iy >= Ny) return;
size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
size_t state_stride = static_cast< size_t >(Gx) * plane_stride;
size_t idx = static_cast< size_t >(b) * state_stride
+ (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
+ static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
+ static_cast< size_t >(iz + HaloRadius);
dst[idx] = src[idx];
}
template < typename FieldT, int HaloRadius >
HIP_GLOBAL void combine_diis_history_kernel(
const FieldT* const* HIP_RESTRICT hist_ptrs,
const typename PrecisionTraits< FieldT >::LowT* HIP_RESTRICT coeffs,
FieldT* HIP_RESTRICT dst,
int hist_len,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz) {
int iz = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + static_cast< int >(threadIdx.x);
int iy = static_cast< int >(blockIdx.y) * static_cast< int >(blockDim.y) + static_cast< int >(threadIdx.y);
int ix = static_cast< int >(blockIdx.z) % Nx;
int b = static_cast< int >(blockIdx.z) / Nx;
if (iz >= Nz || iy >= Ny) return;
size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
size_t state_stride = static_cast< size_t >(Gx) * plane_stride;
size_t idx = static_cast< size_t >(b) * state_stride
+ (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
+ static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
+ static_cast< size_t >(iz + HaloRadius);
std::complex< double > acc{0.0, 0.0};
for (int m = 0; m < hist_len; ++m) {
std::complex< double > cm = rayleigh_to_complex64(coeffs[static_cast< size_t >(b * hist_len + m)]);
std::complex< double > vm = rayleigh_to_complex64(hist_ptrs[m][idx]);
acc += cm * vm;
}
dst[idx] = rayleigh_cast_from_complex64< FieldT >(acc);
}
template < typename LowT, int HaloRadius >
inline std::vector< std::vector< std::complex< double > > > compute_diis_gram_matrix(
const std::array< LowT*, 3 >& d_hist_r,
int hist_len,
std::complex< double >* d_dot_partials,
int my_band_states,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz,
int num_blocks = 256,
hipStream_t stream = nullptr) {
std::vector< std::vector< std::complex< double > > > grams(
static_cast< size_t >(my_band_states),
std::vector< std::complex< double > >(static_cast< size_t >(hist_len * hist_len), std::complex< double >(0.0, 0.0)));
if (my_band_states <= 0 || hist_len <= 1) return grams;
auto full_m3 = compute_batched_gram_m3< LowT, HaloRadius >(
d_hist_r.at(0),
d_hist_r.at(hist_len > 1 ? 1 : 0),
d_hist_r.at(hist_len > 2 ? 2 : 0),
d_dot_partials, my_band_states, Nx, Ny, Nz, Gx, Gy, Gz, num_blocks, stream);
for (int b = 0; b < my_band_states; ++b) {
size_t ub = static_cast< size_t >(b);
for (int i = 0; i < hist_len; ++i) {
for (int j = 0; j < hist_len; ++j) {
grams.at(ub).at(static_cast< size_t >(i * hist_len + j)) =
full_m3.at(ub).at(static_cast< size_t >(i * 3 + j));
}
}
}
return grams;
}

