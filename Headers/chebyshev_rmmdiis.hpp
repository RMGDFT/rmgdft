#pragma once

#include "periodic_halo.hpp"
#include "custom_stencils.hpp"
#include <hip/hip_runtime.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#ifndef HIP_DEVICE
#define HIP_DEVICE __attribute__((device))
#endif

#ifndef HIP_GLOBAL
#define HIP_GLOBAL __attribute__((global))
#endif

#ifndef HIP_RESTRICT
#define HIP_RESTRICT __restrict__
#endif

constexpr int DOT_BLOCK_SIZE = 256;

HIP_DEVICE inline std::complex< double > to_complex64(float x) {
    return std::complex< double >(static_cast< double >(x), 0.0);
}
HIP_DEVICE inline std::complex< double > to_complex64(double x) {
    return std::complex< double >(x, 0.0);
}
HIP_DEVICE inline std::complex< double > to_complex64(std::complex< float > x) {
    return std::complex< double >(static_cast< double >(x.real()), static_cast< double >(x.imag()));
}
HIP_DEVICE inline std::complex< double > to_complex64(std::complex< double > x) {
    return x;
}

template < typename OutT >
__host__ __device__ inline OutT cast_scalar_from_complex64(std::complex< double > v) {
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

template < typename WaveT >
HIP_DEVICE inline WaveT cast_from_low(typename PrecisionTraits< WaveT >::LowT v) {
    if constexpr (std::is_same_v< WaveT, double >) {
        return static_cast< double >(v);
    } else if constexpr (std::is_same_v< WaveT, std::complex< double > >) {
        return std::complex< double >(static_cast< double >(v.real()), static_cast< double >(v.imag()));
    } else {
        return v;
    }
}

// ============================================================================
// 1. CHEBYSHEV SMOOTHER KERNELS
// ============================================================================

template < typename T, int HaloRadius >
HIP_GLOBAL void chebyshev_zero_guess_init_kernel(
    T* HIP_RESTRICT u,
    T* HIP_RESTRICT d,
    const T* HIP_RESTRICT rhs,
    float inv_theta,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz) {
    int iz = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + static_cast< int >(threadIdx.x);
    int iy = static_cast< int >(blockIdx.y) * static_cast< int >(blockDim.y) + static_cast< int >(threadIdx.y);
    int ix = static_cast< int >(blockIdx.z) % Nx;
    int b  = static_cast< int >(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny) return;

    size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
    size_t state_stride = static_cast< size_t >(Gx) * plane_stride;

    size_t idx = static_cast< size_t >(b) * state_stride
               + (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
               +  static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
               +  static_cast< size_t >(iz + HaloRadius);

    T val = rhs[idx] * inv_theta;
    u[idx] = val;
    d[idx] = val;
}

template < typename T, int HaloRadius >
HIP_GLOBAL void chebyshev_recurrence_update_kernel(
    T* HIP_RESTRICT u,
    T* HIP_RESTRICT d,
    const T* HIP_RESTRICT rhs,
    const T* HIP_RESTRICT Au,
    float alpha,
    float beta,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz) {
    int iz = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + static_cast< int >(threadIdx.x);
    int iy = static_cast< int >(blockIdx.y) * static_cast< int >(blockDim.y) + static_cast< int >(threadIdx.y);
    int ix = static_cast< int >(blockIdx.z) % Nx;
    int b  = static_cast< int >(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny) return;

    size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
    size_t state_stride = static_cast< size_t >(Gx) * plane_stride;

    size_t idx = static_cast< size_t >(b) * state_stride
               + (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
               +  static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
               +  static_cast< size_t >(iz + HaloRadius);

    T r_new = rhs[idx] - Au[idx];
    T d_new = d[idx] * beta + r_new * alpha;
    d[idx]  = d_new;
    u[idx]  = u[idx] + d_new;
}

template < typename LowT, int HaloRadius, int MaxOffTaps >
inline void chebyshev_smooth_batched(
    LowT* d_u,
    LowT* d_d,
    LowT* d_Au,
    const LowT* d_rhs,
    const float* d_vloc,
    const float* d_eigvals,
    int my_band_states,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz,
    const CustomStencil< typename PrecisionTraits< LowT >::WeightT, MaxOffTaps >& stencil,
    float lambda_min,
    float lambda_max,
    int degree,
    bool zero_initial_guess = true,
    bool periodic = true,
    hipStream_t stream = nullptr) {
    if (my_band_states <= 0 || degree <= 0) return;

    float theta = 0.5f * (lambda_max + lambda_min);
    float delta = 0.5f * (lambda_max - lambda_min);
    if (theta <= 1e-12f) theta = 1.0f;

    float inv_theta = 1.0f / theta;
    float sigma = (delta > 1e-12f) ? (theta / delta) : 1.0f;
    float rho_prev = 1.0f / sigma;

    dim3 threads(64, 4, 1);
    dim3 blocks(static_cast< unsigned int >((Nz + 63) / 64),
                static_cast< unsigned int >((Ny + 3) / 4),
                static_cast< unsigned int >(my_band_states * Nx));

    int start_k = 0;
    if (zero_initial_guess) {
        chebyshev_zero_guess_init_kernel< LowT, HaloRadius ><<<blocks, threads, 0, stream>>>(
            d_u, d_d, d_rhs, inv_theta, Nx, Ny, Nz, Gx, Gy, Gz);
        start_k = 1;
    }

    for (int k = start_k; k < degree; ++k) {
        launch_shifted_hamiltonian< LowT, HaloRadius, MaxOffTaps >(
            d_u, d_Au, d_vloc, d_eigvals, my_band_states, Nx, Ny, Nz, Gx, Gy, Gz, stencil, periodic, stream);

        float alpha = 0.0f;
        float beta  = 0.0f;
        if (k == 0) {
            alpha = inv_theta;
            beta  = 0.0f;
        } else {
            float rho_curr = 1.0f / (2.0f * sigma - rho_prev);
            alpha = 2.0f * rho_curr / delta;
            beta  = rho_curr * rho_prev;
            rho_prev = rho_curr;
        }

        chebyshev_recurrence_update_kernel< LowT, HaloRadius ><<<blocks, threads, 0, stream>>>(
            d_u, d_d, d_rhs, d_Au, alpha, beta, Nx, Ny, Nz, Gx, Gy, Gz);
    }
}

// ============================================================================
// 2. BATCHED GRAM MATRIX & SMALL DIIS LINEAR SOLVER
// ============================================================================

template < typename LowT, int HaloRadius >
HIP_GLOBAL void batched_gram_m3_kernel(
    const LowT* HIP_RESTRICT r0,
    const LowT* HIP_RESTRICT r1,
    const LowT* HIP_RESTRICT r2,
    std::complex< double >* HIP_RESTRICT block_partials,
    int Nx, int Ny, int Nz,
    int Gx, int Gy, int Gz) {
    int b   = static_cast< int >(blockIdx.y);
    int tid = static_cast< int >(threadIdx.x);

    size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
    size_t state_stride = static_cast< size_t >(Gx) * plane_stride;
    size_t state_offset = static_cast< size_t >(b) * state_stride;

    const LowT* p0 = r0 + state_offset;
    const LowT* p1 = r1 + state_offset;
    const LowT* p2 = r2 + state_offset;

    double g00 = 0.0, g11 = 0.0, g22 = 0.0;
    std::complex< double > g01{0.0, 0.0}, g02{0.0, 0.0}, g12{0.0, 0.0};

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

        std::complex< double > v0 = to_complex64(p0[g_idx]);
        std::complex< double > v1 = to_complex64(p1[g_idx]);
        std::complex< double > v2 = to_complex64(p2[g_idx]);

        g00 += v0.real() * v0.real() + v0.imag() * v0.imag();
        g11 += v1.real() * v1.real() + v1.imag() * v1.imag();
        g22 += v2.real() * v2.real() + v2.imag() * v2.imag();

        g01 += std::conj(v0) * v1;
        g02 += std::conj(v0) * v2;
        g12 += std::conj(v1) * v2;
    }

    __shared__ double s_re[6 * DOT_BLOCK_SIZE];
    __shared__ double s_im[6 * DOT_BLOCK_SIZE];

    s_re[0 * DOT_BLOCK_SIZE + tid] = g00;        s_im[0 * DOT_BLOCK_SIZE + tid] = 0.0;
    s_re[1 * DOT_BLOCK_SIZE + tid] = g11;        s_im[1 * DOT_BLOCK_SIZE + tid] = 0.0;
    s_re[2 * DOT_BLOCK_SIZE + tid] = g22;        s_im[2 * DOT_BLOCK_SIZE + tid] = 0.0;
    s_re[3 * DOT_BLOCK_SIZE + tid] = g01.real(); s_im[3 * DOT_BLOCK_SIZE + tid] = g01.imag();
    s_re[4 * DOT_BLOCK_SIZE + tid] = g02.real(); s_im[4 * DOT_BLOCK_SIZE + tid] = g02.imag();
    s_re[5 * DOT_BLOCK_SIZE + tid] = g12.real(); s_im[5 * DOT_BLOCK_SIZE + tid] = g12.imag();
    __syncthreads();

    for (int s = DOT_BLOCK_SIZE / 2; s > 0; s >>= 1) {
        if (tid < s) {
            for (int k = 0; k < 6; ++k) {
                s_re[k * DOT_BLOCK_SIZE + tid] += s_re[k * DOT_BLOCK_SIZE + tid + s];
                s_im[k * DOT_BLOCK_SIZE + tid] += s_im[k * DOT_BLOCK_SIZE + tid + s];
            }
        }
        __syncthreads();
    }

    if (tid == 0) {
        size_t out_base = (static_cast< size_t >(b) * static_cast< size_t >(gridDim.x) + static_cast< size_t >(blockIdx.x)) * 6;
        for (int k = 0; k < 6; ++k) {
            block_partials[out_base + static_cast< size_t >(k)] =
std::complex< double >(s_re[k * DOT_BLOCK_SIZE], s_im[k * DOT_BLOCK_SIZE]);
}
}
}
template < typename LowT, int HaloRadius >
inline std::vector< std::array< std::complex< double >, 9 > > compute_batched_gram_m3(
const LowT* d_r0,
const LowT* d_r1,
const LowT* d_r2,
std::complex< double >* d_dot_partials,
int my_band_states,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz,
int num_blocks = 256,
hipStream_t stream = nullptr) {
std::vector< std::array< std::complex< double >, 9 > > result(static_cast< size_t >(std::max(0, my_band_states)));
if (my_band_states <= 0) return result;
dim3 threads(DOT_BLOCK_SIZE, 1, 1);
dim3 blocks(static_cast< unsigned int >(num_blocks), static_cast< unsigned int >(my_band_states), 1);
batched_gram_m3_kernel< LowT, HaloRadius ><<<blocks, threads, 0, stream>>>(
d_r0, d_r1, d_r2, d_dot_partials, Nx, Ny, Nz, Gx, Gy, Gz);
size_t total_partials = static_cast< size_t >(my_band_states) * static_cast< size_t >(num_blocks) * 6;
std::vector< std::complex< double > > h_partials(total_partials);
(void)hipMemcpyAsync(h_partials.data(), d_dot_partials,
total_partials * sizeof(std::complex< double >), hipMemcpyDeviceToHost, stream);
(void)hipStreamSynchronize(stream);
for (int b = 0; b < my_band_states; ++b) {
size_t ub = static_cast< size_t >(b);
std::array< std::complex< double >, 6 > sums{};
for (int blk = 0; blk < num_blocks; ++blk) {
size_t base = (ub * static_cast< size_t >(num_blocks) + static_cast< size_t >(blk)) * 6;
for (int k = 0; k < 6; ++k) {
sums.at(static_cast< size_t >(k)) += h_partials.at(base + static_cast< size_t >(k));
}
}
auto& G = result.at(ub);
G.at(0) = std::complex< double >(sums.at(0).real(), 0.0);
G.at(4) = std::complex< double >(sums.at(1).real(), 0.0);
G.at(8) = std::complex< double >(sums.at(2).real(), 0.0);
G.at(1) = sums.at(3);
G.at(3) = std::conj(sums.at(3));
G.at(2) = sums.at(4);
G.at(6) = std::conj(sums.at(4));
G.at(5) = sums.at(5);
G.at(7) = std::conj(sums.at(5));
}
return result;
}
inline void solve_small_diis_system(
int M,
const std::vector< std::complex< double > >& G_mat,
std::vector< std::complex< double > >& c_out) {
c_out.assign(static_cast< size_t >(M), std::complex< double >(0.0, 0.0));
if (M <= 0) return;
if (M == 1) {
c_out.at(0) = std::complex< double >(1.0, 0.0);
return;
}
int N = M + 1;
std::vector< std::complex< double > > A(static_cast< size_t >(N * N), std::complex< double >(0.0, 0.0));
std::vector< std::complex< double > > rhs(static_cast< size_t >(N), std::complex< double >(0.0, 0.0));
double max_diag = 0.0;
for (int i = 0; i < M; ++i) {
max_diag = std::max(max_diag, std::abs(G_mat.at(static_cast< size_t >(i * M + i))));
}
double reg = std::max(1e-18, 1e-3 * max_diag);
for (int i = 0; i < M; ++i) {
size_t ui = static_cast< size_t >(i);
for (int j = 0; j < M; ++j) {
size_t uj = static_cast< size_t >(j);
std::complex< double > val = G_mat.at(ui * static_cast< size_t >(M) + uj);
if (i == j) val += std::complex< double >(reg, 0.0);
A.at(ui * static_cast< size_t >(N) + uj) = val;
}
A.at(ui * static_cast< size_t >(N) + static_cast< size_t >(M)) = std::complex< double >(1.0, 0.0);
A.at(static_cast< size_t >(M) * static_cast< size_t >(N) + ui) = std::complex< double >(1.0, 0.0);
}
rhs.at(static_cast< size_t >(M)) = std::complex< double >(1.0, 0.0);
for (int k = 0; k < N; ++k) {
size_t uk = static_cast< size_t >(k);
int piv = k;
double max_v = std::abs(A.at(uk * static_cast< size_t >(N) + uk));
for (int i = k + 1; i < N; ++i) {
double v = std::abs(A.at(static_cast< size_t >(i) * static_cast< size_t >(N) + uk));
if (v > max_v) {
max_v = v;
piv = i;
}
}
if (!(max_v > 1e-30)) {
c_out.assign(static_cast< size_t >(M), std::complex< double >(0.0, 0.0));
c_out.at(static_cast< size_t >(M - 1)) = std::complex< double >(1.0, 0.0);
return;
}
if (piv != k) {
for (int j = k; j < N; ++j) {
std::swap(A.at(uk * static_cast< size_t >(N) + static_cast< size_t >(j)),
A.at(static_cast< size_t >(piv) * static_cast< size_t >(N) + static_cast< size_t >(j)));
}
std::swap(rhs.at(uk), rhs.at(static_cast< size_t >(piv)));
}
std::complex< double > diag = A.at(uk * static_cast< size_t >(N) + uk);
if (!(std::abs(diag) > 1e-30)) {
c_out.assign(static_cast< size_t >(M), std::complex< double >(0.0, 0.0));
c_out.at(static_cast< size_t >(M - 1)) = std::complex< double >(1.0, 0.0);
return;
}
for (int i = k + 1; i < N; ++i) {
size_t ui = static_cast< size_t >(i);
std::complex< double > factor = A.at(ui * static_cast< size_t >(N) + uk) / diag;
for (int j = k + 1; j < N; ++j) {
size_t uj = static_cast< size_t >(j);
A.at(ui * static_cast< size_t >(N) + uj) -= factor * A.at(uk * static_cast< size_t >(N) + uj);
}
rhs.at(ui) -= factor * rhs.at(uk);
}
}
std::vector< std::complex< double > > sol(static_cast< size_t >(N), std::complex< double >(0.0, 0.0));
for (int i = N - 1; i >= 0; --i) {
size_t ui = static_cast< size_t >(i);
std::complex< double > sum = rhs.at(ui);
for (int j = i + 1; j < N; ++j) {
size_t uj = static_cast< size_t >(j);
sum -= A.at(ui * static_cast< size_t >(N) + uj) * sol.at(uj);
}
std::complex< double > diag_i = A.at(ui * static_cast< size_t >(N) + ui);
if (!(std::abs(diag_i) > 1e-30)) {
c_out.assign(static_cast< size_t >(M), std::complex< double >(0.0, 0.0));
c_out.at(static_cast< size_t >(M - 1)) = std::complex< double >(1.0, 0.0);
return;
}
sol.at(ui) = sum / diag_i;
}
double max_c = 0.0;
for (int m = 0; m < M; ++m) {
c_out.at(static_cast< size_t >(m)) = sol.at(static_cast< size_t >(m));
max_c = std::max(max_c, std::abs(sol.at(static_cast< size_t >(m))));
}
if (!(max_c <= 4.0)) {
c_out.assign(static_cast< size_t >(M), std::complex< double >(0.0, 0.0));
c_out.at(static_cast< size_t >(M - 1)) = std::complex< double >(1.0, 0.0);
}
}
// ============================================================================
// 3. WAVEFUNCTION CORRECTION & RENORMALIZATION KERNELS
// ============================================================================
template < typename WaveT, int HaloRadius >
HIP_GLOBAL void apply_fp32_correction_kernel(
WaveT* HIP_RESTRICT psi,
const typename PrecisionTraits< WaveT >::LowT* HIP_RESTRICT corr_32,
float timestep,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz) {
using RealT = typename PrecisionTraits< WaveT >::RealT;
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
psi[idx] += static_cast< RealT >(timestep) * cast_from_low< WaveT >(corr_32[idx]);
}
template < typename WaveT, int HaloRadius >
inline void launch_apply_fp32_correction(
WaveT* d_psi,
const typename PrecisionTraits< WaveT >::LowT* d_corr_32,
int my_band_states,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz,
hipStream_t stream = nullptr,
float timestep = 0.6f) {
if (my_band_states <= 0) return;
dim3 threads(64, 4, 1);
dim3 blocks(static_cast< unsigned int >((Nz + 63) / 64),
static_cast< unsigned int >((Ny + 3) / 4),
static_cast< unsigned int >(my_band_states * Nx));
apply_fp32_correction_kernel< WaveT, HaloRadius ><<<blocks, threads, 0, stream>>>(
d_psi, d_corr_32, timestep, Nx, Ny, Nz, Gx, Gy, Gz);
}
template < typename WaveT, int HaloRadius >
HIP_GLOBAL void batched_norm_sq_kernel(
const WaveT* HIP_RESTRICT psi,
double* HIP_RESTRICT block_partials,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz) {
int b = static_cast< int >(blockIdx.y);
int tid = static_cast< int >(threadIdx.x);
size_t plane_stride = static_cast< size_t >(Gy) * static_cast< size_t >(Gz);
size_t state_stride = static_cast< size_t >(Gx) * plane_stride;
const WaveT* pb = psi + static_cast< size_t >(b) * state_stride;
double local_sum = 0.0;
int total_inner = Nx * Ny * Nz;
for (int idx = static_cast< int >(blockIdx.x) * static_cast< int >(blockDim.x) + tid;
idx < total_inner;
idx += static_cast< int >(gridDim.x) * static_cast< int >(blockDim.x)) {
int iz = idx % Nz;
int tmp = idx / Nz;
int iy = tmp % Ny;
int ix = tmp / Ny;
size_t g_idx = (static_cast< size_t >(ix + HaloRadius) * static_cast< size_t >(Gy)
+ static_cast< size_t >(iy + HaloRadius)) * static_cast< size_t >(Gz)
+ static_cast< size_t >(iz + HaloRadius);
std::complex< double > pv = to_complex64(pb[g_idx]);
local_sum += pv.real() * pv.real() + pv.imag() * pv.imag();
}
__shared__ double sdata[DOT_BLOCK_SIZE];
sdata[tid] = local_sum;
__syncthreads();
for (int s = DOT_BLOCK_SIZE / 2; s > 0; s >>= 1) {
if (tid < s) sdata[tid] += sdata[tid + s];
__syncthreads();
}
if (tid == 0) {
size_t out_idx = static_cast< size_t >(b) * static_cast< size_t >(gridDim.x) + static_cast< size_t >(blockIdx.x);
block_partials[out_idx] = sdata[0];
}
}
template < typename WaveT, int HaloRadius >
HIP_GLOBAL void batched_scale_wavefunction_kernel(
WaveT* HIP_RESTRICT psi,
const typename PrecisionTraits< WaveT >::RealT* HIP_RESTRICT inv_norms,
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
psi[idx] *= inv_norms[b];
}
template < typename WaveT, int HaloRadius >
inline void launch_renormalize_wavefunctions(
WaveT* d_psi,
double* d_norm_partials,
typename PrecisionTraits< WaveT >::RealT* d_inv_norms,
double vel,
int my_band_states,
int Nx, int Ny, int Nz,
int Gx, int Gy, int Gz,
int num_blocks = 256,
hipStream_t stream = nullptr) {
if (my_band_states <= 0) return;
using RealT = typename PrecisionTraits< WaveT >::RealT;
dim3 red_threads(DOT_BLOCK_SIZE, 1, 1);
dim3 red_blocks(static_cast< unsigned int >(num_blocks), static_cast< unsigned int >(my_band_states), 1);
batched_norm_sq_kernel< WaveT, HaloRadius ><<<red_blocks, red_threads, 0, stream>>>(
d_psi, d_norm_partials, Nx, Ny, Nz, Gx, Gy, Gz);
size_t total_partials = static_cast< size_t >(my_band_states) * static_cast< size_t >(num_blocks);
std::vector< double > h_partials(total_partials);
(void)hipMemcpyAsync(h_partials.data(), d_norm_partials,
total_partials * sizeof(double), hipMemcpyDeviceToHost, stream);
(void)hipStreamSynchronize(stream);
double scale_vol = (vel > 0.0) ? vel : 1.0;
std::vector< RealT > h_inv_norms(static_cast< size_t >(my_band_states), RealT{1});
for (int b = 0; b < my_band_states; ++b) {
double sum = 0.0;
size_t base = static_cast< size_t >(b) * static_cast< size_t >(num_blocks);
for (int blk = 0; blk < num_blocks; ++blk) {
sum += h_partials.at(base + static_cast< size_t >(blk));
}
if (sum > 0.0) {
h_inv_norms.at(static_cast< size_t >(b)) = static_cast< RealT >(1.0 / std::sqrt(sum * scale_vol));
}
}
(void)hipMemcpyAsync(d_inv_norms, h_inv_norms.data(),
static_cast< size_t >(my_band_states) * sizeof(RealT), hipMemcpyHostToDevice, stream);
dim3 scale_threads(64, 4, 1);
dim3 scale_blocks(static_cast< unsigned int >((Nz + 63) / 64),
static_cast< unsigned int >((Ny + 3) / 4),
static_cast< unsigned int >(my_band_states * Nx));
batched_scale_wavefunction_kernel< WaveT, HaloRadius ><<<scale_blocks, scale_threads, 0, stream>>>(
d_psi, d_inv_norms, Nx, Ny, Nz, Gx, Gy, Gz);
}

