#pragma once

#include <cassert>
#include <complex>
#include <cstdint>
#include <hip/hip_runtime.h>

template <typename LowT>
__device__ __forceinline__ LowT scale_real_fp32(float s, LowT v);

template <>
__device__ __forceinline__ float scale_real_fp32<float>(float s, float v) {
    return s * v;
}

template <>
__device__ __forceinline__ std::complex<float>
scale_real_fp32<std::complex<float>>(float s, std::complex<float> v) {
    return std::complex<float>(s * v.real(), s * v.imag());
}

__device__ __forceinline__ int wrap_periodic_idx(int i, int n) {
    int r = i % n;
    return (r < 0) ? (r + n) : r;
}

// 1. 27-Point Full-Weighting Restriction Kernel (Fine -> Coarse)
// Wraps all 26 neighbors (6 faces + 12 edges + 8 corners) periodically when
// periodic == true, or zero-pads exterior boundaries when periodic == false.
template <typename LowT, int HaloFine, int HaloCoarse>
__global__ void restrict_full_weighting_3d_kernel(
    const LowT *__restrict__ r_fine, LowT *__restrict__ r_coarse, int Nx_f,
    int Ny_f, int Nz_f, int Gx_f, int Gy_f, int Gz_f, int Nx_c, int Ny_c,
    int Nz_c, int Gx_c, int Gy_c, int Gz_c, bool periodic) {
    int iz_c = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
               static_cast<int>(threadIdx.x);
    int iy_c = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
               static_cast<int>(threadIdx.y);
    int ix_c = static_cast<int>(blockIdx.z) % Nx_c;
    int b = static_cast<int>(blockIdx.z) / Nx_c;

    if (iz_c >= Nz_c || iy_c >= Ny_c)
        return;

    size_t plane_f = static_cast<size_t>(Gy_f) * static_cast<size_t>(Gz_f);
    size_t state_f = static_cast<size_t>(Gx_f) * plane_f;

    size_t plane_c = static_cast<size_t>(Gy_c) * static_cast<size_t>(Gz_c);
    size_t state_c = static_cast<size_t>(Gx_c) * plane_c;

    const LowT *rf_b = r_fine + static_cast<size_t>(b) * state_f;
    LowT *rc_b = r_coarse + static_cast<size_t>(b) * state_c;

    int ix0 = 2 * ix_c;
    int iy0 = 2 * iy_c;
    int iz0 = 2 * iz_c;

    LowT acc{};

// 3x3x3 full-weighting stencil: weights 1/4, 1/2, 1/4 along each axis
#pragma unroll
    for (int dx = -1; dx <= 1; ++dx) {
        int ix = ix0 + dx;
        if (periodic) {
            ix = wrap_periodic_idx(ix, Nx_f);
        } else if (ix < 0 || ix >= Nx_f) {
            continue;
        }
        int fx = ix + HaloFine;
        float wx = (dx == 0) ? 0.5f : 0.25f;

#pragma unroll
        for (int dy = -1; dy <= 1; ++dy) {
            int iy = iy0 + dy;
            if (periodic) {
                iy = wrap_periodic_idx(iy, Ny_f);
            } else if (iy < 0 || iy >= Ny_f) {
                continue;
            }
            int fy = iy + HaloFine;
            float wxy = wx * ((dy == 0) ? 0.5f : 0.25f);

            size_t row_base =
                static_cast<size_t>(fx) * plane_f +
                static_cast<size_t>(fy) * static_cast<size_t>(Gz_f);

#pragma unroll
            for (int dz = -1; dz <= 1; ++dz) {
                int iz = iz0 + dz;
                if (periodic) {
                    iz = wrap_periodic_idx(iz, Nz_f);
                } else if (iz < 0 || iz >= Nz_f) {
                    continue;
                }
                int fz = iz + HaloFine;
                float wxyz = wxy * ((dz == 0) ? 0.5f : 0.25f);

                acc += scale_real_fp32<LowT>(
                    wxyz, rf_b[row_base + static_cast<size_t>(fz)]);
            }
        }
    }

    size_t c_idx =
        (static_cast<size_t>(ix_c + HaloCoarse) * static_cast<size_t>(Gy_c) +
         static_cast<size_t>(iy_c + HaloCoarse)) *
            static_cast<size_t>(Gz_c) +
        static_cast<size_t>(iz_c + HaloCoarse);

    rc_b[c_idx] = acc;
}

// 2. Trilinear Prolongation + In-Place Addition Kernel (Coarse -> Fine)
// Wraps faces, edges, and corners periodically when periodic == true,
// or zero-pads exterior boundaries when periodic == false.
template <typename LowT, int HaloFine, int HaloCoarse>
__global__ void prolongate_add_3d_kernel(const LowT *__restrict__ e_coarse,
                                         LowT *__restrict__ e_fine, int Nx_f,
                                         int Ny_f, int Nz_f, int Gx_f, int Gy_f,
                                         int Gz_f, int Nx_c, int Ny_c, int Nz_c,
                                         int Gx_c, int Gy_c, int Gz_c,
                                         bool periodic) {
    int iz_f = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
               static_cast<int>(threadIdx.x);
    int iy_f = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
               static_cast<int>(threadIdx.y);
    int ix_f = static_cast<int>(blockIdx.z) % Nx_f;
    int b = static_cast<int>(blockIdx.z) / Nx_f;

    if (iz_f >= Nz_f || iy_f >= Ny_f)
        return;

    size_t plane_c = static_cast<size_t>(Gy_c) * static_cast<size_t>(Gz_c);
    size_t state_c = static_cast<size_t>(Gx_c) * plane_c;

    size_t plane_f = static_cast<size_t>(Gy_f) * static_cast<size_t>(Gz_f);
    size_t state_f = static_cast<size_t>(Gx_f) * plane_f;

    const LowT *ec_b = e_coarse + static_cast<size_t>(b) * state_c;
    LowT *ef_b = e_fine + static_cast<size_t>(b) * state_f;

    int ix_c0 = ix_f / 2;
    int iy_c0 = iy_f / 2;
    int iz_c0 = iz_f / 2;

    int px = ix_f & 1;
    int py = iy_f & 1;
    int pz = iz_f & 1;

    LowT interp{};

#pragma unroll
    for (int dx = 0; dx <= px; ++dx) {
        int ixc = ix_c0 + dx;
        if (periodic) {
            ixc = wrap_periodic_idx(ixc, Nx_c);
        } else if (ixc >= Nx_c) {
            continue;
        }
        int cx = ixc + HaloCoarse;
        float wx = (px == 0) ? 1.0f : 0.5f;

#pragma unroll
        for (int dy = 0; dy <= py; ++dy) {
            int iyc = iy_c0 + dy;
            if (periodic) {
                iyc = wrap_periodic_idx(iyc, Ny_c);
            } else if (iyc >= Ny_c) {
                continue;
            }
            int cy = iyc + HaloCoarse;
            float wxy = wx * ((py == 0) ? 1.0f : 0.5f);
            size_t row_c = static_cast<size_t>(cx) * plane_c +
                           static_cast<size_t>(cy) * static_cast<size_t>(Gz_c);

#pragma unroll
            for (int dz = 0; dz <= pz; ++dz) {
                int izc = iz_c0 + dz;
                if (periodic) {
                    izc = wrap_periodic_idx(izc, Nz_c);
                } else if (izc >= Nz_c) {
                    continue;
                }
                int cz = izc + HaloCoarse;
                float wxyz = wxy * ((pz == 0) ? 1.0f : 0.5f);
                interp += scale_real_fp32<LowT>(
                    wxyz, ec_b[row_c + static_cast<size_t>(cz)]);
            }
        }
    }

    size_t f_idx =
        (static_cast<size_t>(ix_f + HaloFine) * static_cast<size_t>(Gy_f) +
         static_cast<size_t>(iy_f + HaloFine)) *
            static_cast<size_t>(Gz_f) +
        static_cast<size_t>(iz_f + HaloFine);

    ef_b[f_idx] += interp;
}

// Host Launchers
template <typename LowT, int HaloFine, int HaloCoarse>
void launch_restrict_3d(const LowT *d_r_fine, LowT *d_r_coarse,
                        int my_band_states, int Nx_f, int Ny_f, int Nz_f,
                        int Gx_f, int Gy_f, int Gz_f, int Nx_c, int Ny_c,
                        int Nz_c, int Gx_c, int Gy_c, int Gz_c,
                        bool periodic = true, hipStream_t stream = nullptr) {
    if (my_band_states <= 0)
        return;
    assert(Gx_f == Nx_f + 2 * HaloFine);
    assert(Gy_f == Ny_f + 2 * HaloFine);
    assert(Gz_f == Nz_f + 2 * HaloFine);
    assert(Gx_c == Nx_c + 2 * HaloCoarse);
    assert(Gy_c == Ny_c + 2 * HaloCoarse);
    assert(Gz_c == Nz_c + 2 * HaloCoarse);

    dim3 threads(32, 8, 1);
    dim3 blocks(static_cast<unsigned int>((Nz_c + 31) / 32),
                static_cast<unsigned int>((Ny_c + 7) / 8),
                static_cast<unsigned int>(my_band_states * Nx_c));

    restrict_full_weighting_3d_kernel<LowT, HaloFine, HaloCoarse>
        <<<blocks, threads, 0, stream>>>(d_r_fine, d_r_coarse, Nx_f, Ny_f, Nz_f,
                                         Gx_f, Gy_f, Gz_f, Nx_c, Ny_c, Nz_c,
                                         Gx_c, Gy_c, Gz_c, periodic);
}

template <typename LowT, int HaloFine, int HaloCoarse>
void launch_prolongate_add_3d(const LowT *d_e_coarse, LowT *d_e_fine,
                              int my_band_states, int Nx_f, int Ny_f, int Nz_f,
                              int Gx_f, int Gy_f, int Gz_f, int Nx_c, int Ny_c,
                              int Nz_c, int Gx_c, int Gy_c, int Gz_c,
                              bool periodic = true,
                              hipStream_t stream = nullptr) {
    if (my_band_states <= 0)
        return;

    dim3 threads(32, 8, 1);
    dim3 blocks(static_cast<unsigned int>((Nz_f + 31) / 32),
                static_cast<unsigned int>((Ny_f + 7) / 8),
                static_cast<unsigned int>(my_band_states * Nx_f));

    prolongate_add_3d_kernel<LowT, HaloFine, HaloCoarse>
        <<<blocks, threads, 0, stream>>>(d_e_coarse, d_e_fine, Nx_f, Ny_f, Nz_f,
                                         Gx_f, Gy_f, Gz_f, Nx_c, Ny_c, Nz_c,
                                         Gx_c, Gy_c, Gz_c, periodic);
}
