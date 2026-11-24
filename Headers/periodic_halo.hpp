#pragma once

#include <cassert>
#include <cstdint>
#include <hip/hip_runtime.h>

#ifndef HIP_GLOBAL
#define HIP_GLOBAL __attribute__((global))
#endif

#ifndef HIP_RESTRICT
#define HIP_RESTRICT __restrict
#endif

// Sweep 1: Periodic X wrap
template <typename FieldT, int HaloRadius>
HIP_GLOBAL void fill_periodic_x_kernel(FieldT *HIP_RESTRICT field, int Nx,
                                       int Ny, int Nz, int Gx, int Gy, int Gz) {
    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int h = static_cast<int>(blockIdx.z) % HaloRadius;
    int b = static_cast<int>(blockIdx.z) / HaloRadius;

    if (iz >= Nz || iy >= Ny)
        return;

    size_t stride_y = static_cast<size_t>(Gz);
    size_t stride_x = static_cast<size_t>(Gy) * stride_y;
    size_t state_stride = static_cast<size_t>(Gx) * stride_x;

    FieldT *fb = field + static_cast<size_t>(b) * state_stride;
    size_t yz_off = static_cast<size_t>(iy + HaloRadius) * stride_y +
                    static_cast<size_t>(iz + HaloRadius);

    fb[static_cast<size_t>(h) * stride_x + yz_off] =
        fb[static_cast<size_t>(Nx + h) * stride_x + yz_off];
    fb[static_cast<size_t>(Nx + HaloRadius + h) * stride_x + yz_off] =
        fb[static_cast<size_t>(HaloRadius + h) * stride_x + yz_off];
}

// Sweep 2: Periodic Y wrap (across full Gx -> populates XY edges)
template <typename FieldT, int HaloRadius>
HIP_GLOBAL void fill_periodic_y_kernel(FieldT *HIP_RESTRICT field, int Ny,
                                       int Nz, int Gx, int Gy, int Gz) {
    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int gx = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int h = static_cast<int>(blockIdx.z) % HaloRadius;
    int b = static_cast<int>(blockIdx.z) / HaloRadius;

    if (iz >= Nz || gx >= Gx)
        return;

    size_t stride_y = static_cast<size_t>(Gz);
    size_t stride_x = static_cast<size_t>(Gy) * stride_y;
    size_t state_stride = static_cast<size_t>(Gx) * stride_x;

    FieldT *fb = field + static_cast<size_t>(b) * state_stride;
    size_t xz_off = static_cast<size_t>(gx) * stride_x +
                    static_cast<size_t>(iz + HaloRadius);

    fb[xz_off + static_cast<size_t>(h) * stride_y] =
        fb[xz_off + static_cast<size_t>(Ny + h) * stride_y];
    fb[xz_off + static_cast<size_t>(Ny + HaloRadius + h) * stride_y] =
        fb[xz_off + static_cast<size_t>(HaloRadius + h) * stride_y];
}

// Sweep 3: Periodic Z wrap (across full Gx, Gy -> populates XZ, YZ edges & XYZ
// corners)
template <typename FieldT, int HaloRadius>
HIP_GLOBAL void fill_periodic_z_kernel(FieldT *HIP_RESTRICT field, int Nz,
                                       int Gx, int Gy, int Gz) {
    int gy = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int gx = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int h = static_cast<int>(blockIdx.z) % HaloRadius;
    int b = static_cast<int>(blockIdx.z) / HaloRadius;

    if (gy >= Gy || gx >= Gx)
        return;

    size_t stride_y = static_cast<size_t>(Gz);
    size_t stride_x = static_cast<size_t>(Gy) * stride_y;
    size_t state_stride = static_cast<size_t>(Gx) * stride_x;

    FieldT *fb = field + static_cast<size_t>(b) * state_stride;
    size_t xy_off =
        static_cast<size_t>(gx) * stride_x + static_cast<size_t>(gy) * stride_y;

    fb[xy_off + static_cast<size_t>(h)] =
        fb[xy_off + static_cast<size_t>(Nz + h)];
    fb[xy_off + static_cast<size_t>(Nz + HaloRadius + h)] =
        fb[xy_off + static_cast<size_t>(HaloRadius + h)];
}

template <typename FieldT, int HaloRadius>
inline void launch_fill_periodic_halo(FieldT *d_field, int my_band_states,
                                      int Nx, int Ny, int Nz, int Gx, int Gy,
                                      int Gz, bool periodic = true,
                                      hipStream_t stream = nullptr) {
    if (my_band_states <= 0 || !periodic)
        return;
    assert(Gx == Nx + 2 * HaloRadius);
    assert(Gy == Ny + 2 * HaloRadius);
    assert(Gz == Nz + 2 * HaloRadius);

    {
        dim3 threads(64, 4, 1);
        dim3 blocks(static_cast<unsigned int>((Nz + 63) / 64),
                    static_cast<unsigned int>((Ny + 3) / 4),
                    static_cast<unsigned int>(my_band_states * HaloRadius));
        fill_periodic_x_kernel<FieldT, HaloRadius>
            <<<blocks, threads, 0, stream>>>(d_field, Nx, Ny, Nz, Gx, Gy, Gz);
    }
    {
        dim3 threads(64, 4, 1);
        dim3 blocks(static_cast<unsigned int>((Nz + 63) / 64),
                    static_cast<unsigned int>((Gx + 3) / 4),
                    static_cast<unsigned int>(my_band_states * HaloRadius));
        fill_periodic_y_kernel<FieldT, HaloRadius>
            <<<blocks, threads, 0, stream>>>(d_field, Ny, Nz, Gx, Gy, Gz);
    }
    {
        dim3 threads(64, 4, 1);
        dim3 blocks(static_cast<unsigned int>((Gy + 63) / 64),
                    static_cast<unsigned int>((Gx + 3) / 4),
                    static_cast<unsigned int>(my_band_states * HaloRadius));
        fill_periodic_z_kernel<FieldT, HaloRadius>
            <<<blocks, threads, 0, stream>>>(d_field, Nz, Gx, Gy, Gz);
    }
}
