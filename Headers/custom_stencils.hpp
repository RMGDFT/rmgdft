#pragma once

#include "periodic_halo.hpp"
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <hip/hip_runtime.h>
#include <type_traits>
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

// ============================================================================
// 1. PRECISION TRAITS & COMPLEX HELPERS
// ============================================================================

template <typename T> struct PrecisionTraits;

template <> struct PrecisionTraits<float> {
    using RealT = float;
    using LowT = float;
    using WeightT = float;
    HIP_DEVICE static inline LowT to_low(float x) { return x; }
};

template <> struct PrecisionTraits<double> {
    using RealT = double;
    using LowT = float;
    using WeightT = double;
    HIP_DEVICE static inline LowT to_low(double x) {
        return static_cast<float>(x);
    }
};

template <> struct PrecisionTraits<std::complex<float>> {
    using RealT = float;
    using LowT = std::complex<float>;
    using WeightT = std::complex<float>;
    HIP_DEVICE static inline LowT to_low(std::complex<float> x) { return x; }
};

template <> struct PrecisionTraits<std::complex<double>> {
    using RealT = double;
    using LowT = std::complex<float>;
    using WeightT = std::complex<double>;
    HIP_DEVICE static inline LowT to_low(std::complex<double> x) {
        return std::complex<float>(static_cast<float>(x.real()),
                                   static_cast<float>(x.imag()));
    }
};

HIP_DEVICE inline float stencil_mul_add_real(float w, float v) { return w * v; }
HIP_DEVICE inline double stencil_mul_add_real(double w, double v) {
    return w * v;
}
HIP_DEVICE inline float stencil_mul_add_real(std::complex<float> w, float v) {
    return w.real() * v;
}
HIP_DEVICE inline double stencil_mul_add_real(std::complex<double> w,
                                              double v) {
    return w.real() * v;
}

HIP_DEVICE inline std::complex<float>
stencil_mul_add_complex(std::complex<float> w, std::complex<float> v) {
    return std::complex<float>(w.real() * v.real() - w.imag() * v.imag(),
                               w.real() * v.imag() + w.imag() * v.real());
}

HIP_DEVICE inline std::complex<double>
stencil_mul_add_complex(std::complex<double> w, std::complex<double> v) {
    return std::complex<double>(w.real() * v.real() - w.imag() * v.imag(),
                                w.real() * v.imag() + w.imag() * v.real());
}

HIP_DEVICE inline std::complex<float>
stencil_mul_add_complex(float w, std::complex<float> v) {
    return std::complex<float>(w * v.real(), w * v.imag());
}

HIP_DEVICE inline std::complex<double>
stencil_mul_add_complex(double w, std::complex<double> v) {
    return std::complex<double>(w * v.real(), w * v.imag());
}

// ============================================================================
// 2. STATELESS 3D STENCIL STRUCTURES
// ============================================================================

template <typename WeightT> struct StencilTap {
    int dx = 0;
    int dy = 0;
    int dz = 0;
    WeightT weight{};
};

template <typename WeightT, int MaxOffTaps> struct CustomStencil {
    WeightT center_weight{};
    int num_off_taps = 0;
    int dx[MaxOffTaps]{};
    int dy[MaxOffTaps]{};
    int dz[MaxOffTaps]{};
    WeightT weights[MaxOffTaps]{};
};

template <typename WeightT> using FineStencil8th = CustomStencil<WeightT, 104>;

template <typename WeightT> using CoarseStencil2nd = CustomStencil<WeightT, 20>;

template <typename OutWeightT>
inline OutWeightT convert_tap_weight(std::complex<double> w) {
    if constexpr (std::is_same_v<OutWeightT, float>) {
        return static_cast<float>(w.real());
    } else if constexpr (std::is_same_v<OutWeightT, double>) {
        return w.real();
    } else if constexpr (std::is_same_v<OutWeightT, std::complex<float>>) {
        return std::complex<float>(static_cast<float>(w.real()),
                                   static_cast<float>(w.imag()));
    } else {
        return w;
    }
}

template <typename WeightT, int HaloRadius, int MaxOffTaps>
inline CustomStencil<WeightT, MaxOffTaps> make_custom_stencil(
    const std::vector<StencilTap<std::complex<double>>> &user_taps) {
    CustomStencil<WeightT, MaxOffTaps> st{};
    int off_count = 0;

    for (const auto &tap : user_taps) {
        assert(tap.dx >= -HaloRadius && tap.dx <= HaloRadius);
        assert(tap.dy >= -HaloRadius && tap.dy <= HaloRadius);
        assert(tap.dz >= -HaloRadius && tap.dz <= HaloRadius);

        WeightT w_conv = convert_tap_weight<WeightT>(tap.weight);
        if (tap.dx == 0 && tap.dy == 0 && tap.dz == 0) {
            st.center_weight += w_conv;
        } else {
            assert(off_count < MaxOffTaps);
            st.dx[off_count] = tap.dx;
            st.dy[off_count] = tap.dy;
            st.dz[off_count] = tap.dz;
            st.weights[off_count] = w_conv;
            ++off_count;
        }
    }
    st.num_off_taps = off_count;
    return st;
}

template <typename FieldT, typename WeightT, int MaxOffTaps>
HIP_DEVICE inline FieldT
eval_stateless_stencil_3d(const FieldT *HIP_RESTRICT field_state, int gx,
                          int gy, int gz, int Gy, int Gz,
                          const CustomStencil<WeightT, MaxOffTaps> &stencil) {
    size_t c_idx = (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                    static_cast<size_t>(gy)) *
                       static_cast<size_t>(Gz) +
                   static_cast<size_t>(gz);

    FieldT acc{};
    if constexpr (std::is_same_v<FieldT, float> ||
                  std::is_same_v<FieldT, double>) {
        acc = stencil_mul_add_real(stencil.center_weight, field_state[c_idx]);
        for (int t = 0; t < stencil.num_off_taps; ++t) {
            int nx = gx + stencil.dx[t];
            int ny = gy + stencil.dy[t];
            int nz = gz + stencil.dz[t];
            size_t n_idx = (static_cast<size_t>(nx) * static_cast<size_t>(Gy) +
                            static_cast<size_t>(ny)) *
                               static_cast<size_t>(Gz) +
                           static_cast<size_t>(nz);
            acc += stencil_mul_add_real(stencil.weights[t], field_state[n_idx]);
        }
    } else {
        acc =
            stencil_mul_add_complex(stencil.center_weight, field_state[c_idx]);
        for (int t = 0; t < stencil.num_off_taps; ++t) {
            int nx = gx + stencil.dx[t];
            int ny = gy + stencil.dy[t];
            int nz = gz + stencil.dz[t];
            size_t n_idx = (static_cast<size_t>(nx) * static_cast<size_t>(Gy) +
                            static_cast<size_t>(ny)) *
                               static_cast<size_t>(Gz) +
                           static_cast<size_t>(nz);
            acc +=
                stencil_mul_add_complex(stencil.weights[t], field_state[n_idx]);
        }
    }
    return acc;
}

// ============================================================================
// 3. SHIFTED HAMILTONIAN & FINE-GRID RESIDUAL KERNELS
// ============================================================================

template <typename FieldT, int HaloRadius, int MaxOffTaps>
HIP_GLOBAL void apply_shifted_hamiltonian_kernel(
    const FieldT *HIP_RESTRICT u, FieldT *HIP_RESTRICT Au,
    const float *HIP_RESTRICT vloc, const float *HIP_RESTRICT eigvals, int Nx,
    int Ny, int Nz, int Gx, int Gy, int Gz,
    CustomStencil<typename PrecisionTraits<FieldT>::WeightT, MaxOffTaps>
        stencil) {
    using WeightT = typename PrecisionTraits<FieldT>::WeightT;

    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int ix = static_cast<int>(blockIdx.z) % Nx;
    int b = static_cast<int>(blockIdx.z) / Nx;

    if (iz >= Nz || iy >= Ny)
        return;

    int gx = ix + HaloRadius;
    int gy = iy + HaloRadius;
    int gz = iz + HaloRadius;

    size_t plane_stride = static_cast<size_t>(Gy) * static_cast<size_t>(Gz);
    size_t state_stride = static_cast<size_t>(Gx) * plane_stride;

    size_t grid_idx = (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                       static_cast<size_t>(gy)) *
                          static_cast<size_t>(Gz) +
                      static_cast<size_t>(gz);
    size_t idx = static_cast<size_t>(b) * state_stride + grid_idx;

    const FieldT *u_state = u + static_cast<size_t>(b) * state_stride;
    FieldT tu = eval_stateless_stencil_3d<FieldT, WeightT, MaxOffTaps>(
        u_state, gx, gy, gz, Gy, Gz, stencil);

    float shift = vloc[grid_idx] - eigvals[b];
    Au[idx] = tu + u[idx] * shift;
}

template <typename FieldT, int HaloRadius, int MaxOffTaps>
inline void launch_shifted_hamiltonian(
    FieldT *d_u, FieldT *d_Au, const float *d_vloc, const float *d_eigvals,
    int my_band_states, int Nx, int Ny, int Nz, int Gx, int Gy, int Gz,
    const CustomStencil<typename PrecisionTraits<FieldT>::WeightT, MaxOffTaps>
        &stencil,
    bool periodic = true, hipStream_t stream = nullptr) {
    if (my_band_states <= 0)
        return;

    launch_fill_periodic_halo<FieldT, HaloRadius>(
        d_u, my_band_states, Nx, Ny, Nz, Gx, Gy, Gz, periodic, stream);

    dim3 threads(64, 4, 1);
    dim3 blocks(static_cast<unsigned int>((Nz + 63) / 64),
                static_cast<unsigned int>((Ny + 3) / 4),
                static_cast<unsigned int>(my_band_states * Nx));

    apply_shifted_hamiltonian_kernel<FieldT, HaloRadius, MaxOffTaps>
        <<<blocks, threads, 0, stream>>>(d_u, d_Au, d_vloc, d_eigvals, Nx, Ny,
                                         Nz, Gx, Gy, Gz, stencil);
}

template <typename WaveT, int HaloRadius>
HIP_GLOBAL void fine_residual_and_downcast_kernel(
    const WaveT *HIP_RESTRICT psi, const WaveT *HIP_RESTRICT vnl,
    const typename PrecisionTraits<WaveT>::RealT *HIP_RESTRICT vloc,
    const typename PrecisionTraits<WaveT>::RealT *HIP_RESTRICT eigvals,
    typename PrecisionTraits<WaveT>::LowT *HIP_RESTRICT r_32,
    WaveT *HIP_RESTRICT r_wave, int Nx, int Ny, int Nz, int Gx, int Gy, int Gz,
    CustomStencil<typename PrecisionTraits<WaveT>::WeightT, 104> stencil) {
    using RealT = typename PrecisionTraits<WaveT>::RealT;
    using WeightT = typename PrecisionTraits<WaveT>::WeightT;
    int iz = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
             static_cast<int>(threadIdx.x);
    int iy = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
             static_cast<int>(threadIdx.y);
    int ix = static_cast<int>(blockIdx.z) % Nx;
    int b = static_cast<int>(blockIdx.z) / Nx;
    if (iz >= Nz || iy >= Ny)
        return;
    int gx = ix + HaloRadius;
    int gy = iy + HaloRadius;
    int gz = iz + HaloRadius;
    size_t plane_stride = static_cast<size_t>(Gy) * static_cast<size_t>(Gz);
    size_t state_stride = static_cast<size_t>(Gx) * plane_stride;
    size_t grid_idx = (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                       static_cast<size_t>(gy)) *
                          static_cast<size_t>(Gz) +
                      static_cast<size_t>(gz);
    size_t idx = static_cast<size_t>(b) * state_stride + grid_idx;
    const WaveT *psi_state = psi + static_cast<size_t>(b) * state_stride;
    WaveT tpsi = eval_stateless_stencil_3d<WaveT, WeightT, 104>(
        psi_state, gx, gy, gz, Gy, Gz, stencil);
    RealT v_val = vloc[grid_idx];
    RealT eig_val = eigvals[b];
    WaveT hpsi = tpsi + psi[idx] * v_val + vnl[idx];
    WaveT res = psi[idx] * eig_val - hpsi;
    if (r_wave != nullptr) {
        r_wave[idx] = res;
    }
    if (r_32 != nullptr) {
        r_32[idx] = PrecisionTraits<WaveT>::to_low(res);
    }
}
template <typename WaveT>
inline void launch_fine_residual_and_downcast(
    const WaveT *d_psi, const WaveT *d_vnl,
    const typename PrecisionTraits<WaveT>::RealT *d_vloc_wave,
    const typename PrecisionTraits<WaveT>::RealT *d_eig_wave,
    typename PrecisionTraits<WaveT>::LowT *d_r_32, int my_band_states, int Nx,
    int Ny, int Nz, int Gx, int Gy, int Gz,
    const std::vector<StencilTap<std::complex<double>>> &fine_taps,
    hipStream_t stream = nullptr, WaveT *d_r_wave = nullptr) {
    if (my_band_states <= 0)
        return;
    using WeightT = typename PrecisionTraits<WaveT>::WeightT;
    auto stencil = make_custom_stencil<WeightT, 4, 104>(fine_taps);
    dim3 threads(64, 4, 1);
    dim3 blocks(static_cast<unsigned int>((Nz + 63) / 64),
                static_cast<unsigned int>((Ny + 3) / 4),
                static_cast<unsigned int>(my_band_states * Nx));
    fine_residual_and_downcast_kernel<WaveT, 4><<<blocks, threads, 0, stream>>>(
        d_psi, d_vnl, d_vloc_wave, d_eig_wave, d_r_32, d_r_wave, Nx, Ny, Nz, Gx,
        Gy, Gz, stencil);
}
