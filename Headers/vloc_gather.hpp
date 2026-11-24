#pragma once

#include "gpu_solver_bridge.hpp"
#include <cstdint>
#include <mpi.h>
#include <type_traits>
#include <vector>

template <typename T> inline MPI_Datatype vloc_mpi_type() {
    if constexpr (std::is_same_v<T, double>)
        return MPI_DOUBLE;
    else
        return MPI_FLOAT;
}

// Gathers distributed local vtot sub-boxes [my_box.nx][my_box.ny][my_box.nz]
// across all MPI ranks into global padded h_vloc_padded [(Nx + 8) * (Ny + 8) *
// (Nz + 8)] with 4-cell periodic wrap-around. Set periodic = false for CLUSTER
// / Dirichlet boundary conditions (leaves 4-cell halo at 0).
template <typename VlocT>
inline void
gather_vloc_padded(MPI_Comm comm,
                   const VlocT *local_vtot, // Local domain sub-box [my_box.nx *
                                            // my_box.ny * my_box.nz], Z fastest
                   GpuSubBox3D my_box, int Nx, int Ny, int Nz,
                   std::vector<VlocT> &h_vloc_padded, bool periodic = true) {
    constexpr int HALO = 4;
    int Gx = Nx + 2 * HALO;
    int Gy = Ny + 2 * HALO;
    int Gz = Nz + 2 * HALO;

    size_t total_padded = static_cast<size_t>(Gx) * static_cast<size_t>(Gy) *
                          static_cast<size_t>(Gz);
    std::vector<VlocT> local_contrib(total_padded, VlocT{0});
    h_vloc_padded.assign(total_padded, VlocT{0});

    // 1. Place this rank's local sub-box into its interior coordinates in the
    // global padded grid
    for (int ix = 0; ix < my_box.nx; ++ix) {
        size_t gx = static_cast<size_t>(my_box.x0 + ix + HALO);
        for (int iy = 0; iy < my_box.ny; ++iy) {
            size_t gy = static_cast<size_t>(my_box.y0 + iy + HALO);
            for (int iz = 0; iz < my_box.nz; ++iz) {
                size_t gz = static_cast<size_t>(my_box.z0 + iz + HALO);

                size_t lidx =
                    (static_cast<size_t>(ix) * static_cast<size_t>(my_box.ny) +
                     static_cast<size_t>(iy)) *
                        static_cast<size_t>(my_box.nz) +
                    static_cast<size_t>(iz);
                size_t gidx = (gx * static_cast<size_t>(Gy) + gy) *
                                  static_cast<size_t>(Gz) +
                              gz;

                local_contrib.at(gidx) = *(local_vtot + lidx);
            }
        }
    }

    // 2. Global assembly across all domain ranks (< 1.5 ms for 256^3)
    MPI_Allreduce(local_contrib.data(), h_vloc_padded.data(),
                  static_cast<int>(total_padded), vloc_mpi_type<VlocT>(),
                  MPI_SUM, comm);

    if (!periodic)
        return;

    // 3. Fill 4-cell periodic wrap-around halos (X -> Y -> Z sweeps
    // automatically fill edges & corners)
    auto wrap_idx = [](int coord, int N) -> int {
        int r = coord % N;
        return (r < 0) ? (r + N) : r;
    };

    // X-halos (interior Y, Z)
    for (int gx = 0; gx < Gx; ++gx) {
        if (gx >= HALO && gx < Nx + HALO)
            continue;
        int src_gx = wrap_idx(gx - HALO, Nx) + HALO;
        for (int gy = HALO; gy < Ny + HALO; ++gy) {
            for (int gz = HALO; gz < Nz + HALO; ++gz) {
                size_t dst =
                    (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                     static_cast<size_t>(gy)) *
                        static_cast<size_t>(Gz) +
                    static_cast<size_t>(gz);
                size_t src =
                    (static_cast<size_t>(src_gx) * static_cast<size_t>(Gy) +
                     static_cast<size_t>(gy)) *
                        static_cast<size_t>(Gz) +
                    static_cast<size_t>(gz);
                h_vloc_padded.at(dst) = h_vloc_padded.at(src);
            }
        }
    }

    // Y-halos (all X including X-halos, interior Z)
    for (int gx = 0; gx < Gx; ++gx) {
        for (int gy = 0; gy < Gy; ++gy) {
            if (gy >= HALO && gy < Ny + HALO)
                continue;
            int src_gy = wrap_idx(gy - HALO, Ny) + HALO;
            for (int gz = HALO; gz < Nz + HALO; ++gz) {
                size_t dst =
                    (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                     static_cast<size_t>(gy)) *
                        static_cast<size_t>(Gz) +
                    static_cast<size_t>(gz);
                size_t src =
                    (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                     static_cast<size_t>(src_gy)) *
                        static_cast<size_t>(Gz) +
                    static_cast<size_t>(gz);
                h_vloc_padded.at(dst) = h_vloc_padded.at(src);
            }
        }
    }

    // Z-halos (all X and Y -> fills all edges and corners)
    for (int gx = 0; gx < Gx; ++gx) {
        for (int gy = 0; gy < Gy; ++gy) {
            for (int gz = 0; gz < Gz; ++gz) {
                if (gz >= HALO && gz < Nz + HALO)
                    continue;
                int src_gz = wrap_idx(gz - HALO, Nz) + HALO;
                size_t dst =
                    (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                     static_cast<size_t>(gy)) *
                        static_cast<size_t>(Gz) +
                    static_cast<size_t>(gz);
                size_t src =
                    (static_cast<size_t>(gx) * static_cast<size_t>(Gy) +
                     static_cast<size_t>(gy)) *
                        static_cast<size_t>(Gz) +
                    static_cast<size_t>(src_gz);
                h_vloc_padded.at(dst) = h_vloc_padded.at(src);
            }
        }
    }
}
