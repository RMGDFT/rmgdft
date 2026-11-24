#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <hip/hip_runtime.h>
#include <mpi.h>
#include <vector>

#ifndef HIP_GLOBAL
#define HIP_GLOBAL __attribute__((global))
#endif

#ifndef HIP_RESTRICT
#define HIP_RESTRICT __restrict
#endif

template <typename T> inline MPI_Datatype get_mpi_type();

template <> inline MPI_Datatype get_mpi_type<double>() { return MPI_DOUBLE; }

template <> inline MPI_Datatype get_mpi_type<float>() { return MPI_FLOAT; }

template <> inline MPI_Datatype get_mpi_type<std::complex<double>>() {
#if defined(MPI_CXX_DOUBLE_COMPLEX)
    return MPI_CXX_DOUBLE_COMPLEX;
#elif defined(MPI_C_DOUBLE_COMPLEX)
    return MPI_C_DOUBLE_COMPLEX;
#else
    static MPI_Datatype complex_type = MPI_DATATYPE_NULL;
    if (complex_type == MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(2, MPI_DOUBLE, &complex_type);
        MPI_Type_commit(&complex_type);
    }
    return complex_type;
#endif
}

template <> inline MPI_Datatype get_mpi_type<std::complex<float>>() {
#if defined(MPI_CXX_FLOAT_COMPLEX)
    return MPI_CXX_FLOAT_COMPLEX;
#elif defined(MPI_C_FLOAT_COMPLEX)
    return MPI_C_FLOAT_COMPLEX;
#else
    static MPI_Datatype complex_float_type = MPI_DATATYPE_NULL;
    if (complex_float_type == MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(2, MPI_FLOAT, &complex_float_type);
        MPI_Type_commit(&complex_float_type);
    }
    return complex_float_type;
#endif
}

struct SubBox3D {
    int x0 = 0, y0 = 0, z0 = 0; // Interior global start coordinates (0-based,
                                // excluding global ghost layers)
    int nx = 0, ny = 0, nz = 0; // Interior local sub-box dimensions
};

struct Halo3D {
    int hx = 0; // Ghost/boundary layers on each side of X
    int hy = 0; // Ghost/boundary layers on each side of Y
    int hz = 0; // Ghost/boundary layers on each side of Z
};

// Single-array 3D Cartesian scatter/gather kernel (z-axis unit stride, global
// ghost layers)
template <typename T>
HIP_GLOBAL void transpose_3d_cartesian_halo_kernel(
    const T *HIP_RESTRICT src, T *HIP_RESTRICT dst,
    const SubBox3D *HIP_RESTRICT boxes, const size_t *HIP_RESTRICT elem_offsets,
    int n_peers, int b_local, int Gx, int Gy, int Gz, int hx, int hy, int hz,
    bool to_cartesian) {
    int p = static_cast<int>(blockIdx.y);
    int b = static_cast<int>(blockIdx.z);
    if (p >= n_peers || b >= b_local)
        return;

    SubBox3D box = *(boxes + p);
    size_t p_count = static_cast<size_t>(box.nx) * static_cast<size_t>(box.ny) *
                     static_cast<size_t>(box.nz);
    size_t idx =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    if (idx >= p_count)
        return;

    int z_loc = static_cast<int>(idx % static_cast<size_t>(box.nz));
    size_t tmp = idx / static_cast<size_t>(box.nz);
    int y_loc = static_cast<int>(tmp % static_cast<size_t>(box.ny));
    int x_loc = static_cast<int>(tmp / static_cast<size_t>(box.ny));

    size_t gx = static_cast<size_t>(box.x0 + x_loc + hx);
    size_t gy = static_cast<size_t>(box.y0 + y_loc + hy);
    size_t gz = static_cast<size_t>(box.z0 + z_loc + hz);

    size_t cart_idx = ((static_cast<size_t>(b) * static_cast<size_t>(Gx) + gx) *
                           static_cast<size_t>(Gy) +
                       gy) *
                          static_cast<size_t>(Gz) +
                      gz;
    size_t pack_idx = (*(elem_offsets + p)) * static_cast<size_t>(b_local) +
                      static_cast<size_t>(b) * p_count + idx;

    if (to_cartesian) {
        *(dst + cart_idx) = *(src + pack_idx);
    } else {
        *(dst + pack_idx) = *(src + cart_idx);
    }
}

// Dual-array 3D Cartesian scatter kernel: unpacks interleaved [psi_peer_block |
// vnl_peer_block] into two GPU buffers
template <typename T>
HIP_GLOBAL void scatter_3d_cartesian_pair_kernel(
    const T *HIP_RESTRICT src_pair, T *HIP_RESTRICT dst_a,
    T *HIP_RESTRICT dst_b, const SubBox3D *HIP_RESTRICT boxes,
    const size_t *HIP_RESTRICT elem_offsets, int n_peers, int b_local, int Gx,
    int Gy, int Gz, int hx, int hy, int hz) {
    int p = static_cast<int>(blockIdx.y);
    int b = static_cast<int>(blockIdx.z);
    if (p >= n_peers || b >= b_local)
        return;

    SubBox3D box = *(boxes + p);
    size_t p_count = static_cast<size_t>(box.nx) * static_cast<size_t>(box.ny) *
                     static_cast<size_t>(box.nz);
    size_t idx =
        static_cast<size_t>(blockIdx.x) * static_cast<size_t>(blockDim.x) +
        static_cast<size_t>(threadIdx.x);
    if (idx >= p_count)
        return;

    int z_loc = static_cast<int>(idx % static_cast<size_t>(box.nz));
    size_t tmp = idx / static_cast<size_t>(box.nz);
    int y_loc = static_cast<int>(tmp % static_cast<size_t>(box.ny));
    int x_loc = static_cast<int>(tmp / static_cast<size_t>(box.ny));

    size_t gx = static_cast<size_t>(box.x0 + x_loc + hx);
    size_t gy = static_cast<size_t>(box.y0 + y_loc + hy);
    size_t gz = static_cast<size_t>(box.z0 + z_loc + hz);

    size_t cart_idx = ((static_cast<size_t>(b) * static_cast<size_t>(Gx) + gx) *
                           static_cast<size_t>(Gy) +
                       gy) *
                          static_cast<size_t>(Gz) +
                      gz;

    size_t peer_base = (*(elem_offsets + p)) * static_cast<size_t>(2 * b_local);
    size_t slice_off = static_cast<size_t>(b) * p_count + idx;
    size_t half_peer = static_cast<size_t>(b_local) * p_count;

    *(dst_a + cart_idx) = *(src_pair + peer_base + slice_off);
    *(dst_b + cart_idx) = *(src_pair + peer_base + half_peer + slice_off);
}

template <typename T> class CartesianWavefunctionTransposer {
  public:
    CartesianWavefunctionTransposer(MPI_Comm parent_comm, int total_states,
                                    int Nx, int Ny, int Nz, SubBox3D my_box,
                                    Halo3D gpu_global_halo = Halo3D{0, 0, 0},
                                    int desired_batch_states = 1024,
                                    bool enable_pair_buffers = true)
        : total_states_(total_states), Nx_(Nx), Ny_(Ny), Nz_(Nz),
          my_box_(my_box), halo_(gpu_global_halo),
          pair_enabled_(enable_pair_buffers) {
        MPI_Comm_dup(parent_comm, &comm_);
        MPI_Comm_size(comm_, &n_ranks_);
        MPI_Comm_rank(comm_, &rank_);
        assert(total_states_ > 0);
        assert(desired_batch_states > 0);

        if (n_ranks_ <= desired_batch_states && n_ranks_ <= total_states_) {
            b_local_ = std::max(1, desired_batch_states / n_ranks_);
            batch_states_ = b_local_ * n_ranks_;
            n_uniform_batches_ = total_states_ / batch_states_;
            rem_states_ = total_states_ % batch_states_;
        } else {
            b_local_ = 0;
            batch_states_ = total_states_;
            n_uniform_batches_ = 0;
            rem_states_ = total_states_;
        }

        my_g_local_ = static_cast<size_t>(my_box_.nx) *
                      static_cast<size_t>(my_box_.ny) *
                      static_cast<size_t>(my_box_.nz);
        interior_grid_ = static_cast<size_t>(Nx_) * static_cast<size_t>(Ny_) *
                         static_cast<size_t>(Nz_);

        Gx_ = Nx_ + 2 * halo_.hx;
        Gy_ = Ny_ + 2 * halo_.hy;
        Gz_ = Nz_ + 2 * halo_.hz;
        gpu_allocated_grid_ = static_cast<size_t>(Gx_) *
                              static_cast<size_t>(Gy_) *
                              static_cast<size_t>(Gz_);

        boxes_.resize(static_cast<size_t>(n_ranks_));
        MPI_Allgather(&my_box_, static_cast<int>(sizeof(SubBox3D)), MPI_BYTE,
                      boxes_.data(), static_cast<int>(sizeof(SubBox3D)),
                      MPI_BYTE, comm_);

        elem_offsets_.resize(static_cast<size_t>(n_ranks_ + 1), 0);
        max_g_local_ = 0;
        for (int p = 0; p < n_ranks_; ++p) {
            size_t up = static_cast<size_t>(p);
            size_t cnt = static_cast<size_t>(boxes_.at(up).nx) *
                         static_cast<size_t>(boxes_.at(up).ny) *
                         static_cast<size_t>(boxes_.at(up).nz);
            elem_offsets_.at(up + 1) = elem_offsets_.at(up) + cnt;
            max_g_local_ = std::max(max_g_local_, cnt);
        }
        assert(elem_offsets_.at(static_cast<size_t>(n_ranks_)) ==
               interior_grid_);

        (void)hipMalloc(&d_boxes_,
                        static_cast<size_t>(n_ranks_) * sizeof(SubBox3D));
        (void)hipMemcpy(d_boxes_, boxes_.data(),
                        static_cast<size_t>(n_ranks_) * sizeof(SubBox3D),
                        hipMemcpyHostToDevice);

        (void)hipMalloc(&d_elem_offsets_,
                        static_cast<size_t>(n_ranks_ + 1) * sizeof(size_t));
        (void)hipMemcpy(d_elem_offsets_, elem_offsets_.data(),
                        static_cast<size_t>(n_ranks_ + 1) * sizeof(size_t),
                        hipMemcpyHostToDevice);

        build_alltoallv_tables();

        size_t max_stage_states_domain = static_cast<size_t>(
            std::max(n_uniform_batches_ > 0 ? batch_states_ : 0, rem_states_));
        int max_b_local_band = std::max(b_local_, my_rem_local_);
        size_t mult = pair_enabled_ ? 2 : 1;

        size_t domain_stage_bytes = std::max<size_t>(
            1, mult * max_stage_states_domain * my_g_local_ * sizeof(T));
        size_t band_stage_bytes =
            std::max<size_t>(1, mult * static_cast<size_t>(max_b_local_band) *
                                    interior_grid_ * sizeof(T));

        size_t num_slots = (n_uniform_batches_ > 1) ? 2 : 1;
        for (size_t i = 0; i < 2; ++i) {
            if (i < num_slots) {
                (void)hipMalloc(&d_domain_stage_.at(i), domain_stage_bytes);
                (void)hipMalloc(&d_band_stage_.at(i), band_stage_bytes);
            }
            (void)hipStreamCreate(&copy_stream_.at(i));
            (void)hipStreamCreate(&work_stream_.at(i));
        }
    }
    CartesianWavefunctionTransposer() {
        (void)hipDeviceSynchronize();
        (void)hipFree(d_boxes_);
        (void)hipFree(d_elem_offsets_);
        for (size_t i = 0; i < 2; ++i) {
            if (d_domain_stage_.at(i) != nullptr)
                (void)hipFree(d_domain_stage_.at(i));
            if (d_band_stage_.at(i) != nullptr)
                (void)hipFree(d_band_stage_.at(i));
            (void)hipStreamDestroy(copy_stream_.at(i));
            (void)hipStreamDestroy(work_stream_.at(i));
        }
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized && comm_ != MPI_COMM_NULL) {
            MPI_Comm_free(&comm_);
        }
    }
    CartesianWavefunctionTransposer(const CartesianWavefunctionTransposer &) =
        delete;
    CartesianWavefunctionTransposer &
    operator=(const CartesianWavefunctionTransposer &) = delete;
    void forward_pair(const T *h_pinned_a, const T *h_pinned_b, T *d_gpu_band_a,
                      T *d_gpu_band_b) {
        assert(pair_enabled_);
        size_t domain_batch_elems =
            static_cast<size_t>(batch_states_) * my_g_local_;
        if (n_uniform_batches_ > 0 && domain_batch_elems > 0) {
            stage_pair_h2d(h_pinned_a, h_pinned_b, d_domain_stage_.at(0),
                           b_local_, copy_stream_.at(0));
        }
        for (int step = 0; step < n_uniform_batches_; ++step) {
            size_t cur = static_cast<size_t>(step % 2);
            size_t nxt = static_cast<size_t>((step + 1) % 2);
            (void)hipStreamSynchronize(copy_stream_.at(cur));
            if (step + 1 < n_uniform_batches_ && domain_batch_elems > 0) {
                (void)hipStreamSynchronize(work_stream_.at(nxt));
                size_t next_off =
                    static_cast<size_t>(step + 1) * domain_batch_elems;
                stage_pair_h2d(h_pinned_a + next_off, h_pinned_b + next_off,
                               d_domain_stage_.at(nxt), b_local_,
                               copy_stream_.at(nxt));
            }
            MPI_Alltoallv(d_domain_stage_.at(cur), u2_sendcounts_.data(),
                          u2_sdispls_.data(), get_mpi_type<T>(),
                          d_band_stage_.at(cur), u2_recvcounts_.data(),
                          u2_rdispls_.data(), get_mpi_type<T>(), comm_);
            size_t band_off =
                static_cast<size_t>(step * b_local_) * gpu_allocated_grid_;
            launch_3d_scatter_pair(
                d_band_stage_.at(cur), d_gpu_band_a + band_off,
                d_gpu_band_b + band_off, b_local_, work_stream_.at(cur));
        }
        (void)hipDeviceSynchronize();
        if (rem_states_ > 0) {
            forward_pair_remainder(h_pinned_a, h_pinned_b, d_gpu_band_a,
                                   d_gpu_band_b);
            (void)hipDeviceSynchronize();
        }
    }
    void forward(const T *h_pinned_domain, T *d_gpu_cartesian_band) {
        size_t domain_batch_elems =
            static_cast<size_t>(batch_states_) * my_g_local_;
        if (n_uniform_batches_ > 0 && domain_batch_elems > 0) {
            (void)hipMemcpyAsync(d_domain_stage_.at(0), h_pinned_domain,
                                 domain_batch_elems * sizeof(T),
                                 hipMemcpyHostToDevice, copy_stream_.at(0));
        }
        for (int step = 0; step < n_uniform_batches_; ++step) {
            size_t cur = static_cast<size_t>(step % 2);
            size_t nxt = static_cast<size_t>((step + 1) % 2);
            (void)hipStreamSynchronize(copy_stream_.at(cur));
            if (step + 1 < n_uniform_batches_ && domain_batch_elems > 0) {
                (void)hipStreamSynchronize(work_stream_.at(nxt));
                const T *h_nxt =
                    h_pinned_domain +
                    static_cast<size_t>(step + 1) * domain_batch_elems;
                (void)hipMemcpyAsync(d_domain_stage_.at(nxt), h_nxt,
                                     domain_batch_elems * sizeof(T),
                                     hipMemcpyHostToDevice,
                                     copy_stream_.at(nxt));
            }
            MPI_Alltoallv(d_domain_stage_.at(cur), u_sendcounts_.data(),
                          u_sdispls_.data(), get_mpi_type<T>(),
                          d_band_stage_.at(cur), u_recvcounts_.data(),
                          u_rdispls_.data(), get_mpi_type<T>(), comm_);
            T *d_dst =
                d_gpu_cartesian_band +
                static_cast<size_t>(step * b_local_) * gpu_allocated_grid_;
            launch_3d_transpose(d_band_stage_.at(cur), d_dst, b_local_, true,
                                work_stream_.at(cur));
        }
        (void)hipDeviceSynchronize();
        if (rem_states_ > 0) {
            forward_remainder(h_pinned_domain, d_gpu_cartesian_band);
            (void)hipDeviceSynchronize();
        }
    }
    void inverse(const T *d_gpu_cartesian_band, T *h_pinned_domain) {
        size_t domain_batch_elems =
            static_cast<size_t>(batch_states_) * my_g_local_;
        if (n_uniform_batches_ > 0) {
            launch_3d_transpose(d_gpu_cartesian_band, d_band_stage_.at(0),
                                b_local_, false, work_stream_.at(0));
        }
        for (int step = 0; step < n_uniform_batches_; ++step) {
            size_t cur = static_cast<size_t>(step % 2);
            size_t nxt = static_cast<size_t>((step + 1) % 2);
            (void)hipStreamSynchronize(work_stream_.at(cur));
            if (step + 1 < n_uniform_batches_) {
                (void)hipStreamSynchronize(copy_stream_.at(nxt));
                const T *d_nxt = d_gpu_cartesian_band +
                                 static_cast<size_t>((step + 1) * b_local_) *
                                     gpu_allocated_grid_;
                launch_3d_transpose(d_nxt, d_band_stage_.at(nxt), b_local_,
                                    false, work_stream_.at(nxt));
            }
            MPI_Alltoallv(d_band_stage_.at(cur), u_recvcounts_.data(),
                          u_rdispls_.data(), get_mpi_type<T>(),
                          d_domain_stage_.at(cur), u_sendcounts_.data(),
                          u_sdispls_.data(), get_mpi_type<T>(), comm_);
            if (domain_batch_elems > 0) {
                T *h_dst = h_pinned_domain +
                           static_cast<size_t>(step) * domain_batch_elems;
                (void)hipMemcpyAsync(h_dst, d_domain_stage_.at(cur),
                                     domain_batch_elems * sizeof(T),
                                     hipMemcpyDeviceToHost,
                                     copy_stream_.at(cur));
            }
        }
        (void)hipDeviceSynchronize();
        if (rem_states_ > 0) {
            inverse_remainder(d_gpu_cartesian_band, h_pinned_domain);
            (void)hipDeviceSynchronize();
        }
    }
    int my_band_states() const {
        return (total_states_ / n_ranks_) +
               (rank_ < (total_states_ % n_ranks_) ? 1 : 0);
    }
    size_t gpu_band_buffer_elements() const {
        return static_cast<size_t>(my_band_states()) * gpu_allocated_grid_;
    }
    int Gx() const { return Gx_; }
    int Gy() const { return Gy_; }
    int Gz() const { return Gz_; }

  private:
    void stage_pair_h2d(const T *h_a, const T *h_b, T *d_pair_dst,
                        int states_per_peer, hipStream_t stream) {
        size_t peer_bytes =
            static_cast<size_t>(states_per_peer) * my_g_local_ * sizeof(T);
        if (peer_bytes == 0)
            return;
        size_t dst_pitch = 2 * peer_bytes;
        size_t height = static_cast<size_t>(n_ranks_);
        (void)hipMemcpy2DAsync(d_pair_dst, dst_pitch, h_a, peer_bytes,
                               peer_bytes, height, hipMemcpyHostToDevice,
                               stream);
        (void)hipMemcpy2DAsync(
            reinterpret_cast<char *>(d_pair_dst) + peer_bytes, dst_pitch, h_b,
            peer_bytes, peer_bytes, height, hipMemcpyHostToDevice, stream);
    }
    void launch_3d_transpose(const T *src, T *dst, int b_loc, bool to_cartesian,
                             hipStream_t stream) {
        if (b_loc <= 0 || max_g_local_ == 0)
            return;
        dim3 threads(256);
        unsigned int grid_x = static_cast<unsigned int>(
            (max_g_local_ + threads.x - 1) / threads.x);
        dim3 blocks(grid_x, static_cast<unsigned int>(n_ranks_),
                    static_cast<unsigned int>(b_loc));
        transpose_3d_cartesian_halo_kernel<T><<<blocks, threads, 0, stream>>>(
            src, dst, d_boxes_, d_elem_offsets_, n_ranks_, b_loc, Gx_, Gy_, Gz_,
            halo_.hx, halo_.hy, halo_.hz, to_cartesian);
    }
    void launch_3d_scatter_pair(const T *src_pair, T *dst_a, T *dst_b,
                                int b_loc, hipStream_t stream) {
        if (b_loc <= 0 || max_g_local_ == 0)
            return;
        dim3 threads(256);
        unsigned int grid_x = static_cast<unsigned int>(
            (max_g_local_ + threads.x - 1) / threads.x);
        dim3 blocks(grid_x, static_cast<unsigned int>(n_ranks_),
                    static_cast<unsigned int>(b_loc));
        scatter_3d_cartesian_pair_kernel<T><<<blocks, threads, 0, stream>>>(
            src_pair, dst_a, dst_b, d_boxes_, d_elem_offsets_, n_ranks_, b_loc,
            Gx_, Gy_, Gz_, halo_.hx, halo_.hy, halo_.hz);
    }
    void build_alltoallv_tables() {
        size_t nr = static_cast<size_t>(n_ranks_);
        if (n_uniform_batches_ > 0) {
            u_sendcounts_.resize(nr);
            u_sdispls_.resize(nr);
            u_recvcounts_.resize(nr);
            u_rdispls_.resize(nr);
            u2_sendcounts_.resize(nr);
            u2_sdispls_.resize(nr);
            u2_recvcounts_.resize(nr);
            u2_rdispls_.resize(nr);
            for (size_t p = 0; p < nr; ++p) {
                size_t p_cnt = elem_offsets_.at(p + 1) - elem_offsets_.at(p);
                u_sendcounts_.at(p) = b_local_ * static_cast<int>(my_g_local_);
                u_sdispls_.at(p) = static_cast<int>(p) * u_sendcounts_.at(p);
                u_recvcounts_.at(p) = b_local_ * static_cast<int>(p_cnt);
                u_rdispls_.at(p) =
                    b_local_ * static_cast<int>(elem_offsets_.at(p));
                u2_sendcounts_.at(p) = 2 * u_sendcounts_.at(p);
                u2_sdispls_.at(p) = 2 * u_sdispls_.at(p);
                u2_recvcounts_.at(p) = 2 * u_recvcounts_.at(p);
                u2_rdispls_.at(p) = 2 * u_rdispls_.at(p);
            }
        }
        if (rem_states_ == 0) {
            my_rem_local_ = 0;
            return;
        }
        int base = rem_states_ / n_ranks_;
        int extra = rem_states_ % n_ranks_;
        my_rem_local_ = base + (rank_ < extra ? 1 : 0);
        r_sendcounts_.resize(nr);
        r_sdispls_.resize(nr);
        r_recvcounts_.resize(nr);
        r_rdispls_.resize(nr);
        r2_sendcounts_.resize(nr);
        r2_sdispls_.resize(nr);
        r2_recvcounts_.resize(nr);
        r2_rdispls_.resize(nr);
        int s_off = 0;
        for (size_t p = 0; p < nr; ++p) {
            int p_states = base + (static_cast<int>(p) < extra ? 1 : 0);
            size_t p_cnt = elem_offsets_.at(p + 1) - elem_offsets_.at(p);
            r_sendcounts_.at(p) = p_states * static_cast<int>(my_g_local_);
            r_sdispls_.at(p) = s_off;
            s_off += r_sendcounts_.at(p);
            r_recvcounts_.at(p) = my_rem_local_ * static_cast<int>(p_cnt);
            r_rdispls_.at(p) =
                my_rem_local_ * static_cast<int>(elem_offsets_.at(p));
            r2_sendcounts_.at(p) = 2 * r_sendcounts_.at(p);
            r2_sdispls_.at(p) = 2 * r_sdispls_.at(p);
            r2_recvcounts_.at(p) = 2 * r_recvcounts_.at(p);
            r2_rdispls_.at(p) = 2 * r_rdispls_.at(p);
        }
    }
    void forward_pair_remainder(const T *h_pinned_a, const T *h_pinned_b,
                                T *d_gpu_band_a, T *d_gpu_band_b) {
        size_t offset_elems =
            static_cast<size_t>(n_uniform_batches_ * batch_states_) *
            my_g_local_;
        size_t band_offset =
            static_cast<size_t>(n_uniform_batches_ * b_local_) *
            gpu_allocated_grid_;
        for (size_t p = 0; p < static_cast<size_t>(n_ranks_); ++p) {
            size_t h_off = offset_elems + static_cast<size_t>(r_sdispls_.at(p));
            size_t d_off = static_cast<size_t>(r2_sdispls_.at(p));
            size_t cnt = static_cast<size_t>(r_sendcounts_.at(p));
            if (cnt > 0) {
                (void)hipMemcpyAsync(d_domain_stage_.at(0) + d_off,
                                     h_pinned_a + h_off, cnt * sizeof(T),
                                     hipMemcpyHostToDevice, copy_stream_.at(0));
                (void)hipMemcpyAsync(d_domain_stage_.at(0) + d_off + cnt,
                                     h_pinned_b + h_off, cnt * sizeof(T),
                                     hipMemcpyHostToDevice, copy_stream_.at(0));
            }
        }
        (void)hipStreamSynchronize(copy_stream_.at(0));
        MPI_Alltoallv(d_domain_stage_.at(0), r2_sendcounts_.data(),
                      r2_sdispls_.data(), get_mpi_type<T>(),
                      d_band_stage_.at(0), r2_recvcounts_.data(),
                      r2_rdispls_.data(), get_mpi_type<T>(), comm_);
        if (my_rem_local_ > 0) {
            launch_3d_scatter_pair(
                d_band_stage_.at(0), d_gpu_band_a + band_offset,
                d_gpu_band_b + band_offset, my_rem_local_, work_stream_.at(0));
        }
    }
    void forward_remainder(const T *h_pinned_domain, T *d_gpu_cartesian_band) {
        size_t offset_elems =
            static_cast<size_t>(n_uniform_batches_ * batch_states_) *
            my_g_local_;
        size_t rem_elems = static_cast<size_t>(rem_states_) * my_g_local_;
        size_t band_offset =
            static_cast<size_t>(n_uniform_batches_ * b_local_) *
            gpu_allocated_grid_;
        if (rem_elems > 0) {
            (void)hipMemcpyAsync(d_domain_stage_.at(0),
                                 h_pinned_domain + offset_elems,
                                 rem_elems * sizeof(T), hipMemcpyHostToDevice,
                                 copy_stream_.at(0));
            (void)hipStreamSynchronize(copy_stream_.at(0));
        }
        MPI_Alltoallv(d_domain_stage_.at(0), r_sendcounts_.data(),
                      r_sdispls_.data(), get_mpi_type<T>(), d_band_stage_.at(0),
                      r_recvcounts_.data(), r_rdispls_.data(),
                      get_mpi_type<T>(), comm_);
        if (my_rem_local_ > 0) {
            launch_3d_transpose(d_band_stage_.at(0),
                                d_gpu_cartesian_band + band_offset,
                                my_rem_local_, true, work_stream_.at(0));
        }
    }
    void inverse_remainder(const T *d_gpu_cartesian_band, T *h_pinned_domain) {
        size_t offset_elems =
            static_cast<size_t>(n_uniform_batches_ * batch_states_) *
            my_g_local_;
        size_t rem_elems = static_cast<size_t>(rem_states_) * my_g_local_;
        size_t band_offset =
            static_cast<size_t>(n_uniform_batches_ * b_local_) *
            gpu_allocated_grid_;
        if (my_rem_local_ > 0) {
            launch_3d_transpose(d_gpu_cartesian_band + band_offset,
                                d_band_stage_.at(0), my_rem_local_, false,
                                work_stream_.at(0));
            (void)hipStreamSynchronize(work_stream_.at(0));
        }
        MPI_Alltoallv(d_band_stage_.at(0), r_recvcounts_.data(),
                      r_rdispls_.data(), get_mpi_type<T>(),
                      d_domain_stage_.at(0), r_sendcounts_.data(),
                      r_sdispls_.data(), get_mpi_type<T>(), comm_);
        if (rem_elems > 0) {
            (void)hipMemcpyAsync(h_pinned_domain + offset_elems,
                                 d_domain_stage_.at(0), rem_elems * sizeof(T),
                                 hipMemcpyDeviceToHost, copy_stream_.at(0));
        }
    }
    MPI_Comm comm_ = MPI_COMM_NULL;
    int n_ranks_ = 0;
    int rank_ = 0;
    int total_states_ = 0;
    int batch_states_ = 0;
    int b_local_ = 0;
    int n_uniform_batches_ = 0;
    int rem_states_ = 0;
    int my_rem_local_ = 0;
    int Nx_ = 0, Ny_ = 0, Nz_ = 0;
    int Gx_ = 0, Gy_ = 0, Gz_ = 0;
    SubBox3D my_box_{};
    Halo3D halo_{};
    bool pair_enabled_ = true;
    size_t my_g_local_ = 0;
    size_t max_g_local_ = 0;
    size_t interior_grid_ = 0;
    size_t gpu_allocated_grid_ = 0;
    std::vector<SubBox3D> boxes_;
    std::vector<size_t> elem_offsets_;
    SubBox3D *d_boxes_ = nullptr;
    size_t *d_elem_offsets_ = nullptr;
    std::array<T *, 2> d_domain_stage_ = {{nullptr, nullptr}};
    std::array<T *, 2> d_band_stage_ = {{nullptr, nullptr}};
    std::array<hipStream_t, 2> copy_stream_ = {{nullptr, nullptr}};
    std::array<hipStream_t, 2> work_stream_ = {{nullptr, nullptr}};
    std::vector<int> u_sendcounts_, u_sdispls_, u_recvcounts_, u_rdispls_;
    std::vector<int> u2_sendcounts_, u2_sdispls_, u2_recvcounts_, u2_rdispls_;
    std::vector<int> r_sendcounts_, r_sdispls_, r_recvcounts_, r_rdispls_;
    std::vector<int> r2_sendcounts_, r2_sdispls_, r2_recvcounts_, r2_rdispls_;
};
