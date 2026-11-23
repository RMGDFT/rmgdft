// Self-contained loader for RMGDFT's FiniteDiff::fd_combined_coeffs:
// Automatically computes the central coefficient c0 from the off-center taps + k^2
// for both fine (order = 8) and coarse (order = 2) grid levels.
inline std::vector< GpuStencilTap > build_rmg_taps_from_fd(
    FiniteDiff& fd,
    int order,
    double hxgrid,
    double hygrid,
    double hzgrid,
    const double* kvec,
    double key) {
    int radius = order / 2;
    // Allocate extra padding (+ 4) in case fd_combined_coeffs writes up to order elements
    std::vector< std::complex< double > > cm(static_cast< size_t >(order + 4), std::complex< double >(0.0, 0.0));
    std::vector< std::complex< double > > cp(static_cast< size_t >(order + 4), std::complex< double >(0.0, 0.0));
    std::vector< GpuStencilTap > taps;
    taps.reserve(static_cast< size_t >(1 + 6 * radius));

    // Safe k-vector copy (handles kvec == nullptr for Gamma-only runs)
    double k_safe[3] = {0.0, 0.0, 0.0};
    if (kvec != nullptr) {
        k_safe[0] = *(kvec + 0);
        k_safe[1] = *(kvec + 1);
        k_safe[2] = *(kvec + 2);
    }

    std::complex< double > neighbor_sum{0.0, 0.0};
    const double h_axes[3] = {hxgrid, hygrid, hzgrid};

    for (int axis = 0; axis < 3; ++axis) {
        std::fill(cm.begin(), cm.end(), std::complex< double >(0.0, 0.0));
        std::fill(cp.begin(), cp.end(), std::complex< double >(0.0, 0.0));

        fd.fd_combined_coeffs< std::complex< double > >(
            order, key, axis, cm.data(), cp.data(), k_safe);

        for (int idx = 0; idx < radius; ++idx) {
            int m = idx + 1;
            std::complex< double > wm = cm.at(static_cast< size_t >(idx));
            std::complex< double > wp = cp.at(static_cast< size_t >(idx));
            neighbor_sum += wm + wp;

            taps.push_back(GpuStencilTap{(axis == 0 ? -m : 0), (axis == 1 ? -m : 0), (axis == 2 ? -m : 0), wm});
            taps.push_back(GpuStencilTap{(axis == 0 ?  m : 0), (axis == 1 ?  m : 0), (axis == 2 ?  m : 0), wp});
        }
    }

    double kin_scale = (neighbor_sum.real() < 0.0) ? 0.5 : -0.5;
    for (auto& tap : taps) {
        tap.weight *= kin_scale;
    }

    double k2 = k_safe[0] * k_safe[0] + k_safe[1] * k_safe[1] + k_safe[2] * k_safe[2];
    std::complex< double > w0 = (-kin_scale * neighbor_sum) + std::complex< double >(0.5 * k2, 0.0);
    taps.insert(taps.begin(), GpuStencilTap{0, 0, 0, w0});

    return taps;
}

