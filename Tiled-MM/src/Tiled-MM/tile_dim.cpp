#include "tile_dim.hpp"

using namespace gpu;

tile_dim::tile_dim(int n_rows, int n_cols):
    n_rows_(n_rows), n_cols_(n_cols) {
    size_ = n_rows * n_cols;
}

int tile_dim::rows() {
    return n_rows_;
}

int tile_dim::cols() {
    return n_cols_;
}

void tile_dim::set_rows(int n_rows) {
    n_rows_ = n_rows;
}

void tile_dim::set_cols(int n_cols) {
    n_cols_ = n_cols;
}

int tile_dim::size() {
    return size_;
}

