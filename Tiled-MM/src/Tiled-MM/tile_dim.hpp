#pragma once
namespace gpu {

class tile_dim {
public:
    tile_dim() = default;
    tile_dim(int n_rows, int n_cols);

    int rows();

    int cols();

    void set_rows(int n_rows);

    void set_cols(int n_cols);

    int size();

private:
    int n_rows_;
    int n_cols_;
    int size_;
};

}
