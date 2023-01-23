#pragma once
namespace gpu {

struct tile_coord {
public:
    tile_coord() = default;
    tile_coord(int tile_idx_row, int tile_idx_col);

    int row_index();
    int col_index();

private:
    int tile_index_row;
    int tile_index_col;
};

}
