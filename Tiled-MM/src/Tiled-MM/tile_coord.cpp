#include "tile_coord.hpp"

namespace gpu {

tile_coord::tile_coord(int tile_idx_row, int tile_idx_col):
        tile_index_row(tile_idx_row), tile_index_col(tile_idx_col) {}

int tile_coord::row_index() {
    return tile_index_row;
}

int tile_coord::col_index() {
    return tile_index_col;
}

}
