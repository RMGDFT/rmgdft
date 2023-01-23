#pragma once

#include "mm_handle.hpp"

/*
  **************************************
         TILE-COPYING METHODS
  **************************************

  Matrix in host memory is given in column-major form.
  Matrix from host memory is split into tiles and each
  tile is copied to device memory.

  A tile is a small block of the host memory,
  which is also in column major order.

  **************************************
   Example:
  **************************************
  Matrix:
  1 2 5 4
  4 2 3 5
  1 6 7 4
  9 0 1 2

  Host:
  pointer to -> 1 4 1 9 2 2 6 0 5 3 7 1 4 5 4 2

  Device (assuming tile dimensions are 2x2):
  Tile 00: pointer to -> 1 4 2 2
  Tile 10: pointer to -> 1 9 6 0
  Tile 01: pointer to -> 5 3 4 5
  Tile 11: pointer to -> 7 1 4 2
       ^
  Tile id = <row_id><col_id>

  **************************************
      MEMORY WITH MULTIPLE-STREAMS
  **************************************
  On the device, N_STREAMS*tile_size memory is preallocated,
  such that each stream has a separate piece of memory.

  device memory (assuming 3 streams)
            -----------> stream offset
            ___________________________________
  array:   | TILE SIZE | TILE SIZE | TILE SIZE |
            ___________________________________
              stream 1    stream 2    stream 3
          ^
    device pointer

  Observe that the device pointer points to the beginning
  of the pre-allocated device buffer. Therefore, each stream
  has to add offset stream_id * TILE_SIZE to get the pointer
  to device tile it is in charge of.
  */


/*
    copy_tile methods copy a single tile between host<->device.
*/
namespace gpu {

blas_api::OperationType get_blas_operation(char trans);

// **************************************
//        TILED-GEMM ON GPU
// **************************************
template<typename Scalar>
void gemm(
	mm_handle<Scalar>& handle,
	char trans_a, char trans_b, 
        int m, int n, int k,
	Scalar alpha, 
	Scalar* a, int ld_a,
	Scalar* b, int ld_b,
	Scalar beta, 
        Scalar* c, int ld_c,
	bool pin_host_buffers = true, bool copy_c_back=true);
}
