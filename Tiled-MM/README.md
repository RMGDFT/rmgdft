


## Table of Contents
- [Overview](#overview)
- [Performance](#performance)
- [Features](#features)
- [Building and Installing](#building-and-installing)
- [Minimal Working Example](#minimal-working-example)
- [Running the Benchmarks](#running-the-benchmarks)
- [Testing](#testing)
- [Author](#author)


## Overview

Tiled-MM is a very fast and easy-to-use library for multiplying matrices on GPU. As opposed to NVIDIA's `cublas`, this library takes pointer from the host side (CPU), splits the matrices into tiles, pipelines them efficiently to the GPU and copies the result back to the CPU. It can serve as almost a drop-in replacement for `cublasXt`, and is ported to both NVIDIA and AMD gpus.

It offers more features than the standard cublas API. For example, the user can specify the number of gpu streams to be used, as well as the tile size for each dimension separately, which is not possible with the standard cublas API.

Tiled-MM is used in production as a backend of the [COSMA](https://github.com/eth-cscs/COSMA) algorithm and is thus well-tested.

## Performance

The benchmarks were performed on a single node of Piz Daint Supercomputer (Cray XC50), equipped with a `P100` NVIDIA GPU. We compared the performance of our library `Tiled-MM` with the vanilla version of `cublasXt` and also with the manually tuned version of `cublasXt`, where we manually set the tile size to `4000` and enabled the pinned memory mode. `Tiled-MM` was substantially faster than the vanilla version of `cublasXt`, and achieved similar performance as the manually tuned version of `cublasXt`, as can be seen from the results below.
<p align="center"><img src="https://github.com/eth-cscs/Tiled-MM/blob/master/docs/performance.svg" width="90%"></p>

In the benchmark, we used `double precision`, `square matrices` given in `column-major` ordering, and `alpha = beta = 1.0`.

## Features:

- The user can specify the tile size of each dimension separately.
- The user can specify the number of streams to be used.
- The user can reuse the same context (and thus the same device memory) for many multiplications which can lead to significant performance improvements.
- Fully templatized, supporting arbitrary data types.
- Ported to both `NVIDIA` and `AMD` GPUs.

## Building and Installing

Assuming that you want to use the `gcc 8` compiler, you can build the project as follows:
```bash
# clone the repo
git clone --recursive https://github.com/eth-cscs/Tiled-MM
cd Tiled-MM
mkdir build
cd build

# build
CC=gcc-8 CXX=g++-8 cmake -DTILEDMM_GPU_BACKEND=CUDA ..

# compile
make -j 4
```

The option `-DTILEDMM_GPU_BACKEND` can have the following values:
- `CUDA`: for NVIDIA GPUs
- `ROCM`: for AMD GPUs

## Minimal Working Example

Using the library is very simple, just include `#include <tiled_mm.hpp>` and use it as follows:
```cpp
// A dimensions: m x k
auto a_host = gpu::malloc_pinned<double>(m * k, 1);
// B dimensions: k x n
auto b_host = gpu::malloc_pinned<double>(k * n, 1);
// C dimensions: m x n
auto c_host = gpu::malloc_pinned<double>(m * n, 0);

double alpha = 1.0;
double beta = 0.0;

// preallocates device buffers and other CUDA stuff
// the context does not have to be created explicitly
// so the user can omit this part
auto ctx = gpu::make_context();

// compute c = alpha * a * b + beta * c
// There is also a version without ctx, in case the user
// does not want to create the context explicitly
gpu::gemm(*ctx,
          trans_a, trans_b,
          m, n, k,
          alpha,
          a_host, ld_a,
          b_host, ld_b,
          beta,
          c_host, ld_c);

// optionally, we can set the following two boolean flags
bool pin_buffers = false; // since a_host, b_host and c_host are already pinned, gpu::dgemm should not pin them
bool copy_c_back = true;  // if we want to copy the result back to the host or leave it on the gpu
gpu::gemm(*ctx,
          trans_a, trans_b,
          m, n, k,
          alpha,
          a_host, ld_a,
          b_host, ld_b,
          beta,
          c_host, ld_c,
          pin_buffers, copy_c_back);

// if copy_c_back == false, the result is stored on the device with the following pointer:
double* c_device = ctx->get_full_device_buffer_c().data()
```
When creating the context, the user can specify tile dimensions and the number of streams to be used as:
```cpp
int tile_size_m = 5000;
int tile_size_n = 5000;
int tile_size_k = 5000;
int n_streams = 2;

auto ctx = gpu::make_context(n_streams, tile_size_m, tile_size_n, tile_size_k);
```
## Running the Benchmarks

For detailed benchmarking, there is a miniapp that takes the host pointers for A, B and C and computes `C = beta * C + alpha * A * B` outputing the time-to-solution, as well as the throughput.

The miniapp consists of the executable `./build/examples/multiply` which can be run with the following command line (assuming we are in the root folder of the project):
```bash
./build/examples/multiply -m 10000 -n 10000 -k 10000 -r 1
```
The overview of all supported options is given below:
Option Flags | POSSIBLE VALUES | DESCRIPTION
| :------------------- | :------------------- |:------------------- |
`m (--m_dim)` | positive integer | Number of rows of `C`
`n (--n_dim)` | positive integer | Number of columns of `C`
`k (--k_dim)` | positive integer | size of the shared dimension between matrices `A` and `B`
`--tile_m` | positive integer | tile size for dimension `m`
`--tile_n` | positive integer | tile size for dimension `n`
`--tile_k` | positive integer | tile size for dimension `k`
`--ld_a` | positive integer | leading dimension of matrix `A`
`--ld_b` | positive integer | leading dimension of matrix `B`
`--ld_c` | positive integer | leading dimension of matrix `C`
`-t (--transpose)` | a string XY, where X, Y can be one of {N, T, C} | transpose flags for matrices A and B
`--alpha` | real value (double) | the `alpha` in `C = beta * C + alpha * A * B`
`--beta` | real value (double) | the `beta` in `C = beta * C + alpha * A * B`

For example, running with the following flags:
```bash
./build/examples/multiply -m 1000 -n 1000 -k 1000 --transpose=TN -r 1
```
should produce the following output:
```bash
==================================================
                Benchmarking Tiled-MM
==================================================
         MATRIX SIZES
=============================
 A = (1000, 1000)
 B = (1000, 1000)
 C = (1000, 1000)
=============================
         LEADING DIMS
=============================
 LD_A = 1000
 LD_B = 1000
 LD_C = 1000
=============================
      SCALING CONSTANTS
=============================
 alpha = 1
 beta  = 1
=============================
      TRANSPOSE FLAGS
=============================
 trans_a = T
 trans_b = N
=============================
         TILE SIZES
=============================
 tile_m = 5000
 tile_n = 5000
 tile_k = 5000
=============================
      ADDITIONAL OPTIONS
=============================
 num. of gpu streams = 2
 num. of repetitions = 1
=============================

==================================================
         Results of benchmarking Tiled-MM
==================================================
 1) The version with copying C to back to host:
    -> Avg Time [ms] = 11
    -> Throughput [Gflops] = 181.818
==================================================
 2) The version without copying C to back to host:
    -> Avg Time [ms] = 10
    -> Throughput [Gflops] = 200
==================================================
```

## Testing

For testing purposes, there is a testing miniapp that generates random matrices A, B and C, computes `C = beta * C + alpha * A * B` with Tiled-MM as well as with blas and outputs whether the results are correct.

The miniapp consists of the executable `./build/tests/test-multiply` **supports the same parameters** as the benchmarking miniapp (see above). It can be run e.g. with the following command line (assuming we are in the root folder of the project):
```bash
./build/tests/test-multiply -m 1000 -n 1000 -k 1000 --transpose=TN
```
which should produce the following output:
```bash 
==================================================
                Benchmarking Tiled-MM
==================================================
         MATRIX SIZES
=============================
 A = (1000, 1000)
 B = (1000, 1000)
 C = (1000, 1000)
=============================
         LEADING DIMS
=============================
 LD_A = 1000
 LD_B = 1000
 LD_C = 1000
=============================
      SCALING CONSTANTS
=============================
 alpha = 1
 beta  = 1
=============================
      TRANSPOSE FLAGS
=============================
 trans_a = T
 trans_b = N
=============================
         TILE SIZES
=============================
 tile_m = 5000
 tile_n = 5000
 tile_k = 5000
=============================
      ADDITIONAL OPTIONS
=============================
 num. of gpu streams = 2
 num. of repetitions = 1
=============================
Time [ms] with copying C back: 11
Time [ms] without copying C back: 10
The result is CORRECT
```
Running `make test` will few default tests.

## Author
Marko Kabic (marko.kabic@inf.ethz.ch)
