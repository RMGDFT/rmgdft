#ifndef CUTRANSPOSE_H_
#define CUTRANSPOSE_H_

#define TILE_SIZE 16
#define BRICK_SIZE 8

#if HIP_ENABLED
#define  gpuStream_t  hipStream_t
#endif
#if CUDA_ENABLED
#define  gpuStream_t  cudaStream_t
#endif

/********************************************
 * Public function prototypes               *
 ********************************************/

int cut_transpose3d( float*       output,
                     const float* input,
                     const int*    size,
                     const int*    permutation,
                     int           elements_per_thread,
                     gpuStream_t stream );

int cut_transpose3d( double*       output,
                     const double* input,
                     const int*    size,
                     const int*    permutation,
                     int           elements_per_thread,
                     gpuStream_t stream );

#endif /* CUTRANSPOSE_H_ */
