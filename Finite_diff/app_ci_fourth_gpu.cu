#include "make_conf.h"
#if GPU_FD_ENABLED
#include "fixed_dims.h"
#ifdef FD_XSIZE


__global__ void app_cil_fourth_kernel1(const double * __restrict__ psi, 
                                                double *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const double cc,
                                                const double fcx,
                                                const double ecxz)
{

    __shared__ double slice[2*FIXED_YDIM/3 + 2][2*FIXED_ZDIM/3 + 2];
    double accm_b1, accm, accm_a1;
    double acc, acc_a1;
    double x_b1, x_0, x_a1;

    // iz and iy map to the x and y coordinates of the thread
    // within a block
    int iz = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 1;

    // thread z-index into shared memory tile
    int tz = threadIdx.x + 1;
    // thread y-index into shared memory tile
    int ty = threadIdx.y + 1;
    int ix=0;

    int incx = (dimz + 2) * (dimy + 2);
    int incy = dimz + 2;
    int incxr = dimz * dimy;
    int incyr = dimz;

    int ixs, ixps;
    int iys, iyms, iyps;
    int izs, izms, izps;

    ixs = ix * incx;
    ixps = (ix + 1) * incx;

    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;
    izs = iz;
    izms = (iz - 1);
    izps = (iz + 1);

    x_0 = 0.0;
    accm = 0.0;
    acc = 0.0;

    x_a1 = psi[izs + iys + ixs];

    acc_a1 = cc * psi[izs + iys + ixs] +
            ecxz * psi[izms + iyms + ixs] +
            ecxz * psi[izms + iyps + ixs] +
            ecxz * psi[izps + iyms + ixs] +
            ecxz * psi[izps + iyps + ixs];


    accm_a1 = psi[izms + iys + ixs] +
              psi[izps + iys + ixs] +
              psi[izs + iyms + ixs] +
              psi[izs + iyps + ixs];


    for (ix = 0; ix < dimx + 1; ix++)
    {

        x_b1 = x_0;
        x_0 = x_a1;
         
        accm_b1 = accm;
        accm = accm_a1;

        acc = acc_a1;

        __syncthreads();
        // Update the data slice in shared memory
        if(threadIdx.x == 0) {
            slice[ty][0] =
                            psi[(ix + 1)*incx + iy*incy + (blockIdx.x*blockDim.x)];
            slice[ty][threadIdx.x + blockDim.x + 1] =
                            psi[(ix + 1)*incx + iy*incy + (blockDim.x + 1 + blockIdx.x*blockDim.x)];
        }
        if(threadIdx.y == 0) {
            slice[0][tz] =
                            psi[(ix + 1)*incx + (blockIdx.y*blockDim.y)*incy + iz];
            slice[blockDim.y + 1][tz] =
                            psi[(ix + 1)*incx + (blockDim.y + 1 + blockIdx.y*blockDim.y)*incy + iz];
        }


        if((threadIdx.x == 2) && (threadIdx.y == 2)) {
            slice[0][0] = psi[(ix + 1)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x];
        }

        if((threadIdx.x == 2) && (threadIdx.y == 3)) {
            slice[0][blockDim.x + 1] = psi[(ix + 1)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 1];
        }

        if((threadIdx.x == 3) && (threadIdx.y == 2)) {
            slice[blockDim.y + 1][0] = psi[(ix + 1)*incx + (blockIdx.y*blockDim.y + blockDim.y + 1)*incy + blockIdx.x*blockDim.x];
        }

        if((threadIdx.x == 3) && (threadIdx.y == 3)) {
            slice[blockDim.y + 1][blockDim.x + 1] = psi[(ix + 1)*incx + (blockIdx.y*blockDim.y + blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 1];
        }

        // Put the yz tile for the leading x (+1) index into shared memory
        slice[ty][tz] = psi[(ix + 1)*incx + iy*incy + iz];

        __syncthreads();

        x_a1 = slice[ty][tz];

        // Center + corner points 1 ahead become edge points
        acc_a1 = cc * x_a1 +
            ecxz * (slice[(ty - 1)][(tz - 1)] +
                    slice[(ty + 1)][(tz - 1)] +
                    slice[(ty - 1)][(tz + 1)] +
                    slice[(ty + 1)][(tz + 1)]);

        // edge points 1 ahead become face points in center
        accm_a1 =  (slice[ty][(tz - 1)] +
                    slice[ty][(tz + 1)] +
                    slice[(ty - 1)][tz] +
                    slice[(ty + 1)][tz]);


        // Write back the results
        if(ix >= 1) {
            b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] = acc + 
                                                                fcx * (accm + x_b1 + x_a1) +
                                                                ecxz * (accm_b1 + accm_a1);

        }

    }                   /* end for */

}

__global__ void app_cil_fourth_kernel(const double * __restrict__ psi, 
                                                double *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const double cc,
                                                const double fcx,
                                                const double ecxy)
{

    // iz and iy map to the x and y coordinates of the thread
    // within a block
    double xm, xs, xp, x0m, x0, x0p;
    int iz = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int ix=1;

    int incx = (dimz + 2) * (dimy + 2);
    int incy = dimz + 2;
    int incxr = dimz * dimy;
    int incyr = dimz;

    int ixs, ixps;
    int iys, iyms, iyps;

    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;

    x0m = psi[iys + iz];
    x0 = psi[incx + iys + iz];
    x0p = psi[2*incx + iys + iz];
    xm = psi[iys + iz - 1] +
         psi[iyms + iz] +
         psi[iyps + iz] +
         psi[iys + iz + 1];

    xs = psi[incx + iyms + iz] +
         psi[incx + iyps + iz] +
         psi[incx + iys + (iz - 1)] + 
         psi[incx + iys + (iz + 1)];

    xp = psi[2*incx + iys + iz - 1] +
         psi[2*incx + iyms + iz] +
         psi[2*incx + iyps + iz] +
         psi[2*incx + iys + iz + 1];



    for (ix = 1; ix < dimx + 1; ix++)
    {
                ixs = ix * incx;
                ixps = (ix + 1) * incx;

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    cc * x0 +
                    fcx * (x0m + x0p + xs);

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                    ecxy * (xm  + xp +
                            psi[ixs + iyms + iz - 1] +
                            psi[ixs + iyps + iz - 1] +
                            psi[ixs + iyms + iz + 1] + 
                            psi[ixs + iyps + iz + 1]);

                x0m = x0;
                x0 = x0p;
                x0p = psi[ixps + incx + iys + iz];
                xm = xs;
                xs = xp;
                xp = psi[ixps + incx + iys + iz - 1] +
                     psi[ixps + incx + iyms + iz] +
                     psi[ixps + incx + iyps + iz] +
                     psi[ixps + incx + iys + iz + 1];

    }

}

__global__ void app_cir_fourth_kernel(const double * __restrict__ psi, 
                                      double *b, 
                                      const int dimx,
                                      const int dimy,
                                      const int dimz,
                                      const double c000,
                                      const double c100)
{

    // iz and iy map to the x and y coordinates of the thread
    // within a block
    double psi_xm, psi_x, psi_xp;
    int iz = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int ix;

    int incx = (dimz + 2) * (dimy + 2);
    int incy = dimz + 2;
    int incxr = dimz * dimy;
    int incyr = dimz;

    int ixs, ixps, ixms;
    int iys, iyms, iyps;
//    int ixs;
//    int iys, iyms, iyps;

    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;

    psi_xm = psi[iys + iz];
    psi_x = psi[incx + iys + iz];
    psi_xp = psi[2*incx + iys + iz];


    for (ix = 1; ix < dimx + 1; ix++)
    {
		ixs = ix * incx;
		ixms = (ix - 1) * incx;
		ixps = (ix + 1) * incx;


                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    
#if 1
                    c000 *  psi[ixs + iys + iz] +
                    c100 * (psi[ixs + iys + (iz - 1)] +
                            psi[ixs + iys + (iz + 1)] +
                            psi[ixms + iys + iz] +
                            psi[ixps + iys + iz] +
                            psi[ixs + iyms + iz] +
                            psi[ixs + iyps + iz]);
#else
                    c000 * psi_x +
                    c100 * (psi[ixs + iys + (iz - 1)] +
                            psi[ixs + iys + (iz + 1)] +
                            psi_xm +
                            psi_xp +
                            psi[ixs + iyms + iz] +
                            psi[ixs + iyps + iz]);

                    psi_xm = psi_x;
                    psi_x = psi_xp;
                    psi_xp = psi[(ix+2)*incx + iys + iz];
#endif        
    }                   /* end for */


}




// C wrapper functions that call the cuda kernels above
extern "C" double app_cil_fourth_gpu(const double *psi, 
                                                double *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const double gridhx,
                                                const double gridhy,
                                                const double gridhz,
                                                const double xside,
                                                const double yside,
                                                const double zside,
                                                cudaStream_t cstream)
{


  dim3 Grid, Block;
  Grid.x = 3;
  Grid.y = 3;
  Block.x = dimz/Grid.x;
  Block.y = dimy/Grid.y;
  double ihx = 1.0 / (gridhx * gridhx * xside * xside);
  double cc = (-4.0 / 3.0) * (ihx + ihx + ihx);
  double fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
  double ecxy = (1.0 / 12.0) * (ihx + ihx);
//  cudaFuncSetCacheConfig(&app_cil_fourth_kernel,cudaFuncCachePreferL1);

  app_cil_fourth_kernel1<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,    
                                                   cc,
                                                   fcx,
                                                   ecxy);

  return cc;
}

extern "C" void app_cir_fourth_gpu(const double *psi, 
                                                double *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                cudaStream_t cstream)
{

  dim3 Grid, Block;
  double c000 = 0.5;
  double c100 = 1.0 / 12.0;
  Grid.x = 2;
  Grid.y = 2;
  Block.x = dimz/Grid.x;
  Block.y = dimy/Grid.y;

//  cudaFuncSetCacheConfig(&app_cir_fourth_kernel,cudaFuncCachePreferL1);
  app_cir_fourth_kernel<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,
                                                   c000,
                                                   c100);
}
#endif
#endif
