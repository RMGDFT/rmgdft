#if GPU_FD_ENABLED
#include "make_conf.h"
#include "fixed_dims.h"
#ifdef FD_XSIZE
__global__ void app_cil_fourth_f_kernel(const float * __restrict__ psi, 
                                                float *b, 
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

#endif    // end of #ifdef FD_XSIZE - left hand operator needs fixed grid but right does not

__global__ void app_cir_fourth_f_kernel(const float * __restrict__ psi, 
                                      float *b, 
                                      const int dimx,
                                      const int dimy,
                                      const int dimz,
                                      const float c000,
                                      const float c100)
{

    // iz and iy map to the x and y coordinates of the thread
    // within a block
    float psi_xm, psi_x, psi_xp;
    int iz = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int ix;

    int incx = (dimz + 2) * (dimy + 2);
    int incy = dimz + 2;
    int incxr = dimz * dimy;
    int incyr = dimz;

    int ixs;
    int iys, iyms, iyps;

    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;

    psi_xm = psi[iys + iz];
    psi_x = psi[incx + iys + iz];
    psi_xp = psi[2*incx + iys + iz];


    for (ix = 1; ix < dimx + 1; ix++)
    {
		ixs = ix * incx;

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    
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

    }                   /* end for */


}




// C wrapper functions that call the cuda kernels above
extern "C" double app_cil_fourth_f_gpu(const float *psi, 
                                                float *b, 
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
  Grid.x = 1;
  Grid.y = 3;
  Block.x = dimz/Grid.x;
  Block.y = dimy/Grid.y;
  double ihx = 1.0 / (gridhx * gridhx * xside * xside);
  double cc = (-4.0 / 3.0) * (ihx + ihx + ihx);
  double fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
  double ecxy = (1.0 / 12.0) * (ihx + ihx);
//  cudaFuncSetCacheConfig(&app_cil_fourth_f_kernel,cudaFuncCachePreferL1);

  app_cil_fourth_f_kernel<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,    
                                                   cc,
                                                   fcx,
                                                   ecxy);
  return (double)cc;
}

extern "C" void app_cir_fourth_f_gpu(const float *psi, 
                                                float *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                cudaStream_t cstream)
{

  dim3 Grid, Block;
  float c000 = 0.5;
  float c100 = 1.0 / 12.0;
  Grid.x = 1;
  Grid.y = 3;
  Block.x = dimz/Grid.x;
  Block.y = dimy/Grid.y;

//  cudaFuncSetCacheConfig(&app_cir_fourth_f_kernel,cudaFuncCachePreferL1);
  app_cir_fourth_f_kernel<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,
                                                   c000,
                                                   c100);
}
#endif
