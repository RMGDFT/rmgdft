#include "params.h"

__global__ void app_cil_sixth_cuda_kernel(const double *psi, 
                                                double *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const double cc,
                                                const double fcx,
                                                const double ecxz,
                                                const double cor,
                                                const double fc2x,
                                                const double tcx)
{

    __shared__ double slice[PY0_GRID/2 + 4][PZ0_GRID/2 + 4];
    double accp_b2, accp_b1, accp, accp_a1, accp_a2;
    double accm_b2, accm_b1, accm, accm_a1;
    double acc, acc_a1, acc_a2;

    // iz and iy map to the x and y coordinates of the thread
    // withing a block
    int iz = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 2;

    // thread z-index into shared memory tile
    int tz = threadIdx.x + 2;
    // thread y-index into shared memory tile
    int ty = threadIdx.y + 2;
    int ix=0;

    int incx = (dimz + 4) * (dimy + 4);
    int incy = dimz + 4;
    int incxr = dimz * dimy;
    int incyr = dimz;
    
    int ixs, ixps;
    int iys, iyms, iyps, iymms, iypps;
    int izs, izms, izps, izmms, izpps;

    int tizs, tizms, tizps, tizmms, tizpps;
    int tiys, tiyms, tiyps, tiymms, tiypps;

    ixs = ix * incx;
    ixps = (ix + 1) * incx;

    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;
    iymms = (iy - 2) * incy;
    iypps = (iy + 2) * incy;

    izs = iz;
    izms = (iz - 1);
    izps = (iz + 1);
    izmms = (iz - 2);
    izpps = (iz + 2);

    tizs = tz;
    tizms = (tz - 1);
    tizps = (tz + 1);
    tizmms = (tz - 2);
    tizpps = (tz + 2);

    tiys = ty;
    tiyms = (ty - 1);
    tiyps = (ty + 1);
    tiymms = (ty - 2);
    tiypps = (ty + 2);

    accp_b1 = 0.0;

    accp_a1 =
            fc2x * psi[izs + iys + ixs] +
            tcx * psi[izms + iys + ixs] +
            tcx * psi[izps + iys + ixs] +
            tcx * psi[izs + iyps + ixs] +
            tcx * psi[izs + iyms + ixs];


    accm =
            fcx * psi[izs + iys + ixs] +
            ecxz * psi[izms + iys + ixs] +
            ecxz * psi[izps + iys + ixs] +
            ecxz * psi[izs + iyms + ixs] +
            ecxz * psi[izs + iyps + ixs] +
            cor * psi[izms + iyms + ixs] +
            cor * psi[izps + iyms + ixs] +
            cor * psi[izms + iyps + ixs] +
            cor * psi[izps + iyps + ixs] +
            tcx * psi[izpps + iys + ixs] +
            tcx * psi[izmms + iys + ixs] +
            tcx * psi[izs + iypps + ixs] +
            tcx * psi[izs + iymms + ixs];


    acc_a1 = cc * psi[izs + iys + ixs] +
            fcx * psi[izms + iys + ixs] +
            fcx * psi[izps + iys + ixs] +
            fcx * psi[izs + iyms + ixs] +
            fcx * psi[izs + iyps + ixs] +
            ecxz * psi[izms + iyms + ixs] +
            ecxz * psi[izms + iyps + ixs] +
            ecxz * psi[izps + iyms + ixs] +
            ecxz * psi[izps + iyps + ixs] +
            fc2x * psi[izmms + iys + ixs] +
            fc2x * psi[izpps + iys + ixs] +
            fc2x * psi[izs + iymms + ixs] +
            fc2x * psi[izs + iypps + ixs] +
            tcx * psi[izps + iypps + ixs] +
            tcx * psi[izps + iymms + ixs] +
            tcx * psi[izms + iypps + ixs] +
            tcx * psi[izms + iymms + ixs] +
            tcx * psi[izpps + iyps + ixs] +
            tcx * psi[izmms + iyps + ixs] +
            tcx * psi[izpps + iyms + ixs] +
            tcx * psi[izmms + iyms + ixs];

    accp_a2 =
            fc2x * psi[izs + iys + ixps] +
            tcx * psi[izms + iys + ixps] +
            tcx * psi[izps + iys + ixps] +
            tcx * psi[izs + iyps + ixps] +
            tcx * psi[izs + iyms + ixps];

    accm_a1 =
            fcx * psi[izs + iys + ixps] +
            ecxz * psi[izms + iys + ixps] +
            ecxz * psi[izps + iys + ixps] +
            ecxz * psi[izs + iyms + ixps] +
            ecxz * psi[izs + iyps + ixps] +
            cor * psi[izms + iyms + ixps] +
            cor * psi[izps + iyms + ixps] +
            cor * psi[izms + iyps + ixps] +
            cor * psi[izps + iyps + ixps] +
            tcx * psi[izpps + iys + ixps] +
            tcx * psi[izmms + iys + ixps] +
            tcx * psi[izs + iypps + ixps] +
            tcx * psi[izs + iymms + ixps];

    acc_a2 = cc * psi[izs + iys + ixps] +
            fcx * psi[izms + iys + ixps] + 
            fcx * psi[izps + iys + ixps] + 
            fcx * psi[izs + iyms + ixps] + 
            fcx * psi[izs + iyps + ixps] + 
            ecxz * psi[izms + iyms + ixps] + 
            ecxz * psi[izms + iyps + ixps] + 
            ecxz * psi[izps + iyms + ixps] + 
            ecxz * psi[izps + iyps + ixps] + 
            fc2x * psi[izmms + iys + ixps] + 
            fc2x * psi[izpps + iys + ixps] + 
            fc2x * psi[izs + iymms + ixps] + 
            fc2x * psi[izs + iypps + ixps] + 
            tcx * psi[izps + iypps + ixps] + 
            tcx * psi[izps + iymms + ixps] + 
            tcx * psi[izms + iypps + ixps] + 
            tcx * psi[izms + iymms + ixps] + 
            tcx * psi[izpps + iyps + ixps] + 
            tcx * psi[izmms + iyps + ixps] + 
            tcx * psi[izpps + iyms + ixps] + 
            tcx * psi[izmms + iyms + ixps];


    for (ix = 0; ix < dimx + 2; ix++)
    {

        // Advance the slice partial sums
        accp_b2 = accp_b1;
        accp_b1 = accp;
        accp = accp_a1;
        accp_a1 = accp_a2;

        accm_b2 = accm_b1;
        accm_b1 = accm;
        accm = accm_a1;

        acc = acc_a1;
        acc_a1 = acc_a2;

        __syncthreads();
        // Update the data slice in shared memory
        if(threadIdx.x < 2) {
            slice[ty][threadIdx.x] = 
                            psi[(ix + 2)*incx + iy*incy + (threadIdx.x + blockIdx.x*blockDim.x)];
            slice[ty][threadIdx.x + blockDim.x + 2] = 
                            psi[(ix + 2)*incx + iy*incy + (threadIdx.x + blockDim.x + 2 + blockIdx.x*blockDim.x)];
        }
        if(threadIdx.y < 2) {
            slice[threadIdx.y][tz] = 
                            psi[(ix + 2)*incx + (threadIdx.y + blockIdx.y*blockDim.y)*incy + iz];
            slice[threadIdx.y + blockDim.y + 2][tz] = 
                            psi[(ix + 2)*incx + (threadIdx.y + blockDim.y + 2 + blockIdx.y*blockDim.y)*incy + iz];
        }

        if((threadIdx.x == 0) && (threadIdx.y == 0)) {
            slice[0][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x];
            slice[1][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x];
            slice[0][1] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + 1];
            slice[1][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + 1];
        }
        if((threadIdx.x == (blockDim.x-1)) && (threadIdx.y == 0)) {
            slice[0][blockDim.x + 2] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[1][blockDim.x + 3] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
            slice[1][blockDim.x + 2] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[0][blockDim.x + 3] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
        }

        if((threadIdx.x == 0) && (threadIdx.y == (blockDim.y-1))) {
            slice[blockDim.y + 2][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x];
            slice[blockDim.y + 3][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x];
            slice[blockDim.y + 2][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x + 1];
            slice[blockDim.y + 3][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x + 1];
        }

        if((threadIdx.x == (blockDim.x-1)) && (threadIdx.y == (blockDim.y-1))) {
            slice[blockDim.y + 2][blockDim.x + 2] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[blockDim.y + 2][blockDim.x + 3] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
            slice[blockDim.y + 3][blockDim.x + 2] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[blockDim.y + 3][blockDim.x + 3] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
        }

        // Put the xy tile for the leading x (+2) index into shared memory
        slice[ty][tz] = psi[(ix + 2)*incx + iy*incy + iz];

        __syncthreads();

        acc_a2 = cc * slice[tiys][tizs] +
            fcx *  (slice[tiys][tizms] + 
                    slice[tiys][tizps] +
                    slice[tiyms][tizs] + 
                    slice[tiyps][tizs]) +
            ecxz * (slice[tiyms][tizms] + 
                    slice[tiyps][tizms] +
                    slice[tiyms][tizps] + 
                    slice[tiyps][tizps]) + 
            fc2x * (slice[tiys][tizmms] + 
                    slice[tiys][tizpps] +
                    slice[tiymms][tizs] + 
                    slice[tiypps][tizs]) +
            tcx *  (slice[tiypps][tizps] + 
                    slice[tiymms][tizps] +
                    slice[tiypps][tizms] + 
                    slice[tiymms][tizms] +
                    slice[tiyps][tizpps] + 
                    slice[tiyps][tizmms] +
                    slice[tiyms][tizpps] + 
                    slice[tiyms][tizmms]);

        accm_a1 =
            fcx * slice[tiys][tizs] +
            ecxz * (slice[tiys][tizms] + 
                    slice[tiys][tizps] +
                    slice[tiyms][tizs] + 
                    slice[tiyps][tizs]) +
            cor *  (slice[tiyms][tizms] + 
                    slice[tiyms][tizps] +
                    slice[tiyps][tizms] + 
                    slice[tiyps][tizps]) +
            tcx *  (slice[tiys][tizpps] + 
                    slice[tiys][tizmms] +
                    slice[tiypps][tizs] + 
                    slice[tiymms][tizs]);

        accp_a2 = 
            fc2x * slice[tiys][tizs] +
            tcx *  (slice[tiys][tizms] + 
                    slice[tiys][tizps] + 
                    slice[tiyps][tizs] + 
                    slice[tiyms][tizs]);

        // Write back the results
        if(ix >= 2) {
            b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = acc +

                                                            accm_b2 + 
                                                            accm + 

                                                            accp_b2 + 
                                                            accp_a2;
        }
        
    }                   /* end for */


}


// This is the C wrapper function that calls the cuda kernel above
extern "C" void app_cil_sixth_gpu(const double *psi, 
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
  Grid.x = 2;
  Grid.y = 2;
  Block.x = dimz/2;
  Block.y = dimy/2;
  double ihx = 1.0 / (gridhx * gridhx * xside * xside);
  double ihy = 1.0 / (gridhy * gridhy * yside * yside);
  double ihz = 1.0 / (gridhz * gridhz * zside * zside);
  double cc = (-116.0 / 90.0) * (ihx + ihy + ihz);
  double fcx = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihx;
  double ecxz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihz;
  double cor = (1.0 / 144.0) * (ihx + ihy + ihz);
  double fc2x = (1.0 / 120.0) * (ihy + ihz);
  double tcx = (-1.0 / 240.0) * ihx;

  app_cil_sixth_cuda_kernel<<<Grid, Block, 0, cstream>>>(
//  app_cil_sixth_cuda_kernel<<<Grid, Block>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,    
                                                   cc,
                                                   fcx,
                                                   ecxz,
                                                   cor,
                                                   fc2x,
                                                   tcx);
}
