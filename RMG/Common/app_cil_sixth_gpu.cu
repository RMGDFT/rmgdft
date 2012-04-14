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

    __shared__ double slice[PY0_GRID/2 + 4][PX0_GRID/2 + 4];
    double accp_b2, accp_b1, accp, accp_a1, accp_a2;
    double accm_b2, accm_b1, accm, accm_a1;
    double acc, acc_a1, acc_a2;

    int ix = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 2;

    // thread x-index into shared memory tile
    int tx = threadIdx.x + 2;
    // thread y-index into shared memory tile
    int ty = threadIdx.y + 2;

    int iz;
    int incx = (dimz + 4) * (dimy + 4);
    int incy = dimz + 4;
    int incxr = dimz * dimy;
    int incyr = dimz;
    int incx_slice = (PY0_GRID/2 + 4);
    int incy_slice = 1;
    
    int ixs, ixms, ixps, ixmms, ixpps;
    int iys, iyms, iyps, iymms, iypps;

    int tixs, tixms, tixps, tixmms, tixpps;
    int tiys, tiyms, tiyps, tiymms, tiypps;

    ixs = ix * incx;
    ixms = (ix - 1) * incx;
    ixps = (ix + 1) * incx;
    ixmms = (ix - 2) * incx;
    ixpps = (ix + 2) * incx;

    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;
    iymms = (iy - 2) * incy;
    iypps = (iy + 2) * incy;

    tixs = tx;
    tixms = (tx - 1);
    tixps = (tx + 1);
    tixmms = (tx - 2);
    tixpps = (tx + 2);

    tiys = ty;
    tiyms = (ty - 1);
    tiyps = (ty + 1);
    tiymms = (ty - 2);
    tiypps = (ty + 2);

    accp_b1 =
            fc2x * psi[ixs + iys] +
            tcx * psi[ixms + iys] +
            tcx * psi[ixps + iys] +
            tcx * psi[ixs + iyps] +
            tcx * psi[ixs + iyms];

    accp =
            fc2x * psi[ixs + iys + 1] +
            tcx * psi[ixms + iys + 1] +
            tcx * psi[ixps + iys + 1] +
            tcx * psi[ixs + iyps + 1] +
            tcx * psi[ixs + iyms + 1];

    accm_b1 =
            fcx * psi[ixs + iys + 1] +
            ecxz * psi[ixms + iys + 1] +
            ecxz * psi[ixps + iys + 1] +
            ecxz * psi[ixs + iyms + 1] +
            ecxz * psi[ixs + iyps + 1] +
            cor * psi[ixms + iyms + 1] +
            cor * psi[ixps + iyms + 1] +
            cor * psi[ixms + iyps + 1] +
            cor * psi[ixps + iyps + 1] +
            tcx * psi[ixpps + iys + 1] +
            tcx * psi[ixmms + iys + 1] +
            tcx * psi[ixs + iypps + 1] +
            tcx * psi[ixs + iymms + 1];

    accp_a1 =
            fc2x * psi[ixs + iys + 2] +
            tcx * psi[ixms + iys + 2] +
            tcx * psi[ixps + iys + 2] +
            tcx * psi[ixs + iyps + 2] +
            tcx * psi[ixs + iyms + 2];

    accm =
            fcx * psi[ixs + iys + 2] +
            ecxz * psi[ixms + iys + 2] +
            ecxz * psi[ixps + iys + 2] +
            ecxz * psi[ixs + iyms + 2] +
            ecxz * psi[ixs + iyps + 2] +
            cor * psi[ixms + iyms + 2] +
            cor * psi[ixps + iyms + 2] +
            cor * psi[ixms + iyps + 2] +
            cor * psi[ixps + iyps + 2] +
            tcx * psi[ixpps + iys + 2] +
            tcx * psi[ixmms + iys + 2] +
            tcx * psi[ixs + iypps + 2] +
            tcx * psi[ixs + iymms + 2];

    acc_a1 = cc * psi[ixs + iys + 2] +
            fcx * psi[ixms + iys + 2] +
            fcx * psi[ixps + iys + 2] +
            fcx * psi[ixs + iyms + 2] +
            fcx * psi[ixs + iyps + 2] +
            ecxz * psi[ixms + iyms + 2] +
            ecxz * psi[ixms + iyps + 2] +
            ecxz * psi[ixps + iyms + 2] +
            ecxz * psi[ixps + iyps + 2] +
            fc2x * psi[ixmms + iys + 2] +
            fc2x * psi[ixpps + iys + 2] +
            fc2x * psi[ixs + iymms + 2] +
            fc2x * psi[ixs + iypps + 2] +
            tcx * psi[ixps + iypps + 2] +
            tcx * psi[ixps + iymms + 2] +
            tcx * psi[ixms + iypps + 2] +
            tcx * psi[ixms + iymms + 2] +
            tcx * psi[ixpps + iyps + 2] +
            tcx * psi[ixmms + iyps + 2] +
            tcx * psi[ixpps + iyms + 2] +
            tcx * psi[ixmms + iyms + 2];

    accp_a2 =
            fc2x * psi[ixs + iys + 3] +
            tcx * psi[ixms + iys + 3] +
            tcx * psi[ixps + iys + 3] +
            tcx * psi[ixs + iyps + 3] +
            tcx * psi[ixs + iyms + 3];

    accm_a1 =
            fcx * psi[ixs + iys + 3] +
            ecxz * psi[ixms + iys + 3] +
            ecxz * psi[ixps + iys + 3] +
            ecxz * psi[ixs + iyms + 3] +
            ecxz * psi[ixs + iyps + 3] +
            cor * psi[ixms + iyms + 3] +
            cor * psi[ixps + iyms + 3] +
            cor * psi[ixms + iyps + 3] +
            cor * psi[ixps + iyps + 3] +
            tcx * psi[ixpps + iys + 3] +
            tcx * psi[ixmms + iys + 3] +
            tcx * psi[ixs + iypps + 3] +
            tcx * psi[ixs + iymms + 3];

    acc_a2 = cc * psi[ixs + iys + 3] +
            fcx * psi[ixms + iys + 3] + 
            fcx * psi[ixps + iys + 3] + 
            fcx * psi[ixs + iyms + 3] + 
            fcx * psi[ixs + iyps + 3] + 
            ecxz * psi[ixms + iyms + 3] + 
            ecxz * psi[ixms + iyps + 3] + 
            ecxz * psi[ixps + iyms + 3] + 
            ecxz * psi[ixps + iyps + 3] + 
            fc2x * psi[ixmms + iys + 3] + 
            fc2x * psi[ixpps + iys + 3] + 
            fc2x * psi[ixs + iymms + 3] + 
            fc2x * psi[ixs + iypps + 3] + 
            tcx * psi[ixps + iypps + 3] + 
            tcx * psi[ixps + iymms + 3] + 
            tcx * psi[ixms + iypps + 3] + 
            tcx * psi[ixms + iymms + 3] + 
            tcx * psi[ixpps + iyps + 3] + 
            tcx * psi[ixmms + iyps + 3] + 
            tcx * psi[ixpps + iyms + 3] + 
            tcx * psi[ixmms + iyms + 3];


    // Get acc for initial z

    for (iz = 2; iz < dimz + 2; iz++)
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
                            psi[(threadIdx.x + blockIdx.x*blockDim.x)*incx + iy*incy + iz + 2];
            slice[ty][threadIdx.x + blockDim.x + 2] = 
                            psi[(threadIdx.x + blockDim.x + 2 + blockIdx.x*blockDim.x)*incx + iy*incy + iz + 2];
        }
        if(threadIdx.y < 2) {
            slice[threadIdx.y][tx] = 
                            psi[ix*incx + (threadIdx.y + blockIdx.y*blockDim.y)*incy + iz + 2];
            slice[threadIdx.y + blockDim.y + 2][tx] = 
                            psi[ix*incx + (threadIdx.y + blockDim.y + 2 + blockIdx.y*blockDim.y)*incy + iz + 2];
        }

        if((threadIdx.x == 0) && (threadIdx.y == 0)) {
            slice[0][0] =
                            psi[(blockIdx.x*blockDim.x)*incx + (blockIdx.y*blockDim.y)*incy + iz + 2];
            slice[1][0] =
                            psi[(blockIdx.x*blockDim.x)*incx + (blockIdx.y*blockDim.y + 1)*incy + iz + 2];
            slice[0][1] = 
                            psi[(blockIdx.x*blockDim.x + 1)*incx + (blockIdx.y*blockDim.y)*incy + iz + 2];
            slice[1][1] =
                            psi[(blockIdx.x*blockDim.x + 1)*incx + (blockIdx.y*blockDim.y + 1)*incy + iz + 2];
        }
        if((threadIdx.x == (blockDim.x-1)) && (threadIdx.y == 0)) {
            slice[0][blockDim.x + 2] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 2)*incx + (blockIdx.y*blockDim.y)*incy + iz + 2];
            slice[1][blockDim.x + 3] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 3)*incx + (blockIdx.y*blockDim.y + 1)*incy + iz + 2];
            slice[1][blockDim.x + 2] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + iz + 2];
            slice[0][blockDim.x + 3] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 3)*incx + (blockIdx.y*blockDim.y)*incy + iz + 2];
        }

        if((threadIdx.x == 0) && (threadIdx.y == (blockDim.y-1))) {
            slice[blockDim.y + 2][0] =
                            psi[(blockIdx.x*blockDim.x)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + iz + 2];
            slice[blockDim.y + 3][0] =
                            psi[(blockIdx.x*blockDim.x)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + iz + 2];
            slice[blockDim.y + 2][1] =
                            psi[(blockIdx.x*blockDim.x + 1)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + iz + 2];
            slice[blockDim.y + 3][1] =
                            psi[(blockIdx.x*blockDim.x + 1)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + iz + 2];
        }

        if((threadIdx.x == (blockDim.x-1)) && (threadIdx.y == (blockDim.y-1))) {
            slice[blockDim.y + 2][blockDim.x + 2] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + iz + 2];
            slice[blockDim.y + 2][blockDim.x + 3] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 3)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + iz + 2];
            slice[blockDim.y + 3][blockDim.x + 2] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + iz + 2];
            slice[blockDim.y + 3][blockDim.x + 3] = 
                            psi[(blockIdx.x*blockDim.x + blockDim.x + 3)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + iz + 2];
        }

        // Put the xy tile for the leading z (+2) index into shared memory
        slice[ty][tx] = psi[ix*incx + iy*incy + iz + 2];

        __syncthreads();

        acc_a2 = cc * slice[tiys][tixs] +
            fcx *  (slice[tiys][tixms] + 
                    slice[tiys][tixps] +
                    slice[tiyms][tixs] + 
                    slice[tiyps][tixs]) +
            ecxz * (slice[tiyms][tixms] + 
                    slice[tiyps][tixms] +
                    slice[tiyms][tixps] + 
                    slice[tiyps][tixps]) + 
            fc2x * (slice[tiys][tixmms] + 
                    slice[tiys][tixpps] +
                    slice[tiymms][tixs] + 
                    slice[tiypps][tixs]) +
            tcx *  (slice[tiypps][tixps] + 
                    slice[tiymms][tixps] +
                    slice[tiypps][tixms] + 
                    slice[tiymms][tixms] +
                    slice[tiyps][tixpps] + 
                    slice[tiyps][tixmms] +
                    slice[tiyms][tixpps] + 
                    slice[tiyms][tixmms]);

        accm_a1 =
            fcx * slice[tiys][tixs] +
            ecxz * (slice[tiys][tixms] + 
                    slice[tiys][tixps] +
                    slice[tiyms][tixs] + 
                    slice[tiyps][tixs]) +
            cor *  (slice[tiyms][tixms] + 
                    slice[tiyms][tixps] +
                    slice[tiyps][tixms] + 
                    slice[tiyps][tixps]) +
            tcx *  (slice[tiys][tixpps] + 
                    slice[tiys][tixmms] +
                    slice[tiypps][tixs] + 
                    slice[tiymms][tixs]);

        accp_a2 = 
            fc2x * slice[tiys][tixs] +
            tcx *  (slice[tiys][tixms] + 
                    slice[tiys][tixps] + 
                    slice[tiyps][tixs] + 
                    slice[tiyms][tixs]);

        b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = acc +

                                                            accm_b2 + 
                                                            accm + 

                                                            accp_b2 + 
                                                            accp_a2;
        
    }                   /* end for */


}


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
  Block.x = dimx / 2;
  Block.y = dimy / 2;
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
