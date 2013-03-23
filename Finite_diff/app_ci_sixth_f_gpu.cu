#include "make_conf.h"
#if GPU_ENABLED
#include "fixed_dims.h"
#ifdef FD_XSIZE


__global__ void app_cil_sixth_f_kernel(const float * __restrict__ psi, 
                                                float *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const float cc,
                                                const float fcx,
                                                const float ecxz,
                                                const float cor,
                                                const float fc2x,
                                                const float tcx)
{

    __shared__ float slice[FIXED_YDIM/2 + 4][FIXED_ZDIM/2 + 4];
    float accp_b2, accp_b1, accp, accp_a1, accp_a2;
    float accm_b2, accm_b1, accm, accm_a1;
    float acc, acc_a1, acc_a2;

    // iz and iy map to the x and y coordinates of the thread
    // within a block
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

        if((threadIdx.x == 2) && (threadIdx.y == 2)) {
            slice[0][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x];
            slice[1][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x];
            slice[0][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + 1];
            slice[1][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + 1];
        }
        if((threadIdx.x == 2) && (threadIdx.y == 3)) {
            slice[0][blockDim.x + 2] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[1][blockDim.x + 3] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
            slice[1][blockDim.x + 2] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[0][blockDim.x + 3] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
        }

       if((threadIdx.x == 3) && (threadIdx.y == 2)) {
            slice[blockDim.y + 2][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x];
            slice[blockDim.y + 3][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x];
            slice[blockDim.y + 2][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x + 1];
            slice[blockDim.y + 3][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x + 1];
        }

        if((threadIdx.x == 3) && (threadIdx.y == 3)) {
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

        acc_a2 = cc * slice[ty][tz] +
            fcx *  (slice[ty][(tz - 1)] +
                    slice[ty][(tz + 1)] +
                    slice[(ty - 1)][tz] +
                    slice[(ty + 1)][tz]) +
            ecxz * (slice[(ty - 1)][(tz - 1)] +
                    slice[(ty + 1)][(tz - 1)] +
                    slice[(ty - 1)][(tz + 1)] +
                    slice[(ty + 1)][(tz + 1)]) +
            fc2x * (slice[ty][(tz - 2)] +
                    slice[ty][(tz + 2)] +
                    slice[(ty - 2)][tz] +
                    slice[(ty + 2)][tz]) +
            tcx *  (slice[(ty + 2)][(tz + 1)] +
                    slice[(ty - 2)][(tz + 1)] +
                    slice[(ty + 2)][(tz - 1)] +
                    slice[(ty - 2)][(tz - 1)] +
                    slice[(ty + 1)][(tz + 2)] +
                    slice[(ty + 1)][(tz - 2)] +
                    slice[(ty - 1)][(tz + 2)] +
                    slice[(ty - 1)][(tz - 2)]);

        accm_a1 =
            fcx * slice[ty][tz] +
            ecxz * (slice[ty][(tz - 1)] +
                    slice[ty][(tz + 1)] +
                    slice[(ty - 1)][tz] +
                    slice[(ty + 1)][tz]) +
            cor *  (slice[(ty - 1)][(tz - 1)] +
                    slice[(ty - 1)][(tz + 1)] +
                    slice[(ty + 1)][(tz - 1)] +
                    slice[(ty + 1)][(tz + 1)]) +
            tcx *  (slice[ty][(tz + 2)] +
                    slice[ty][(tz - 2)] +
                    slice[(ty + 2)][tz] +
                    slice[(ty - 2)][tz]);

        accp_a2 =
            fc2x * slice[ty][tz] +
            tcx *  (slice[ty][(tz - 1)] +
                    slice[ty][(tz + 1)] +
                    slice[(ty + 1)][tz] +
                    slice[(ty - 1)][tz]);

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

__global__ void app_cir_sixth_f_kernel_small(const float * __restrict__ psi, 
                                                float *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const float c000,
                                                const float c100,
                                                const float c110,
                                                const float c200)
{


    // iz and iy map to the x and y coordinates of the thread
    // withing a block
    int iz = blockIdx.x * blockDim.x + threadIdx.x + 2;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 2;

    int ix=0;

    int incx = (dimz + 4) * (dimy + 4);
    int incy = dimz + 4;
    int incxr = dimz * dimy;
    int incyr = dimz;

    int ixs, ixps, ixmms, ixpps;
    int iys, iyms, iyps, iymms, iypps;
    float sum=0.0;
    float a0m, a0s, a0p, b0m, b0s, b0p;


    iys = iy * incy;
    iyms = (iy - 1) * incy;
    iyps = (iy + 1) * incy;
    iymms = (iy - 2) * incy;
    iypps = (iy + 2) * incy;

    b0m = psi[incx + iys + iz]; 
    b0s = psi[2*incx + iys + iz]; 
    b0p = psi[3*incx + iys + iz]; 

    a0m =   psi[incx + iyps + iz] +
            psi[incx + iyms + iz] +
            psi[incx + iys + (iz + 1)] +
            psi[incx + iys + (iz - 1)];
    a0s =   psi[2*incx + iyps + iz] +
            psi[2*incx + iyms + iz] +
            psi[2*incx + iys + (iz + 1)] +
            psi[2*incx + iys + (iz - 1)];
    a0p =   psi[3*incx + iyps + iz] +
            psi[3*incx + iyms + iz] +
            psi[3*incx + iys + (iz + 1)] +
            psi[3*incx + iys + (iz - 1)];

    for (ix = 2; ix < dimx + 2; ix++)
    {
		ixs = ix * incx;
		ixps = (ix + 1) * incx;
		ixmms = (ix - 2) * incx;
		ixpps = (ix + 2) * incx;

                sum = c000 * b0s +
                      c100 * (a0s + b0m + b0p);

                b0m = b0s;
                b0s = b0p;
                b0p = psi[ixps + incx + iys + iz];

                sum += c110 * (a0m + a0p +
                               psi[ixs + iyps + (iz + 1)] +
                               psi[ixs + iyps + (iz - 1)] +
                               psi[ixs + iyms + (iz + 1)] + 
                               psi[ixs + iyms + (iz - 1)]);



                a0m = a0s;
                a0s = a0p;
                a0p =   psi[ixps + incx + iyps + iz] +
                        psi[ixps + incx + iyms + iz] +
                        psi[ixps + incx + iys + (iz + 1)] +
                        psi[ixps + incx + iys + (iz - 1)];
        
                b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = sum +
                    c200 * (psi[ixs + iys + (iz - 2)] +
                            psi[ixs + iys + (iz + 2)] +
                            psi[ixmms + iys + iz] +
                            psi[ixpps + iys + iz] +
                            psi[ixs + iymms + iz] + 
                            psi[ixs + iypps + iz]);

        
    }                   /* end for */


}


__global__ void app_cir_sixth_f_kernel(const float * __restrict__ psi, 
                                                float *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const float c000,
                                                const float c100,
                                                const float c110,
                                                const float c200)
{

    __shared__ float slice[FIXED_YDIM/2 + 4][FIXED_ZDIM/2 + 4];
    float accp_b2, accp_b1, accp, accp_a1, accp_a2;
    float accm_b1, accm, accm_a1, accm_a2;
    float acc, acc_a1, acc_a2;
    float acce, acce_a1, acce_a2;

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

    accp_b1 = 0.0;
    accm_b1 = 0.0;

    accp_a1 = psi[izs + iys + ixs];

    accm_a1 = psi[izms + iys + ixs] +
              psi[izps + iys + ixs] +
              psi[izs + iyms + ixs] +
              psi[izs + iyps + ixs];

    acc_a1 =  psi[izms + iyms + ixs] + 
              psi[izms + iyps + ixs] + 
              psi[izps + iyms + ixs] + 
              psi[izps + iyps + ixs];

    acce_a1 = psi[izmms + iys + ixs] + 
              psi[izpps + iys + ixs] + 
              psi[izs + iymms + ixs] + 
              psi[izs + iypps + ixs];

    accm_a2 = psi[izms + iys + ixps] +
              psi[izps + iys + ixps] +
              psi[izs + iyms + ixps] +
              psi[izs + iyps + ixps];

    acc_a2 =  psi[izms + iyms + ixps] + 
              psi[izms + iyps + ixps] + 
              psi[izps + iyms + ixps] + 
              psi[izps + iyps + ixps];

    acce_a2 = psi[izmms + iys + ixps] + 
              psi[izpps + iys + ixps] + 
              psi[izs + iymms + ixps] + 
              psi[izs + iypps + ixps];

    accp_a2 = psi[izs + iys + ixps];

    for (ix = 0; ix < dimx + 2; ix++)
    {

        // Advance the slice partial sums
        accp_b2 = accp_b1;
        accp_b1 = accp;
        accp = accp_a1;
        accp_a1 = accp_a2;

        accm_b1 = accm;
        accm = accm_a1;
        accm_a1 = accm_a2;

        acce = acce_a1;
        acce_a1 = acce_a2;

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

        if((threadIdx.x == 2) && (threadIdx.y == 2)) {
            slice[0][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x];
            slice[1][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x];
            slice[0][1] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + 1];
            slice[1][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + 1];
        }
        if((threadIdx.x == 3) && (threadIdx.y == 2)) {
            slice[0][blockDim.x + 2] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[1][blockDim.x + 3] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
            slice[1][blockDim.x + 2] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + 1)*incy + blockIdx.x*blockDim.x + blockDim.x + 2];
            slice[0][blockDim.x + 3] = 
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y)*incy + blockIdx.x*blockDim.x + blockDim.x + 3];
        }

        if((threadIdx.x == 2) && (threadIdx.y == 3)) {
            slice[blockDim.y + 2][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x];
            slice[blockDim.y + 3][0] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x];
            slice[blockDim.y + 2][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 2)*incy + blockIdx.x*blockDim.x + 1];
            slice[blockDim.y + 3][1] =
                            psi[(ix + 2)*incx + (blockIdx.y*blockDim.y + blockDim.y + 3)*incy + blockIdx.x*blockDim.x + 1];
        }

        if((threadIdx.x == 3) && (threadIdx.y == 3)) {
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

        acce_a2 =   slice[ty][(tz - 2)] + 
                    slice[ty][(tz + 2)] + 
                    slice[(ty - 2)][tz] + 
                    slice[(ty + 2)][tz];

        acc_a2 =    slice[(ty - 1)][(tz - 1)] + 
                    slice[(ty + 1)][(tz - 1)] +
                    slice[(ty - 1)][(tz + 1)] + 
                    slice[(ty + 1)][(tz + 1)];

        // Edge points 1 ahead 
        accm_a2 =   slice[ty][(tz - 1)] + 
                    slice[ty][(tz + 1)] +
                    slice[(ty - 1)][tz] + 
                    slice[(ty + 1)][tz];

        // Central point 2 ahead
        accp_a2 = slice[ty][tz];

        // Write back the results
        if(ix >= 2) {
            b[(ix - 2) * incxr + (iy - 2) * incyr + (iz - 2)] = 

                                                            c110 * acc +

                                                            c100 * accm +
                                                            c110 * accm_a1 + 
                                                            c110 * accm_b1 + 

                                                            // c200 points in plane
                                                            c200 * acce +

                                                            // central points
                                                            c200 * accp_b2 +
                                                            c100 * accp_b1 +
                                                            c000 * accp +
                                                            c100 * accp_a1 +
                                                            c200 * accp_a2;
        }
        
    }                   /* end for */


}

#if 1
__global__ void app_cil_sixth_f_batch_kernel(const float * __restrict__ psi,
                                                float *b,
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                const float cc,
                                                const float fcx,
                                                const float ecxz,
                                                const float cor,
                                                const float fc2x,
                                                const float tcx,
                                                int nthreads)
{

    dim3 Grid, Block;
    Grid.x = 2;
    Grid.y = 2;
    Block.x = dimz/Grid.x;
    Block.y = dimy/Grid.y;

    int i, sbasis, pbasis;
    sbasis = (dimx + 4) * (dimy + 4) * (dimz + 4); 
    pbasis = dimx*dimy*dimz;
    cudaStream_t s[16];

    if(threadIdx.x == 0) {
        
        for(i = 0;i < nthreads;i++) {

            cudaStreamCreateWithFlags( &s[i], cudaStreamNonBlocking );
            app_cil_sixth_f_kernel<<<Grid, Block, 0, s[i]>>>(
                                                   &psi[i*sbasis],
                                                   &b[i*pbasis],
                                                   dimx,
                                                   dimy,
                                                   dimz,
                                                   (float)cc,
                                                   (float)fcx,
                                                   (float)ecxz,
                                                   (float)cor,
                                                   (float)fc2x,
                                                   (float)tcx);
       }
    }

    cudaDeviceSynchronize();

    if(threadIdx.x == 0) {
        for(i = 0;i < nthreads;i++) {
           cudaStreamDestroy( s[i] );
        }
    }
}

#endif

#if 1
// C wrapper functions that call the cuda kernels above
extern "C" double app_cil_sixth_f_batch_gpu(const float *psi, 
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
  Grid.y = 1;
  Block.x = 1;
  Block.y = 1;
  double ihx = 1.0 / (gridhx * gridhx * xside * xside);
  double ihy = 1.0 / (gridhy * gridhy * yside * yside);
  double ihz = 1.0 / (gridhz * gridhz * zside * zside);
  double cc = (-116.0 / 90.0) * (ihx + ihy + ihz);
  double fcx = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihx;
  double ecxz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihz;
  double cor = (1.0 / 144.0) * (ihx + ihy + ihz);
  double fc2x = (1.0 / 120.0) * (ihy + ihz);
  double tcx = (-1.0 / 240.0) * ihx;

  app_cil_sixth_f_batch_kernel<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,    
                                                   (float)cc,
                                                   (float)fcx,
                                                   (float)ecxz,
                                                   (float)cor,
                                                   (float)fc2x,
                                                   (float)tcx,
                                                   12);
  return cc;
}
#endif

// C wrapper functions that call the cuda kernels above
extern "C" double app_cil_sixth_f_gpu(const float *psi, 
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
  Grid.x = 2;
  Grid.y = 2;
  Block.x = dimz/Grid.x;
  Block.y = dimy/Grid.y;
  double ihx = 1.0 / (gridhx * gridhx * xside * xside);
  double ihy = 1.0 / (gridhy * gridhy * yside * yside);
  double ihz = 1.0 / (gridhz * gridhz * zside * zside);
  double cc = (-116.0 / 90.0) * (ihx + ihy + ihz);
  double fcx = (31.0 / 232.0) * cc + (49.0 / 60.0) * ihx;
  double ecxz = (-31.0 / 464.0) * cc - (1.0 / 10.0) * ihz;
  double cor = (1.0 / 144.0) * (ihx + ihy + ihz);
  double fc2x = (1.0 / 120.0) * (ihy + ihz);
  double tcx = (-1.0 / 240.0) * ihx;

  app_cil_sixth_f_kernel<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,    
                                                   (float)cc,
                                                   (float)fcx,
                                                   (float)ecxz,
                                                   (float)cor,
                                                   (float)fc2x,
                                                   (float)tcx);
  return cc;
}



extern "C" void app_cir_sixth_f_gpu(const float *psi, 
                                                float *b, 
                                                const int dimx,
                                                const int dimy,
                                                const int dimz,
                                                cudaStream_t cstream)
{

  dim3 Grid, Block;
  Grid.x = 2;
  Grid.y = 2;
  Block.x = dimz/Grid.x;
  Block.y = dimy/Grid.y;

  double c000 = 61.0 / 120.0;
  double c100 = 13.0 / 180.0;
  double c110 = 1.0 / 144.0;
  double c200 = -1.0 / 240.0;
//cudaFuncSetCacheConfig(&app_cir_sixth_f_kernel_small,cudaFuncCachePreferL1);

  app_cir_sixth_f_kernel<<<Grid, Block, 0, cstream>>>(
                                                   psi,
                                                   b,
                                                   dimx,    
                                                   dimy,    
                                                   dimz,
                                                   (float)c000,
                                                   (float)c100,
                                                   (float)c110,
                                                   (float)c200);
}
#endif
#endif
