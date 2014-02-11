/************************** SVN Revision Information **************************
 **    $Id: app_cil_fourth.c 1722 2012-05-07 19:38:37Z ebriggs $    **
******************************************************************************/

#include "const.h"
#include "common_prototypes.h"
#include "rmg_alloc.h"
#include "fixed_dims.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "hybrid.h"

static rmg_double_t app_cil_fourth_global (rmg_double_t * a, rmg_double_t * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
static rmg_double_t app_cil_fourth_standard (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);


rmg_double_t app_cil_fourth (rmg_double_t * a, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int  numgrid, tid, used_alloc=FALSE, ibrav, P0_BASIS;
    rmg_double_t cc;
    rmg_double_t *rptr=NULL;
    rmg_double_t *gpu_a, *gpu_b;

    ibrav = get_ibrav_type();
    P0_BASIS = get_P0_BASIS();

    int pbasis = dimx * dimy * dimz, itid;
    int sbasis = (dimx + 2) * (dimy + 2) * (dimz + 2);
#if GPU_FD_ENABLED
    cudaStream_t *cstream;
    cstream = get_thread_cstream();
#endif

#if HYBRID_MODEL
    tid = get_thread_tid();
    if(tid < 0) tid = 0;  // OK in this case
#else
    tid = 0;
#endif


#if (GPU_FD_ENABLED && FD_XSIZE)
    // cudaMallocHost is painfully slow so we use a pointers into regions that were previously allocated.
    rptr = (rmg_float_t *)&ct.gpu_host_fdbuf1[0];
    rptr += tid*sbasis;
    gpu_a = (rmg_double_t *)&ct.gpu_work1[0];
    gpu_a += tid*sbasis;
    gpu_b = (rmg_double_t *)&ct.gpu_work2[0];
    gpu_b += tid*pbasis;
#endif

    // If rptr is null then we must allocate it here
    if(rptr == NULL) {
        my_malloc (rptr, sbasis + 64, rmg_double_t);
        used_alloc = TRUE;
    }


    if((ibrav != CUBIC_PRIMITIVE) && (ibrav != ORTHORHOMBIC_PRIMITIVE)) {
        error_handler("Grid symmetry not programmed yet in app_cil_fourth.\n");
    }

#if (GPU_FD_ENABLED && FD_XSIZE)
    trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_FD);
    cudaMemcpyAsync( gpu_a, rptr, sbasis * sizeof(rmg_double_t), cudaMemcpyHostToDevice, *cstream);
    cc = app_cil_fourth_gpu (gpu_a, gpu_b, dimx, dimy, dimz,
                              gridhx, gridhy, gridhz,
                              get_xside(), get_yside(), get_zside(), *cstream);
    cudaMemcpyAsync(b, gpu_b, pbasis * sizeof(rmg_double_t), cudaMemcpyDeviceToHost, *cstream);

    return cc;
#endif

    trade_imagesx (a, rptr, dimx, dimy, dimz, 1, FULL_FD);


    // first check for fixed dim case 
    numgrid = dimx * dimy * dimz;
    if(numgrid == P0_BASIS && get_anisotropy() < 1.000001)
    {
        cc = app_cil_fourth_global (rptr, b, gridhx, gridhy, gridhz);
    }
    else {
        cc = app_cil_fourth_standard (rptr, b, dimx, dimy, dimz, gridhx, gridhy, gridhz);
    }


    if(used_alloc)
        my_free(rptr);
    return cc;

}

rmg_double_t app_cil_fourth_standard (rmg_double_t * rptr, rmg_double_t * b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int iz, ix, iy, incx, incy, incxr, incyr;
    int ixs, iys, ixms, ixps, iyms, iyps;
    rmg_double_t ecxy, ecxz, ecyz, cc = 0.0, fcx, fcy, fcz;
    rmg_double_t ihx, ihy, ihz, a1, a2, a3;

    incx = (dimz + 2) * (dimy + 2);
    incy = dimz + 2;
    incxr = dimz * dimy;
    incyr = dimz;

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
    ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

    if (get_anisotropy() < 1.000001)
    {

        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
        cc = (-4.0 / 3.0) * (ihx + ihx + ihx);
        fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
        ecxy = (1.0 / 12.0) * (ihx + ihx);
        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;

        for (ix = 1; ix <= dimx; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            for (iy = 1; iy <= dimy; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;

                for (iz = 1; iz <= dimz; iz++)
                {

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        cc * rptr[ixs + iys + iz] +
                        fcx * (rptr[ixms + iys + iz] +
                                rptr[ixps + iys + iz] +
                                rptr[ixs + iyms + iz] +
                                rptr[ixs + iyps + iz] +
                                rptr[ixs + iys + (iz - 1)] + rptr[ixs + iys + (iz + 1)]);

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        ecxy * (rptr[ixms + iys + iz - 1] +
                                rptr[ixps + iys + iz - 1] +
                                rptr[ixs + iyms + iz - 1] +
                                rptr[ixs + iyps + iz - 1] +
                                rptr[ixms + iyms + iz] +
                                rptr[ixms + iyps + iz] +
                                rptr[ixps + iyms + iz] +
                                rptr[ixps + iyps + iz] +
                                rptr[ixms + iys + iz + 1] +
                                rptr[ixps + iys + iz + 1] +
                                rptr[ixs + iyms + iz + 1] + rptr[ixs + iyps + iz + 1]);


                }           /* end for */

            }               /* end for */

        }                   /* end for */
    }
    else
    {

        /* Compute coefficients for this grid spacing */
        ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
        ihy = 1.0 / (gridhy * gridhy * get_yside() * get_yside());
        ihz = 1.0 / (gridhz * gridhz * get_zside() * get_zside());

        cc = (-4.0 / 3.0) * (ihx + ihy + ihz);

        fcx = (5.0 / 6.0) * ihx + (cc / 8.0);
        fcy = (5.0 / 6.0) * ihy + (cc / 8.0);
        fcz = (5.0 / 6.0) * ihz + (cc / 8.0);

        ecxy = (1.0 / 12.0) * (ihx + ihy);
        ecxz = (1.0 / 12.0) * (ihx + ihz);
        ecyz = (1.0 / 12.0) * (ihy + ihz);


        incy = dimz + 2;
        incx = (dimz + 2) * (dimy + 2);
        incyr = dimz;
        incxr = dimz * dimy;



        for (ix = 1; ix <= dimx; ix++)
        {
            ixs = ix * incx;
            ixms = (ix - 1) * incx;
            ixps = (ix + 1) * incx;
            for (iy = 1; iy <= dimy; iy++)
            {
                iys = iy * incy;
                iyms = (iy - 1) * incy;
                iyps = (iy + 1) * incy;

                for (iz = 1; iz <= dimz; iz++)
                {

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                        cc * rptr[ixs + iys + iz] +
                        fcx * rptr[ixms + iys + iz] +
                        fcx * rptr[ixps + iys + iz] +
                        fcy * rptr[ixs + iyms + iz] +
                        fcy * rptr[ixs + iyps + iz] +
                        fcz * rptr[ixs + iys + (iz - 1)] + fcz * rptr[ixs + iys + (iz + 1)];

                    b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                        ecxz * rptr[ixms + iys + iz - 1] +
                        ecxz * rptr[ixps + iys + iz - 1] +
                        ecyz * rptr[ixs + iyms + iz - 1] +
                        ecyz * rptr[ixs + iyps + iz - 1] +
                        ecxy * rptr[ixms + iyms + iz] +
                        ecxy * rptr[ixms + iyps + iz] +
                        ecxy * rptr[ixps + iyms + iz] +
                        ecxy * rptr[ixps + iyps + iz] +
                        ecxz * rptr[ixms + iys + iz + 1] +
                        ecxz * rptr[ixps + iys + iz + 1] +
                        ecyz * rptr[ixs + iyms + iz + 1] + ecyz * rptr[ixs + iyps + iz + 1];


                }           /* end for */

            }               /* end for */

        }                   /* end for */

    }                       /* end if */


    return cc;

}                               /* end app_cil */



rmg_double_t app_cil_fourth_global (rmg_double_t * rptr, rmg_double_t * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz)
{

    int ix, iy, iz;
    int ixs, iys, ixms, ixps, iyms, iyps;
    int incy, incx;
    int incyr, incxr;
    rmg_double_t rz, rzps, rzms, rzpps;
    rmg_double_t c000, c100, c110;
    rmg_double_t ihx;

    incx = (FIXED_ZDIM + 2) * (FIXED_YDIM + 2);
    incy = FIXED_ZDIM + 2;
    incxr = FIXED_ZDIM * FIXED_YDIM;
    incyr = FIXED_ZDIM;

    ihx = 1.0 / (gridhx * gridhx * get_xside() * get_xside());
    c000 = (-4.0 / 3.0) * (ihx + ihx + ihx);
    c100 = (5.0 / 6.0) * ihx + (c000 / 8.0);
    c110 = (1.0 / 12.0) * (ihx + ihx);


    // Handle the general case first
    for (ix = 1; ix < FIXED_XDIM + 1; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy < FIXED_YDIM + 1; iy++)
        {
            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz < FIXED_ZDIM + 1; iz++)
            {

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] =
                    c100 * (rptr[ixs + iys + (iz - 1)] +
                            rptr[ixs + iys + (iz + 1)] +
                            rptr[ixms + iys + iz] +
                            rptr[ixps + iys + iz] +
                            rptr[ixs + iyms + iz] +
                            rptr[ixs + iyps + iz]) + c000 * rptr[ixs + iys + iz];

                b[(ix - 1) * incxr + (iy - 1) * incyr + (iz - 1)] +=
                    c110 * (rptr[ixps + iyps + iz] +
                            rptr[ixps + iyms + iz] +
                            rptr[ixms + iyps + iz] +
                            rptr[ixms + iyms + iz] +
                            rptr[ixps + iys + (iz + 1)] +
                            rptr[ixps + iys + (iz - 1)] +
                            rptr[ixms + iys + (iz + 1)] +
                            rptr[ixms + iys + (iz - 1)] +
                            rptr[ixs + iyps + (iz + 1)] +
                            rptr[ixs + iyps + (iz - 1)] +
                            rptr[ixs + iyms + (iz + 1)] + rptr[ixs + iyms + (iz - 1)]);

            }                   /* end for */

        }                       /* end for */

    }                           /* end for */

    return c000;

}
