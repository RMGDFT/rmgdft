/************************** SVN Revision Information **************************
 **    $Id: InitNucForward.c 3542 2016-05-03 18:37:56Z luw $    **
******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Atomic.h"

#include "AtomicInterpolate.h"

void InitLocalBackward (double * vnuc_f, double * rhoc_f, double * rhocore_f)
{

    int ix, iy, iz, ixx, iyy, izz;
    int xstart, ystart, zstart, xend, yend, zend;
    int ion, idx, idx1;
    int ix1, iy1, iz1;
    int ilow, jlow, klow, ihi, jhi, khi;
    int dimx, dimy, dimz;
    int FP0_BASIS;
    int FPX0_GRID, FPY0_GRID, FPZ0_GRID;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET;
    int FNX_GRID, FNY_GRID, FNZ_GRID;
    int coarse_size;

    fftw_complex *in, *out;
    fftw_plan p2;


    double r, rc, rcnorm, t1;
    double hxxgrid, hyygrid, hzzgrid;
    SPECIES *sp;
    ION *iptr;

    std::complex<double> *gbptr, *vnuc_ptr, *rhoc_ptr, *rhocore_ptr, *phase_fftw;
    std::complex<double> *fptr;
    double shift[3], crds_shift[3];
   
    hxxgrid = get_hxxgrid();
    hyygrid = get_hyygrid();
    hzzgrid = get_hzzgrid();

    FP0_BASIS = get_FP0_BASIS();
    FPX0_GRID = get_FPX0_GRID();
    FPY0_GRID = get_FPY0_GRID();
    FPZ0_GRID = get_FPZ0_GRID();
    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();
    FNX_GRID = get_FNX_GRID();
    FNY_GRID = get_FNY_GRID();
    FNZ_GRID = get_FNZ_GRID();


    ilow = FPX_OFFSET;
    jlow = FPY_OFFSET;
    klow = FPZ_OFFSET;
    ihi = ilow + FPX0_GRID;
    jhi = jlow + FPY0_GRID;
    khi = klow + FPZ0_GRID;

    /* Loop over ions determine num of ions whose local pp has overlap with this pe */

    pct.num_loc_ions = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        dimx =  sp->ldim;
        dimy =  sp->ldim;
        dimz =  sp->ldim;

        dimx = dimx * 2 + 1;
        dimy = dimy * 2 + 1;
        dimz = dimz * 2 + 1;
        

        xstart = iptr->xtal[0] / hxxgrid - dimx/2;
        xend = xstart + dimx;
        ystart = iptr->xtal[1] / hyygrid - dimy/2;
        yend = ystart + dimy;
        zstart = iptr->xtal[2] / hzzgrid - dimz/2;
        zend = zstart + dimz;


        bool map_x = false, map_y = false, map_z = false;
        for (ix = xstart; ix < xend; ix++)
        {
            // fold the grid into the unit cell
            // maxium fold 20 times, local potential can extend to 20-1 unit cell
            ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
            if(ixx >= ilow && ixx < ihi)
            {
                map_x = true;
                break;
            }
        }

        if(!map_x) continue;

        for (iy = ystart; iy < yend; iy++)
        {
            // fold the grid into the unit cell
            iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
            if(iyy >= jlow && iyy < jhi)
            {
                map_y = true;
                break;
            }
        }

        if(!map_y) continue;


        for (iz = zstart; iz < zend; iz++)
        {
            // fold the grid into the unit cell
            izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
            if(izz >= klow && izz < khi)
            {
                map_z = true;
                break;
            }
        }

        if(!map_z) continue;

        pct.loc_ions_list[pct.num_loc_ions] = ion;
        pct.num_loc_ions ++;
    }

    if(pct.localpp != NULL) 
    {
        free(pct.localpp);
        free(pct.localrhoc);
        free(pct.localrhonlcc);
    }

    pct.localpp = (double *) malloc(pct.num_loc_ions * FP0_BASIS * sizeof(double) + 1024);
    pct.localrhoc = (double *) malloc(pct.num_loc_ions * FP0_BASIS * sizeof(double) + 1024);
    pct.localrhonlcc = (double *) malloc(pct.num_loc_ions * FP0_BASIS * sizeof(double) + 1024);
    for(idx = 0; idx < pct.num_loc_ions *FP0_BASIS; idx++)
    {
        pct.localpp[idx] = 0.0;
        pct.localrhoc[idx] = 0.0;
        pct.localrhonlcc[idx] = 0.0;
    }


    int ion1;
    for (ion1 = 0; ion1 < pct.num_loc_ions; ion1++)
    {
        ion = pct.loc_ions_list[ion1];
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];
        rc = sp->rc;
        rcnorm = rc * rc * rc * pow (PI, 1.5);
        rcnorm = 1.0 / rcnorm;

        dimx = sp->ldim;
        dimy = sp->ldim;
        dimz = sp->ldim;


        xstart = iptr->xtal[0] / hxxgrid - dimx/2;
        xend = xstart + dimx;
        ystart = iptr->xtal[1] / hyygrid - dimy/2;
        yend = ystart + dimy;
        zstart = iptr->xtal[2] / hzzgrid - dimz/2;
        zend = zstart + dimz;


        coarse_size = dimx * dimy * dimz;

        in = (fftw_complex *)fftw_malloc(sizeof(std::complex<double>) * dimx * dimy * dimz);
        out = (fftw_complex *)fftw_malloc(sizeof(std::complex<double>) * dimx * dimy * dimz);


        gbptr = new std::complex<double>[coarse_size];
        vnuc_ptr = new std::complex<double>[coarse_size];
        rhoc_ptr = new std::complex<double>[coarse_size];
        rhocore_ptr = new std::complex<double>[coarse_size];
        phase_fftw = new std::complex<double>[coarse_size];

        p2 = fftw_plan_dft_3d (dimx, dimy, dimz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        shift[0] = iptr->xtal[0] - (xstart + dimx/2) * hxxgrid;
        shift[1] = iptr->xtal[1] - (ystart + dimy/2) * hyygrid;
        shift[2] = iptr->xtal[2] - (zstart + dimz/2) * hzzgrid;
        to_cartesian(shift, crds_shift);

        RmgTimer *RT = new RmgTimer("phase");
        FindFftwPhaseLocalpp (sp->ldim, sp->ldim, sp->ldim, crds_shift, phase_fftw, 2);
        delete RT;

        fptr = (std::complex<double> *) sp->forward_vnuc;
        for(idx = 0; idx < coarse_size; idx++) 
            gbptr[idx] = fptr[idx] * std::conj(phase_fftw[idx]);
        fftw_execute_dft (p2,  reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(vnuc_ptr));


        fptr = (std::complex<double> *) sp->forward_rhoc;
        for(idx = 0; idx < coarse_size; idx++) 
            gbptr[idx] = fptr[idx] * std::conj(phase_fftw[idx]);
        fftw_execute_dft (p2,  reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(rhoc_ptr));

        if (sp->nlccflag)
        {
            fptr = (std::complex<double> *) sp->forward_rhocore;
            for(idx = 0; idx < coarse_size; idx++) 
                gbptr[idx] = fptr[idx] * std::conj(phase_fftw[idx]);
            fftw_execute_dft (p2,  reinterpret_cast<fftw_complex*>(gbptr), reinterpret_cast<fftw_complex*>(rhocore_ptr));

        }

        for (ix = xstart; ix < xend; ix++)
        {
            // fold the grid into the unit cell
            // maxium fold 20 times, local potential can extend to 20-1 unit cell
            ix1 = ix - xstart - dimx/2;
            if(ix1 < 0) ix1 += dimx;
            ixx = (ix + 20 * FNX_GRID) % FNX_GRID;
            if(ixx >= ilow && ixx < ihi)
            {

                for (iy = ystart; iy < yend; iy++)
                {
                    // fold the grid into the unit cell
                    iy1 = iy - ystart - dimy/2;
                    if(iy1 < 0) iy1 += dimy;
                    iyy = (iy + 20 * FNY_GRID) % FNY_GRID;
                    if(iyy >= jlow && iyy < jhi)
                    {
                        for (iz = zstart; iz < zend; iz++)
                        {
                            // fold the grid into the unit cell
                            iz1 = iz - zstart - dimz/2;
                            if(iz1 < 0) iz1 += dimz;
                            izz = (iz + 20 * FNZ_GRID) % FNZ_GRID;
                            if(izz >= klow && izz < khi)
                            {

                                idx = (ixx-ilow) * FPY0_GRID * FPZ0_GRID + (iyy-jlow) * FPZ0_GRID + izz-klow;
                                idx1 = ix1 * dimy * dimz + iy1 * dimz + iz1;

                                pct.localpp[ion1 * FP0_BASIS + idx] += std::real(vnuc_ptr[idx1]);
                                pct.localrhoc[ion1 * FP0_BASIS + idx] += std::real(rhoc_ptr[idx1]);


                                if (sp->nlccflag)
                                {

                                    t1 = AtomicInterpolateInline (&sp->rhocorelig[0], r);
                                    pct.localrhonlcc[ion1 * FP0_BASIS + idx] += std::real(rhocore_ptr[idx1]);

                                }

                            }                           /* end for */

                        }
                    }
                }
            }

        }
        delete [] gbptr;
        delete [] vnuc_ptr;
        delete [] rhoc_ptr;
        delete [] rhocore_ptr;
        delete [] phase_fftw;
    }

    /* Check compensating charges */
    ct.crho = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
        ct.crho += rhoc_f[idx];


    ct.crho = ct.crho * get_vel_f();
    ct.crho = real_sum_all (ct.crho, pct.grid_comm);  /* sum over pct.grid_comm  */

    if (ct.verbose)
    {
	if (pct.imgpe==0)
	    printf("\nCompensating charge is %.8e\n", ct.crho);
    }

    t1 = 0.0;
    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        if (rhocore_f[idx] < 0.0)
            rhocore_f[idx] = 0.0;
        t1 += rhocore_f[idx];
    }


    /* Wait until everyone gets here */
    /*my_barrier(); */

    /*   add the saw-tooth potential induced by applied electric field  */
    if(ct.runflag == 5 || ct.runflag == 6 || ct.forceflag == TDDFT) return;
    init_efield (vnuc_f);


}                               /* end init_nuc */

/******/
