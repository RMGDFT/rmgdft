#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "common_prototypes.h"
#include "main.h"

/*Conversion factor to convert dipole in a.u. to debye*/
#define DEBYE_CONVERSION 2.54174618772314

void get_dipole (double * rho, double *rhoc, double *dipole)
{

    double hxgrid, hygrid, hzgrid;
    double xc, yc, zc, ax[3], bx[3], x, y, z, temp, icharge, vel;
    int ix, iy, iz, ion, FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET, dimx, dimy, dimz ;
    ION *iptr;


    for(int i = 0; i < 3; i++) dipole[i] = 0.0;
    vel = get_vel_f();

    hxgrid = get_hxxgrid();
    hygrid = get_hyygrid();
    hzgrid = get_hzzgrid();

    dimx = get_FPX0_GRID();
    dimy = get_FPY0_GRID();
    dimz = get_FPZ0_GRID();

    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();

    double dipole_center_xtal[3];
    to_crystal(dipole_center_xtal, ct.dipole_center);

    /*This gets dipole moment for electrons */

    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {

                /*vector between a grid point and a middle of the cell in crystal coordinates */
                int ixx = ix + FPX_OFFSET;
                int iyy = iy + FPY_OFFSET;
                int izz = iz + FPZ_OFFSET;
                xc = ixx * hxgrid - dipole_center_xtal[0] +0.5;
                yc = iyy * hygrid - dipole_center_xtal[1] +0.5;
                zc = izz * hzgrid - dipole_center_xtal[2] +0.5;

                if(xc > 1.0) xc -= 1.0;
                if(yc > 1.0) yc -= 1.0;
                if(zc > 1.0) zc -= 1.0;
                if(xc < 0.0) xc += 1.0;
                if(yc < 0.0) yc += 1.0;
                if(zc < -.0) zc += 1.0;
                ax[0] = xc ;
                ax[1] = yc ;
                ax[2] = zc ;

                /*Transform ax into vector in cartesian coordinates */
                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                temp = rho[ix * dimy * dimz + iy * dimz + iz];
                temp -= rhoc[ix * dimy * dimz + iy * dimz + iz];

                dipole[0] += x * temp;
                dipole[1] += y * temp;
                dipole[2] += z * temp;


            }                   /* end for */


        }                       /* end for */


    }                           /* end for */





    /* Sum these up over all processors, multiply by volume elemnt so that sum is 
     * integration and invert sign, since electron charge should be negative*/
    dipole[0] = -1.0 * vel * real_sum_all (dipole[0], pct.grid_comm);
    dipole[1] = -1.0 * vel * real_sum_all (dipole[1], pct.grid_comm);
    dipole[2] = -1.0 * vel * real_sum_all (dipole[2], pct.grid_comm);


    /*Now we have dipole moment for electrons, need to add ions now */
    if(!ct.localize_localpp)
    {
        for (ion = 0; ion < ct.num_ions; ion++)
        {

            /*Ionic pointer and ionic charge */
            iptr = &Atoms[ion];
            icharge = Species[iptr->species].zvalence;

            /*Difference vector between center of the cell and ionic position */
            ax[0] = iptr->xtal[0] - dipole_center_xtal[0] +0.5;
            ax[1] = iptr->xtal[1] - dipole_center_xtal[1] +0.5;
            ax[2] = iptr->xtal[2] - dipole_center_xtal[2] +0.5;

            if(ax[0] > 1.0) ax[0] -= 1.0;
            if(ax[1] > 1.0) ax[1] -= 1.0;
            if(ax[2] > 1.0) ax[2] -= 1.0;
            if(ax[0] < 0.0) ax[0] += 1.0;
            if(ax[1] < 0.0) ax[1] += 1.0;
            if(ax[2] < 0.0) ax[2] += 1.0;

            /*Find vector ax in cartesian coordinates */
            to_cartesian (ax, bx);

            /*Ionic contribution to dipole moment */
            dipole[0] += icharge * bx[0];
            dipole[1] += icharge * bx[1];
            dipole[2] += icharge * bx[2];

        }                           /*end for (ion=0; ion<c->num_ions; ion++) */
    }

    /*Now we need to convert to debye units */
    //if (pct.imgpe==0)
    //{
    //   printf("\n\n Dipole moment [Debye]: Absolute value: %.3f, vector: (%.3f,%.3f, %.3f)", 
    //          DEBYE_CONVERSION * sqrt (px * px + py * py + pz * pz), 
    //         DEBYE_CONVERSION *px, 
    //        DEBYE_CONVERSION *py, 
    //       DEBYE_CONVERSION *pz);
    //}




}                               /* end getpoi_bc.c */

/******/
void get_dipole (double * rho, double *dipole)
{

    double hxgrid, hygrid, hzgrid;
    double xc, yc, zc, ax[3], bx[3], x, y, z, temp, icharge, vel;
    int ix, iy, iz, ion, FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET, dimx, dimy, dimz ;
    ION *iptr;


    for(int i = 0; i < 3; i++) dipole[i] = 0.0;
    vel = get_vel_f();

    hxgrid = get_hxxgrid();
    hygrid = get_hyygrid();
    hzgrid = get_hzzgrid();

    dimx = get_FPX0_GRID();
    dimy = get_FPY0_GRID();
    dimz = get_FPZ0_GRID();

    FPX_OFFSET = get_FPX_OFFSET();
    FPY_OFFSET = get_FPY_OFFSET();
    FPZ_OFFSET = get_FPZ_OFFSET();

    double dipole_center_xtal[3];
    to_crystal(dipole_center_xtal, ct.dipole_center);

    /*This gets dipole moment for electrons */

    for (ix = 0; ix < dimx; ix++)
    {

        for (iy = 0; iy < dimy; iy++)
        {

            for (iz = 0; iz < dimz; iz++)
            {

                /*vector between a grid point and a middle of the cell in crystal coordinates */
                int ixx = ix + FPX_OFFSET;
                int iyy = iy + FPY_OFFSET;
                int izz = iz + FPZ_OFFSET;
                xc = ixx * hxgrid - dipole_center_xtal[0] +0.5;
                yc = iyy * hygrid - dipole_center_xtal[1] +0.5;
                zc = izz * hzgrid - dipole_center_xtal[2] +0.5;

                if(xc > 1.0) xc -= 1.0;
                if(yc > 1.0) yc -= 1.0;
                if(zc > 1.0) zc -= 1.0;
                if(xc < 0.0) xc += 1.0;
                if(yc < 0.0) yc += 1.0;
                if(zc < -.0) zc += 1.0;
                ax[0] = xc ;
                ax[1] = yc ;
                ax[2] = zc ;

                /*Transform ax into vector in cartesian coordinates */
                to_cartesian (ax, bx);
                x = bx[0];
                y = bx[1];
                z = bx[2];

                temp = rho[ix * dimy * dimz + iy * dimz + iz];

                dipole[0] += x * temp;
                dipole[1] += y * temp;
                dipole[2] += z * temp;


            }                   /* end for */


        }                       /* end for */


    }                           /* end for */





    /* Sum these up over all processors, multiply by volume elemnt so that sum is 
     * integration and invert sign, since electron charge should be negative*/
    dipole[0] = -1.0 * vel * real_sum_all (dipole[0], pct.grid_comm);
    dipole[1] = -1.0 * vel * real_sum_all (dipole[1], pct.grid_comm);
    dipole[2] = -1.0 * vel * real_sum_all (dipole[2], pct.grid_comm);

    rmg_printf("\n electronic dipole %f %f %f", dipole[0], dipole[1], dipole[2]);

    /*Now we have dipole moment for electrons, need to add ions now */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Ionic pointer and ionic charge */
        iptr = &Atoms[ion];
        icharge = Species[iptr->species].zvalence;

        /*Difference vector between center of the cell and ionic position */
        ax[0] = iptr->xtal[0] - dipole_center_xtal[0] +0.5;
        ax[1] = iptr->xtal[1] - dipole_center_xtal[1] +0.5;
        ax[2] = iptr->xtal[2] - dipole_center_xtal[2] +0.5;

        if(ax[0] > 1.0) ax[0] -= 1.0;
        if(ax[1] > 1.0) ax[1] -= 1.0;
        if(ax[2] > 1.0) ax[2] -= 1.0;
        if(ax[0] < 0.0) ax[0] += 1.0;
        if(ax[1] < 0.0) ax[1] += 1.0;
        if(ax[2] < 0.0) ax[2] += 1.0;

        /*Find vector ax in cartesian coordinates */
        to_cartesian (ax, bx);

        /*Ionic contribution to dipole moment */
        dipole[0] += icharge * bx[0];
        dipole[1] += icharge * bx[1];
        dipole[2] += icharge * bx[2];

    }                           /*end for (ion=0; ion<c->num_ions; ion++) */

    /*Now we need to convert to debye units */
    //if (pct.imgpe==0)
    //{
    //   printf("\n\n Dipole moment [Debye]: Absolute value: %.3f, vector: (%.3f,%.3f, %.3f)", 
    //          DEBYE_CONVERSION * sqrt (px * px + py * py + pz * pz), 
    //         DEBYE_CONVERSION *px, 
    //        DEBYE_CONVERSION *py, 
    //       DEBYE_CONVERSION *pz);
    //}




}                               /* end getpoi_bc.c */

/******/

