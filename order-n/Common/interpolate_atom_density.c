/************************** SVN Revision Information **************************
 **    $Id: interpolate_atom_density.c 648 2006-10-19 17:35:46Z luw $    **
******************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md.h"


void interpolate_atom_density(double *rho_tem, double *rho_out, int ixmin, int ixmax,
                              int iymin, int iymax, int izmin, int izmax,
                              double *crds, double *crds1, double hxgrid,
                              double hygrid, double hzgrid,
                              double hxgrid_new, double hygrid_new, double hzgrid_new)
{

    int grid_x;
    int grid_y;
    int grid_z;
    int i, ix, iy, iz , idx;
    double *data_x, *data_y, *data_z;
    double *data_x_new, *data_y_new, *data_z_new;
    double *frac_x, *frac_y, *frac_z;
    double frac000, frac100, frac010, frac001;
    double frac110, frac101, frac011, frac111;
    int idx000, idx100, idx010, idx001;
    int idx110, idx101, idx011, idx111;
    int *indx, *indy, *indz;
    int incx, incy;

    grid_x = ixmax - ixmin;
    grid_y = iymax - iymin;
    grid_z = izmax - izmin;

    incx = grid_y * grid_z;
    incy = grid_z;

    my_malloc( data_x, grid_x, double );
    my_malloc( data_y, grid_y, double );
    my_malloc( data_z, grid_z, double );
    my_malloc( data_x_new, grid_x, double );
    my_malloc( data_y_new, grid_y, double );
    my_malloc( data_z_new, grid_z, double );

    my_malloc( frac_x, grid_x, double );
    my_malloc( frac_y, grid_y, double );
    my_malloc( frac_z, grid_z, double );

    my_malloc( indx, grid_x, int );
    my_malloc( indy, grid_y, int );
    my_malloc( indz, grid_z, int );

    for (i = 0; i < grid_x; i++)
        data_x[i] = (ixmin + i) * hxgrid - crds[0];
    for (i = 0; i < grid_y; i++)
        data_y[i] = (iymin + i) * hygrid - crds[1];
    for (i = 0; i < grid_z; i++)
        data_z[i] = (izmin + i) * hzgrid - crds[2];



    /* Determine the new ixmin, iymin, izmin  */

    ixmin += (int) (crds1[0] / hxgrid_new) - (int) (crds[0] / hxgrid);
    iymin += (int) (crds1[1] / hygrid_new) - (int) (crds[1] / hygrid);
    izmin += (int) (crds1[2] / hzgrid_new) - (int) (crds[2] / hzgrid);

    for (i = 0; i < grid_x; i++)
        data_x_new[i] = (ixmin + i) * hxgrid_new - crds1[0];
    for (i = 0; i < grid_y; i++)
        data_y_new[i] = (iymin + i) * hygrid_new - crds1[1];
    for (i = 0; i < grid_z; i++)
        data_z_new[i] = (izmin + i) * hzgrid_new - crds1[2];

    for(ix = 0;  ix < grid_x; ix++)
    {
        indx[ix] = -1;

        for(i=0; i < grid_x-1; i++)
        {
            if(data_x_new[ix] >= data_x[i] && data_x_new[ix] < data_x[i+1]) 
            {
                indx[ix] = i;
                frac_x[ix] = (data_x_new[ix] - data_x[i])/hxgrid;
                break;
            }
        }

    }
            
    for(iy = 0;  iy < grid_y; iy++)
    {
        indy[iy] = -1;

        for(i=0; i < grid_y; i++)
        {
            if(data_y_new[iy] >= data_y[i] && data_y_new[iy] < data_y[i+1]) 
            {
                indy[iy] = i;
                frac_y[iy] = (data_y_new[iy] - data_y[i])/hygrid;
                break;
            }
        }
    }


    for(iz = 0;  iz < grid_z; iz++)
    {
        indz[iz] = -1;

        for(i=0; i < grid_z; i++)
        {
            if(data_z_new[iz] >= data_z[i] && data_z_new[iz] < data_z[i+1]) 
            {
                indz[iz] = i;
                frac_z[iz] = (data_z_new[iz] - data_z[i])/hzgrid;
                break;
            }
        }
    }



    for(ix = 0;  ix < grid_x; ix++)
    {
        for(iy = 0;  iy < grid_y; iy++)
        {
            for(iz = 0;  iz < grid_z; iz++)
            {
                idx = ix * grid_y * grid_z + iy * grid_z + iz;

                if(      indx[ix] <0 | indx[ix] >= grid_x 
                        |indy[iy] <0 | indy[iy] >= grid_y 
                        |indz[iz] <0 | indz[iz] >= grid_z)
                {
                   rho_out[idx] = 0.0;
                }

                else 
                {

                    idx000 = (indx[ix]+0) * incx + (indy[iy]+0) * incy + indz[iz]+0;
                    idx100 = (indx[ix]+1) * incx + (indy[iy]+0) * incy + indz[iz]+0;
                    idx010 = (indx[ix]+0) * incx + (indy[iy]+1) * incy + indz[iz]+0;
                    idx001 = (indx[ix]+0) * incx + (indy[iy]+0) * incy + indz[iz]+1;
                    idx110 = (indx[ix]+1) * incx + (indy[iy]+1) * incy + indz[iz]+0;
                    idx011 = (indx[ix]+0) * incx + (indy[iy]+1) * incy + indz[iz]+1;
                    idx101 = (indx[ix]+1) * incx + (indy[iy]+0) * incy + indz[iz]+1;
                    idx111 = (indx[ix]+1) * incx + (indy[iy]+1) * incy + indz[iz]+1;

                    frac000 = (1.0-frac_x[ix]) * (1.0 - frac_y[iy]) * (1.0 - frac_z[iz]);
                    frac100 = (    frac_x[ix]) * (1.0 - frac_y[iy]) * (1.0 - frac_z[iz]);
                    frac010 = (1.0-frac_x[ix]) * (      frac_y[iy]) * (1.0 - frac_z[iz]);
                    frac001 = (1.0-frac_x[ix]) * (1.0 - frac_y[iy]) * (      frac_z[iz]);
                    frac110 = (    frac_x[ix]) * (      frac_y[iy]) * (1.0 - frac_z[iz]);
                    frac011 = (1.0-frac_x[ix]) * (      frac_y[iy]) * (      frac_z[iz]);
                    frac101 = (    frac_x[ix]) * (1.0 - frac_y[iy]) * (      frac_z[iz]);
                    frac111 = (    frac_x[ix]) * (      frac_y[iy]) * (      frac_z[iz]);

                    rho_out[idx] =  frac000 * rho_tem[idx000]
                                  + frac100 * rho_tem[idx100]
                                  + frac010 * rho_tem[idx010]
                                  + frac001 * rho_tem[idx001]
                                  + frac110 * rho_tem[idx110]
                                  + frac011 * rho_tem[idx011]
                                  + frac101 * rho_tem[idx101]
                                  + frac111 * rho_tem[idx111];
                }
            }
        }
    }




    my_free(indx);
    my_free(indy);
    my_free(indz);
    my_free(frac_x);
    my_free(frac_y);
    my_free(frac_z);
    my_free(data_x_new);
    my_free(data_y_new);
    my_free(data_z_new);
    my_free(data_x);
    my_free(data_y);
    my_free(data_z);

    return;

}
