/************************** SVN Revision Information **************************
 **    $Id: write_global_data.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*
                        write_data.c

    Functions to write data to files.
*/


#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "md.h"



/* Writes global data to a file */

void write_global_data (int file_handle, double *data, int NX, int NY, int NZ)
{
    int i, j, k;
    int amode;
    double *x_plane;
    int size;
    int pe_x, pe_y, pe_z;
    int x_off, y_off, z_off;
    int ix, iy, iz;
    int global_index, local_index;

    my_barrier ();
    size = NY * NZ;
    my_malloc( x_plane, size, double );


    pe2xyz (pct.thispe, &pe_x, &pe_y, &pe_z);
    x_off = pe_x * NX / pct.pe_x;
    y_off = pe_y * NY / pct.pe_y;
    z_off = pe_z * NZ / pct.pe_z;

    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY * NZ; j++)
            x_plane[j] = 0.0;

        if ((i >= x_off) && (i < x_off + NX / pct.pe_x))
        {
            for (iy = 0; iy < NY / pct.pe_y; iy++)
            {
                for (iz = 0; iz < NZ / pct.pe_z; iz++)
                {
                    global_index = (iy + y_off) * NZ + iz + z_off;
                    local_index = (i - x_off) * NY * NZ / pct.pe_y / pct.pe_z + iy * NZ / pct.pe_z + iz;
                    x_plane[global_index] = data[local_index];
                }
            }
        }

        /* get one slice global data in a x plane */
        global_sums (x_plane, &size);


        if (pct.thispe == 0)
            write (file_handle, x_plane, NY * NZ * sizeof (double));

    }

    my_free( x_plane );

    my_barrier ();
}



void write_global_data_lead (int file_handle, double *data, int NX, int NY, int NZ)
{
    int i, j, k;
    int amode;
    double *x_plane;
    int size;
    int pe_x, pe_y, pe_z;
    int x_off, y_off, z_off;
    int ix, iy, iz;
    int global_index, local_index;
    double sum;

    my_barrier ();
    size = NY * NZ;
    my_malloc( x_plane, size, double );


    pe2xyz (pct.thispe, &pe_x, &pe_y, &pe_z);
    x_off = pe_x * NX / pct.pe_x;
    y_off = pe_y * NY / pct.pe_y;
    z_off = pe_z * NZ / pct.pe_z;

    for (i = 0; i < NX / 3; i++)
    {
        for (j = 0; j < NY * NZ; j++)
            x_plane[j] = 0.0;

        if ((i >= x_off) && (i < x_off + NX / pct.pe_x))
        {
            for (iy = 0; iy < NY / pct.pe_y; iy++)
            {
                for (iz = 0; iz < NZ / pct.pe_z; iz++)
                {
                    global_index = (iy + y_off) * NZ + iz + z_off;
                    local_index = (i - x_off) * NY * NZ / pct.pe_y / pct.pe_z
                        + iy * NZ / pct.pe_z + iz;
                    x_plane[global_index] = data[local_index];
                }
            }
        }

        /* get one slice global data in a x plane */
        global_sums (x_plane, &size);


        if (pct.thispe == 0)
            write (file_handle, x_plane, NY * NZ * sizeof (double));

    }

    my_free( x_plane );

    my_barrier ();
}
