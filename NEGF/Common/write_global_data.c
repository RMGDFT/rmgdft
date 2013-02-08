/************************** SVN Revision Information **************************
 **    $Id$    **
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


    pe2xyz (pct.gridpe, &pe_x, &pe_y, &pe_z);
    x_off = pct.PX_OFFSET;
    y_off = pct.PY_OFFSET;
    z_off = pct.PZ_OFFSET;

    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY * NZ; j++)
            x_plane[j] = 0.0;

        if ((i >= x_off) && (i < x_off + pct.PX0_GRID))
        {
            for (iy = 0; iy < pct.PY0_GRID; iy++)
            {
                for (iz = 0; iz < pct.PZ0_GRID; iz++)
                {
                    global_index = (iy + y_off) * NZ + iz + z_off;
                    local_index = (i - x_off) * pct.PY0_GRID * pct.PZ0_GRID + iy * pct.PZ0_GRID + iz;
                    x_plane[global_index] = data[local_index];
                }
            }
        }

        /* get one slice global data in a x plane */
        global_sums (x_plane, &size, pct.grid_comm);


        if (pct.gridpe == 0)
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


    x_off = pct.PX_OFFSET;
    y_off = pct.PY_OFFSET;
    z_off = pct.PZ_OFFSET;

    for (i = 0; i < NX / 3; i++)
    {
        for (j = 0; j < NY * NZ; j++)
            x_plane[j] = 0.0;

        if ((i >= x_off) && (i < x_off + pct.PX0_GRID))
        {
            for (iy = 0; iy < pct.PY0_GRID; iy++)
            {
                for (iz = 0; iz < pct.PZ0_GRID; iz++)
                {
                    global_index = (iy + y_off) * NZ + iz + z_off;
                    local_index = (i - x_off) * pct.PY0_GRID * pct.PZ0_GRID
                        + iy * pct.PZ0_GRID + iz;
                    x_plane[global_index] = data[local_index];
                }
            }
        }

        /* get one slice global data in a x plane */
        global_sums (x_plane, &size, pct.grid_comm);


        if (pct.gridpe == 0)
            write (file_handle, x_plane, NY * NZ * sizeof (double));

    }

    my_free( x_plane );

    my_barrier ();
}
