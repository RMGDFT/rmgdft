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
#include "main.h"
#include "init_var_negf.h"
#include "LCR.h"



/* Writes global data to a file */

void write_global_data (int file_handle, double *data, int fnx, int fny, int fnz)
{
    int i, j, k;
    int amode;
    double *x_plane;
    int size;
    int x_off, y_off, z_off;
    int ix, iy, iz;
    int global_index, local_index;

    my_barrier ();
    size = fny * fnz;
    my_malloc( x_plane, size, double );


    x_off = get_FPX_OFFSET();
    y_off = get_FPY_OFFSET();
    z_off = get_FPZ_OFFSET();

    for (i = 0; i < fnx; i++)
    {
        for (j = 0; j < fny * fnz; j++)
            x_plane[j] = 0.0;

        if ((i >= x_off) && (i < x_off + get_FPX0_GRID()))
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                {
                    global_index = (iy + y_off) * fnz + iz + z_off;
                    local_index = (i - x_off) * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz;
                    x_plane[global_index] = data[local_index];
                }
            }
        }

        /* get one slice global data in a x plane */
        global_sums (x_plane, &size, pct.grid_comm);


        if (pct.gridpe == 0)
            write (file_handle, x_plane, fny * fnz * sizeof (double));

    }

    my_free( x_plane );

    my_barrier ();
}



void write_global_data_lead (int file_handle, double *data, int fnx, int fny, int fnz)
{
    int i, j, k;
    int amode;
    double *x_plane;
    int size;
    int x_off, y_off, z_off;
    int ix, iy, iz;
    int global_index, local_index;
    double sum;

    my_barrier ();
    size = fny * fnz;
    my_malloc( x_plane, size, double );


    x_off = get_FPX_OFFSET();
    y_off = get_FPY_OFFSET();
    z_off = get_FPZ_OFFSET();

    for (i = 0; i < fnx / 3; i++)
    {
        for (j = 0; j < fny * fnz; j++)
            x_plane[j] = 0.0;

        if ((i >= x_off) && (i < x_off + get_FPX0_GRID()))
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                {
                    global_index = (iy + y_off) * fnz + iz + z_off;
                    local_index = (i - x_off) * get_FPY0_GRID() * get_FPZ0_GRID()
                        + iy * get_FPZ0_GRID() + iz;
                    x_plane[global_index] = data[local_index];
                }
            }
        }

        /* get one slice global data in a x plane */
        global_sums (x_plane, &size, pct.grid_comm);


        if (pct.gridpe == 0)
            write (file_handle, x_plane, fny * fnz * sizeof (double));

    }

    my_free( x_plane );

    my_barrier ();
}
