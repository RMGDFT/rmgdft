/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

/* Prints density to at points (grid_x, grid_y, z) - first two coordiantes being fixed - to standard output*/


#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"


/*Input parameters: grid_x - global x grid coordinate (0 <= grid_x < PE_X*px0_grid) 
 *                  grid_y - global y grid coordinate (0 <= grid_y < PE_Y*py0_grid) 
 *                  density; array defined on px0_grid*py0_grid*pz0_grid grid
 * 		    px0_grid, py0_grid, pz0_grid : Number of grid points in x, y and z directions (per processor)
 * 		    zside: length of cell in z direction */

void print_density_z_direction (int grid_x, int grid_y, double * density, int px0_grid, int py0_grid,
                                int pz0_grid, double zside)
{
    int i, j, ii, jj, kk, counter;
    MPI_Status mstatus;
    double *temp_buff;
    int p0_basis, xgrid, ygrid;
    int min_grid_x, min_grid_y, max_grid_x, max_grid_y;


    /*Check consistency */
    if ((grid_x >= get_PE_X() * px0_grid) || (grid_x < 0))
        error_handler ("grid_x is specified incorrectly");
    if ((grid_y >= get_PE_Y() * py0_grid) || (grid_y < 0))
        error_handler ("grid_y is specified incorrectly");

    p0_basis = px0_grid * py0_grid * pz0_grid;

    my_malloc (temp_buff, p0_basis, double);

    pe2xyz (pct.gridpe, &ii, &jj, &kk);

    /*Max and minimum global grid for given processor */
    min_grid_x = ii * px0_grid;
    min_grid_y = jj * py0_grid;

    max_grid_x = min_grid_x + px0_grid;
    max_grid_y = min_grid_y + py0_grid;


    /*printf("\nPE %d: ii %d jj %d kk %d",pct.gridpe, ii, jj, kk);
       fflush(NULL);

       MPI_Barrier(pct.grid_comm); */


    for (i = 0; i < get_PE_Z(); i++)
    {

        /*only one processor will be sending data at this time */
        if ((min_grid_x <= grid_x) && (grid_x < max_grid_x)
            && (min_grid_y <= grid_y) && (grid_y < max_grid_y) && (kk == i))
        {

            counter = 0;
            for (j = 0; j < p0_basis; j++)
            {

                /*Absolute grid coordiantes for grid point j */
                xgrid = j / (pz0_grid * py0_grid) + ii * px0_grid;
                ygrid = (j % (pz0_grid * py0_grid)) / pz0_grid + jj * py0_grid;

                if ((xgrid == grid_x) && (ygrid == grid_y))
                {
                    temp_buff[counter] = density[j];
                    counter++;
                }


            }                   /*end for (j=1;j<P0_BASIS;j++) */

            /*printf("\nPE %d: ii %d jj %d kk %d, sending %d data points to 0",pct.gridpe, ii, jj, kk, counter);
               printf("\nPE %d: Sent temp_buff[0]=%e",pct.gridpe, temp_buff[0]); */

            /*check */
            if (counter != pz0_grid)
            {
                printf ("\n\n PE %d: Counter is %d, expected %d !!", pct.gridpe, counter, pz0_grid);
                fflush (NULL);
            }

            /*send data to PE 0 */
            MPI_Send (temp_buff, pz0_grid, MPI_DOUBLE, 0, 100, pct.grid_comm);

        }                       /*end if ((ii == ((double) PE_X)/2.0) ... */

        /*PE 0 will write the data out */
        if (pct.gridpe == 0)
        {

            MPI_Recv (temp_buff, pz0_grid, MPI_DOUBLE, MPI_ANY_SOURCE, 100, pct.grid_comm,
                      &mstatus);

            for (j = 0; j < pz0_grid; j++)
                printf ("\n%f\t%e", (i * pz0_grid + j) / ((double) (pz0_grid * get_PE_Z())) * zside,
                        temp_buff[j]);

            fflush (NULL);

        }                       /*end if(pct.gridpe == 0) */

        MPI_Barrier (pct.grid_comm);

    }                           /*end for (i=0;i<PE_Z;i++) */

    my_free (temp_buff);


}

 /*EOF*/
