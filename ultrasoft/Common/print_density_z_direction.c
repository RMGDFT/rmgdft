/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/* Prints density to at points (grid_x, grid_y, z) - first two coordiantes being fixed - to standard output*/


#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#if MPI


/*Input parameters: grid_x - global x grid coordinate (0 <= grid_x < PE_X*px0_grid) 
 *                  grid_y - global y grid coordinate (0 <= grid_y < PE_Y*py0_grid) 
 *                  density; array defined on px0_grid*py0_grid*pz0_grid grid
 * 		    px0_grid, py0_grid, pz0_grid : Number of grid points in x, y and z directions (per processor)
 * 		    zside: length of cell in z direction */

void print_density_z_direction (int grid_x, int grid_y, REAL * density, int px0_grid, int py0_grid,
                                int pz0_grid, REAL zside)
{
    int i, j, ii, jj, kk, counter;
    MPI_Status mstatus;
    REAL *temp_buff;
    int p0_basis, xgrid, ygrid, zgrid;
    int min_grid_x, min_grid_y, min_grid_z, max_grid_x, max_grid_y, max_grid_z;


    /*Check consistency */
    if ((grid_x >= PE_X * px0_grid) || (grid_x < 0))
        error_handler ("grid_x is specified incorrectly");
    if ((grid_y >= PE_Y * py0_grid) || (grid_y < 0))
        error_handler ("grid_y is specified incorrectly");

    p0_basis = px0_grid * py0_grid * pz0_grid;

    my_malloc (temp_buff, p0_basis, REAL);

    pe2xyz (pct.thispe, &ii, &jj, &kk);

    /*Max and minimum global grid for given processor */
    min_grid_x = ii * px0_grid;
    min_grid_y = jj * py0_grid;
    min_grid_z = kk * pz0_grid;

    max_grid_x = min_grid_x + px0_grid;
    max_grid_y = min_grid_y + py0_grid;
    max_grid_z = min_grid_z + pz0_grid;


    /*printf("\nPE %d: ii %d jj %d kk %d",pct.thispe, ii, jj, kk);
       fflush(NULL);

       MPI_Barrier(pct.grid_comm); */


    for (i = 0; i < PE_Z; i++)
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
                zgrid = j % pz0_grid + kk * pz0_grid;

                if ((xgrid == grid_x) && (ygrid == grid_y))
                {
                    temp_buff[counter] = density[j];
                    counter++;
                }


            }                   /*end for (j=1;j<P0_BASIS;j++) */

            /*printf("\nPE %d: ii %d jj %d kk %d, sending %d data points to 0",pct.thispe, ii, jj, kk, counter);
               printf("\nPE %d: Sent temp_buff[0]=%e",pct.thispe, temp_buff[0]); */

            /*check */
            if (counter != pz0_grid)
            {
                printf ("\n\n PE %d: Counter is %d, expected %d !!", pct.thispe, counter, pz0_grid);
                fflush (NULL);
            }

            /*send data to PE 0 */
            MPI_Send (temp_buff, pz0_grid, MPI_DOUBLE, 0, 100, pct.grid_comm);

        }                       /*end if ((ii == ((double) PE_X)/2.0) ... */

        /*PE 0 will write the data out */
        if (pct.thispe == 0)
        {

            MPI_Recv (temp_buff, pz0_grid, MPI_DOUBLE, MPI_ANY_SOURCE, 100, pct.grid_comm,
                      &mstatus);

            for (j = 0; j < pz0_grid; j++)
                printf ("\n%f\t%e", (i * pz0_grid + j) / ((REAL) (pz0_grid * PE_Z)) * zside,
                        temp_buff[j]);

            fflush (NULL);

        }                       /*end if(pct.thispe == 0) */

        MPI_Barrier (pct.grid_comm);

    }                           /*end for (i=0;i<PE_Z;i++) */

    my_free (temp_buff);


}

 /*EOF*/
#else
void print_density_z_direction (REAL * density, int px0_grid, int py0_grid, int pz0_grid)
{
    error_handler ("This function is implemented for parallel mpi machines only");
}
#endif
