/************************** SVN Revision Information **************************
 **    $Id: print_wave 2013-06-02 15:38:05Z BTAN $    **
******************************************************************************/
/*

print out wave_plot stored in wave_global in coarse grids to file

*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "md.h"
void print_wave(REAL * wave_global, int wave_plot, STATE * states)
{
    int  ix, iy, iz;
    FILE *file, *file1;
    char filename[20];
    char filename1[20];


    if (pct.gridpe == 0) // for printing, only root processor will do it
    {
        double dx = ct.celldm[0] / NX_GRID;
        double dy = ct.celldm[0] * ct.celldm[1] / NY_GRID;
        double dz = ct.celldm[0] * ct.celldm[2] / NZ_GRID;

        int count = 0;

	sprintf(filename, "wave_%d.cube", wave_plot);
	sprintf(filename1, "rho_%d.cube", wave_plot);

        file = fopen (filename, "w"); //create gaussian file to plot state_plot stored in eig_rho_global 
        file1 = fopen (filename1, "w"); //create gaussian file to plot state_plot stored in eig_rho_global 

        fprintf( file, "Cubfile created from PWScf calculation\n" );
        fprintf( file, "Plotted wave with eigenvalue = %f \n", states[wave_plot].eig * Ha_eV);
        fprintf( file, "1     0.000000    0.000000    0.000000 \n" );//hack the cube file by pretending there is only one atom in the gaussian file
        fprintf (file, "%d    %12.9f      0.000000    0.000000 \n", NX_GRID, dx );//dx is the grid spacing in x in bohr
        fprintf (file, "%d    0.000000    %12.9f      0.000000 \n", NY_GRID, dy );
        fprintf (file, "%d    0.000000    0.000000    %12.9f   \n", NZ_GRID, dz );
        fprintf (file, "6     6.000000   10.000000   10.000000   10.000000 \n");//hack file by assigning just one carbon atom at some random position

        for (ix = 0; ix < NX_GRID; ix++)
        {
                for (iy = 0; iy < NY_GRID; iy++)
                {
                        for (iz = 0; iz < NZ_GRID; iz++)
                        {
                                count ++;
                                fprintf ( file , " %18.6e",
                                         wave_global[ix * NY_GRID * NZ_GRID + iy  *  NZ_GRID + iz] );
                                if (count % 6 == 0)
                                        fprintf (file, "\n");
                        }
                }
        }


        fclose (file);

//print out the charge density plot for comparison!

        fprintf( file1, "Cubfile created from PWScf calculation\n" );
        fprintf( file1, "Plotted rho with eigenvalue = %f \n", states[wave_plot].eig * Ha_eV);
        fprintf( file1, "1     0.000000    0.000000    0.000000 \n" );
        fprintf (file1, "%d    %12.9f      0.000000    0.000000 \n", NX_GRID, dx );
        fprintf (file1, "%d    0.000000    %12.9f      0.000000 \n", NY_GRID, dy );
        fprintf (file1, "%d    0.000000    0.000000    %12.9f   \n", NZ_GRID, dz );
        fprintf (file1, "6     6.000000   10.000000   10.000000   10.000000 \n");

        for (ix = 0; ix < NX_GRID; ix++)
        {
                for (iy = 0; iy < NY_GRID; iy++)
                {
                        for (iz = 0; iz < NZ_GRID; iz++)
                        {
                                count ++;
                                fprintf ( file1 , " %18.6e",
                                       wave_global[ix * NY_GRID * NZ_GRID + iy * NZ_GRID + iz] * wave_global[ix * NY_GRID * NZ_GRID + iy * NZ_GRID + iz] );
                                if (count % 6 == 0)
                                        fprintf (file1, "\n");
                        }
                }
        }


        fclose (file1);
    }

}
