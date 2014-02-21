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
#include "main_on.h"
void print_wave(int wave_plot, STATE * states, int coarse_level)
{
    int  ix, iy, iz;
    FILE *file, *file1;
    char filename[20];
    char filename1[20];
    REAL wave_value;


    if (pct.gridpe == 0) // for printing, only root processor will do it
    {
        double dx = ct.celldm[0] / get_NX_GRID();
        double dy = ct.celldm[0] * ct.celldm[1] / get_NY_GRID();
        double dz = ct.celldm[0] * ct.celldm[2] / get_NZ_GRID();

        int count = 0;

	sprintf(filename, "wave_%d.cube", wave_plot);

        file = fopen (filename, "w"); //create gaussian file to plot state_plot stored in eig_rho_global 

        fprintf (file, "Cubfile created from PWScf calculation\n" );
        fprintf (file, "Plotted wave with eigenvalue = %f \n", states[wave_plot].eig[0] * Ha_eV);
        fprintf (file, "1     0.000000    0.000000    0.000000 \n" );//hack the cube file by pretending there is only one atom in the gaussian file
        fprintf (file, "%d    %12.9f      0.000000    0.000000 \n", get_NX_GRID()/coarse_level, dx*coarse_level);//dx is the grid spacing in x in bohr
        fprintf (file, "%d    0.000000    %12.9f      0.000000 \n", get_NY_GRID()/coarse_level, dy*coarse_level);
        fprintf (file, "%d    0.000000    0.000000    %12.9f   \n", get_NZ_GRID()/coarse_level, dz*coarse_level);
        fprintf (file, "6     6.000000   10.000000   10.000000   10.000000 \n");//hack file by assigning just one carbon atom at some random position

        for (ix = 0; ix < get_NX_GRID()/coarse_level; ix++)
        {
                for (iy = 0; iy < get_NY_GRID()/coarse_level; iy++)
                {
                        for (iz = 0; iz < get_NZ_GRID()/coarse_level; iz++)
                        {
                                count ++;
                                wave_value = wave_global[ix*coarse_level * get_NY_GRID() * get_NZ_GRID() + iy*coarse_level * get_NZ_GRID() + iz*coarse_level] ;
                                fprintf ( file , " %18.6e", wave_value );
                                if (count % 6 == 0)
                                        fprintf (file, "\n");
                        }
                }
        }


        fclose (file);

//print out the charge density plot for comparison!
	sprintf(filename1, "rho_%d.cube", wave_plot);
        file1 = fopen (filename1, "w"); //create gaussian file to plot state_plot stored in eig_rho_global 

        fprintf (file1, "Cubfile created from PWScf calculation\n" );
        fprintf (file1, "Plotted rho with eigenvalue = %f \n", states[wave_plot].eig[0] * Ha_eV);
        fprintf (file1, "1     0.000000    0.000000    0.000000 \n" );
        fprintf (file1, "%d    %12.9f      0.000000    0.000000 \n", get_NX_GRID()/coarse_level, dx*coarse_level );
        fprintf (file1, "%d    0.000000    %12.9f      0.000000 \n", get_NY_GRID()/coarse_level, dy*coarse_level );
        fprintf (file1, "%d    0.000000    0.000000    %12.9f   \n", get_NZ_GRID()/coarse_level, dz*coarse_level );
        fprintf (file1, "6     6.000000   10.000000   10.000000   10.000000 \n");

        for (ix = 0; ix < get_NX_GRID()/coarse_level; ix++)
        {
                for (iy = 0; iy < get_NY_GRID()/coarse_level; iy++)
                {
                        for (iz = 0; iz < get_NZ_GRID()/coarse_level; iz++)
                        {
                                count ++;
                                wave_value = wave_global[ix*coarse_level * get_NY_GRID() * get_NZ_GRID() + iy*coarse_level * get_NZ_GRID() + iz*coarse_level] ;
                                fprintf ( file1 , " %18.6e", wave_value * wave_value );
                                if (count % 6 == 0)
                                        fprintf (file1, "\n");
                        }
                }
        }


        fclose (file1);
    }

}
