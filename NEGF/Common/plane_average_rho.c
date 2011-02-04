/************************** SVN Revision Information **************************
 **    $Id: plane_average_rho.c 1242 2011-02-02 18:55:23Z luw $ 
******************************************************************************/
 
/*

   write out plane average of Hatree potential and charge density
   ct.plane[] = (0, 0, 1) for xy plane
   ct.plane[] = (0, 1, 0) for xz plane
   ct.plane[] = (1, 0, 0) for yz plane

   no other plane are programed in this routine
   ct.plane[3] and ct.plane[4] define the arange for average. 


*/


#include <stdio.h>
#include "md.h"




void plane_average_rho (double *rho)
{




    int ix, iy, iz;
    int px, py, pz;
    REAL t1;
    REAL *zvec;
    int pxoff, pyoff, pzoff;
    FILE *file;


    /* Get this processors offset */
    pe2xyz (pct.thispe, &px, &py, &pz);
    pxoff = px * FPX0_GRID;
    pyoff = py * FPY0_GRID;
    pzoff = pz * FPZ0_GRID;

    if(ct.plane[0] + ct.plane[1] + ct.plane[2] >1)
    {
        printf("\n plane is not defined right %d %d %d \n", ct.plane[0], ct.plane[1],ct.plane[2]);
    }

    
    t1 = ct.plane[4] - ct.plane[3];

    if(ct.plane[0] == 1) 
    {

        my_calloc(zvec, FNY_GRID * FNZ_GRID, REAL );


        /* Loop over this processor */
        for (ix = 0; ix < FPX0_GRID; ix++)
        {
            if(ix + pxoff >= ct.plane[3] && ix + pxoff < ct.plane[4])
            {
                for (iy = 0; iy < FPY0_GRID; iy++)
                {
                    for (iz = 0; iz < FPZ0_GRID; iz++)
                    {
                        zvec[ (iy+pyoff) * FNZ_GRID + iz + pzoff]  += rho[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz];
                    }
                }
                
            }
        }



        /* Now sum over all processors */
        ix = FNY_GRID * FNZ_GRID;
        global_sums (zvec, &ix);

        if (pct.thispe == 0)
        {
            file = fopen("rho_yz_plane.dat", "w");
            fprintf (file, "\n average of rho from nx = %d to %d", ct.plane[3], ct.plane[4]);
            for (ix = 0; ix < FNY_GRID * FNZ_GRID; ix++)
            {
                fprintf(file, "\n %e", zvec[ix]/t1);
            }
            fclose(file);
        }

    }  



    if(ct.plane[1] == 1) 
    {

        my_calloc(zvec, FNX_GRID * FNZ_GRID, REAL );


        /* Loop over this processor */
        for (ix = 0; ix < FPX0_GRID; ix++)
        {
            for (iy = 0; iy < FPY0_GRID; iy++)
            {
                if(iy + pyoff >= ct.plane[3] && iy + pyoff < ct.plane[4])
                {
                    for (iz = 0; iz < FPZ0_GRID; iz++)
                    {
                        zvec[ (ix+pxoff) * FNZ_GRID + iz + pzoff]  += rho[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz];
                    }
                }

            }
        }



        /* Now sum over all processors */
        ix = FNX_GRID * FNZ_GRID;
        global_sums (zvec, &ix);

        if (pct.thispe == 0)
        {
            file = fopen("rho_xz_plane.dat", "w");
            fprintf (file, "\n average of rho from ny = %d to %d", ct.plane[3], ct.plane[4]);
            for (ix = 0; ix < FNX_GRID * FNZ_GRID; ix++)
            {
                fprintf(file, "\n %e", zvec[ix]/t1);
            }
            fclose(file);
        }

    }  




    if(ct.plane[2] == 1) 
    {

        my_calloc(zvec, FNX_GRID * FNY_GRID, REAL );


        /* Loop over this processor */
        for (ix = 0; ix < FPX0_GRID; ix++)
        {
            for (iy = 0; iy < FPY0_GRID; iy++)
            {
                for (iz = 0; iz < FPZ0_GRID; iz++)
                {
                    if(iz + pzoff >= ct.plane[3] && iz + pzoff < ct.plane[4])
                    {
                        zvec[ (ix+pxoff) * FNY_GRID + iy + pyoff]  += rho[ix * FPY0_GRID * FPZ0_GRID + iy * FPZ0_GRID + iz];
                    }
                }

            }
        }



        /* Now sum over all processors */
        ix = FNX_GRID * FNY_GRID;
        global_sums (zvec, &ix);

        if (pct.thispe == 0)
        {
            file = fopen("rho_xy_plane.dat", "w");
            fprintf (file, "\n average of rho from nz = %d to %d", ct.plane[3], ct.plane[4]);
            for (ix = 0; ix < FNX_GRID * FNY_GRID; ix++)
            {
                fprintf(file, "\n %e", zvec[ix]/t1);
            }
            fclose(file);
        }

    }  



    my_free( zvec );
}


