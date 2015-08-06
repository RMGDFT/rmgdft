/************************** SVN Revision Information **************************
 **    $Id$ 
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
#include "main.h"
#include "init_var.h"
#include "LCR.h"




void plane_average_rho (double *rho)
{




    int ix, iy, iz;
    double t1;
    double *zvec;
    int pxoff, pyoff, pzoff;
    FILE *file;


    /* Get this processors offset */
    pxoff = get_FPX_OFFSET();
    pyoff = get_FPY_OFFSET();
    pzoff = get_FPZ_OFFSET();

    if(ct.plane[0] + ct.plane[1] + ct.plane[2] >1)
    {
        printf("\n plane is not defined right %d %d %d \n", ct.plane[0], ct.plane[1],ct.plane[2]);
    }

    
    t1 = ct.plane[4] - ct.plane[3];

    if(ct.plane[0] == 1) 
    {

        my_calloc(zvec, get_FNY_GRID() * get_FNZ_GRID(), double );


        /* Loop over this processor */
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
        {
            if(ix + pxoff >= ct.plane[3] && ix + pxoff < ct.plane[4])
            {
                for (iy = 0; iy < get_FPY0_GRID(); iy++)
                {
                    for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                    {
                        zvec[ (iy+pyoff) * get_FNZ_GRID() + iz + pzoff]  += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];
                    }
                }
                
            }
        }



        /* Now sum over all processors */
        ix = get_FNY_GRID() * get_FNZ_GRID();
        global_sums (zvec, &ix, pct.grid_comm);

        if (pct.gridpe == 0)
        {
            file = fopen("rho_yz_plane.dat", "w");
            fprintf (file, "\n average of rho from nx = %d to %d", ct.plane[3], ct.plane[4]);
            for (ix = 0; ix < get_FNY_GRID() * get_FNZ_GRID(); ix++)
            {
                fprintf(file, "\n %e", zvec[ix]/t1);
            }
            fclose(file);
        }

    }  



    if(ct.plane[1] == 1) 
    {

        my_calloc(zvec, get_FNX_GRID() * get_FNZ_GRID(), double );


        /* Loop over this processor */
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                if(iy + pyoff >= ct.plane[3] && iy + pyoff < ct.plane[4])
                {
                    for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                    {
                        zvec[ (ix+pxoff) * get_FNZ_GRID() + iz + pzoff]  += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];
                    }
                }

            }
        }



        /* Now sum over all processors */
        ix = get_FNX_GRID() * get_FNZ_GRID();
        global_sums (zvec, &ix, pct.grid_comm);

        if (pct.gridpe == 0)
        {
            file = fopen("rho_xz_plane.dat", "w");
            fprintf (file, "\n average of rho from ny = %d to %d", ct.plane[3], ct.plane[4]);
            for (ix = 0; ix < get_FNX_GRID() * get_FNZ_GRID(); ix++)
            {
                fprintf(file, "\n %e", zvec[ix]/t1);
            }
            fclose(file);
        }

    }  




    if(ct.plane[2] == 1) 
    {

        my_calloc(zvec, get_FNX_GRID() * get_FNY_GRID(), double );


        /* Loop over this processor */
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)
            {
                for (iz = 0; iz < get_FPZ0_GRID(); iz++)
                {
                    if(iz + pzoff >= ct.plane[3] && iz + pzoff < ct.plane[4])
                    {
                        zvec[ (ix+pxoff) * get_FNY_GRID() + iy + pyoff]  += rho[ix * get_FPY0_GRID() * get_FPZ0_GRID() + iy * get_FPZ0_GRID() + iz];
                    }
                }

            }
        }



        /* Now sum over all processors */
        ix = get_FNX_GRID() * get_FNY_GRID();
        global_sums (zvec, &ix, pct.grid_comm);

        if (pct.gridpe == 0)
        {
            file = fopen("rho_xy_plane.dat", "w");
            fprintf (file, "\n average of rho from nz = %d to %d", ct.plane[3], ct.plane[4]);
            for (ix = 0; ix < get_FNX_GRID() * get_FNY_GRID(); ix++)
            {
                fprintf(file, "\n %e", zvec[ix]/t1);
            }
            fclose(file);
        }

    }  



    my_free( zvec );
}


