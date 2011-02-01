/************************** SVN Revision Information **************************
 **    $Id: init_rho_atom.c 648 2006-10-19 17:35:46Z luw $    **
******************************************************************************/
 
/*
     	Just generates a random start.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "md.h"


void init_rho_atom(double *rho)
{

    int idx1, idx, state, ix, iy, iz;
    int ix1, iy1, iz1;
    int ix0, iy0, iz0;
    int ixx, iyy, izz;
    char newname[MAX_PATH + 200];
    int i, ii;
    int ion, species, ist, fhand, nbytes;

    double hxgrid, hygrid, hzgrid;
    double hxgrid_new, hygrid_new, hzgrid_new;
    double *rho_tem, *rho_out, crds[3], *crds1;
    int ixdim, iydim, izdim;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int pex, pey, pez;
    int *map;

    my_malloc_init( map, ct.num_ions, int );
    hxgrid_new = ct.hxxgrid * ct.xside;
    hygrid_new = ct.hyygrid * ct.yside;
    hzgrid_new = ct.hzzgrid * ct.zside;

    pe2xyz(pct.thispe, &pex, &pey, &pez);

    ix0 = pex * FPX0_GRID;
    iy0 = pey * FPY0_GRID;
    iz0 = pez * FPZ0_GRID;
    ix1 = ix0 + FPX0_GRID;
    iy1 = iy0 + FPY0_GRID;
    iz1 = iz0 + FPZ0_GRID;



    if (pct.thispe == 0)
        printf(" initial rho from atom \n");
    my_barrier();

    ixdim = 0;
    iydim = 0;
    izdim = 0;

    for (idx = 0; idx < FP0_BASIS; idx++)
    {
        rho[idx] = 0.0;
    }


    for (ion = 0; ion < ct.num_ions; ion++)
    {
        species = ct.ions[ion].species;

        sprintf(newname, "%s%s", ct.file_atomic_orbit[species], ".rho_firstatom");
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
        {
            printf("\n unable to open file: %s \n", newname);
            error_handler(" Unable to open file ");
        }

        read(fhand, &ixmin, sizeof(int));
        read(fhand, &ixmax, sizeof(int));
        read(fhand, &iymin, sizeof(int));
        read(fhand, &iymax, sizeof(int));
        read(fhand, &izmin, sizeof(int));
        read(fhand, &izmax, sizeof(int));
        read(fhand, &hxgrid, sizeof(double));
        read(fhand, &hygrid, sizeof(double));
        read(fhand, &hzgrid, sizeof(double));
        read(fhand, &crds[0], 3 * sizeof(double));
        ixdim = max(ixdim, ixmax - ixmin);
        iydim = max(iydim, iymax - iymin);
        izdim = max(izdim, izmax - izmin);

        close(fhand);

        crds1 = &ct.ions[ion].crds[0];
        ixmin += (int) (crds1[0] / hxgrid_new) - (int) (crds[0] / hxgrid);
        iymin += (int) (crds1[1] / hygrid_new) - (int) (crds[1] / hygrid);
        izmin += (int) (crds1[2] / hzgrid_new) - (int) (crds[2] / hzgrid);

        ixmax += (int) (crds1[0] / hxgrid_new) - (int) (crds[0] / hxgrid);
        iymax += (int) (crds1[1] / hygrid_new) - (int) (crds[1] / hygrid);
        izmax += (int) (crds1[2] / hzgrid_new) - (int) (crds[2] / hzgrid);

        map[ion] = 0;
        for (ix = ixmin; ix < ixmax; ix++)
        {
            for (iy = iymin; iy < iymax; iy++)
            {
                for (iz = izmin; iz < izmax; iz++)
                {

                    ixx = ix;
                    iyy = iy;
                    izz = iz;
                    if (ixx < 0)
                        ixx += FNX_GRID;
                    if (iyy < 0)
                        iyy += FNY_GRID;
                    if (izz < 0)
                        izz += FNZ_GRID;
                    if (ixx >= FNX_GRID)
                        ixx -= FNX_GRID;
                    if (iyy >= FNY_GRID)
                        iyy -= FNY_GRID;
                    if (izz >= FNZ_GRID)
                        izz -= FNZ_GRID;

                    if (ixx >= ix0 && ixx < ix1
                        && iyy >= iy0 && iyy < iy1 && izz >= iz0 && izz < iz1)
                        map[ion] = 1;
                }
            }
        }


    }

    my_malloc_init( rho_tem, ixdim * iydim * izdim, REAL );
    my_malloc_init( rho_out, ixdim * iydim * izdim, REAL );



    for (ion = 0; ion < ct.num_ions; ion++)
    {

        if (map[ion] == 1)
        {
            species = ct.ions[ion].species;
            crds1 = &ct.ions[ion].crds[0];

            sprintf(newname, "%s%s", ct.file_atomic_orbit[species], ".rho_firstatom");
            fhand = open(newname, O_RDWR);
            if (fhand < 0)
            {
                printf("\n unable to open file: %s \n", newname);
                error_handler(" Unable to open file ");
            }

            read(fhand, &ixmin, sizeof(int));
            read(fhand, &ixmax, sizeof(int));
            read(fhand, &iymin, sizeof(int));
            read(fhand, &iymax, sizeof(int));
            read(fhand, &izmin, sizeof(int));
            read(fhand, &izmax, sizeof(int));
            read(fhand, &hxgrid, sizeof(double));
            read(fhand, &hygrid, sizeof(double));
            read(fhand, &hzgrid, sizeof(double));
            read(fhand, &crds[0], 3 * sizeof(double));

            ixdim = ixmax - ixmin;
            iydim = iymax - iymin;
            izdim = izmax - izmin;
            idx = ixdim * iydim * izdim * sizeof(double);
            read(fhand, rho_tem, idx);
            close(fhand);

            /* rho_tem: input old one, output  new interpolated one
             * ixmin, ... input for old, output for new
             */


            interpolate_atom_density(rho_tem, rho_out, ixmin, ixmax, iymin, iymax, izmin,
                                 izmax, crds, crds1, hxgrid, hygrid, hzgrid,
                                 hxgrid_new, hygrid_new, hzgrid_new);
            ixmin += (int) (crds1[0] / hxgrid_new) - (int) (crds[0] / hxgrid);
            iymin += (int) (crds1[1] / hygrid_new) - (int) (crds[1] / hygrid);
            izmin += (int) (crds1[2] / hzgrid_new) - (int) (crds[2] / hzgrid);
            ixmax = ixmin + ixdim;
            iymax = iymin + iydim;
            izmax = izmin + izdim;


            for (ix = ixmin; ix < ixmax; ix++)
            {
                for (iy = iymin; iy < iymax; iy++)
                {
                    for (iz = izmin; iz < izmax; iz++)
                    {

                        ixx = ix;
                        iyy = iy;
                        izz = iz;
                        if (ixx < 0)
                            ixx += FNX_GRID;
                        if (iyy < 0)
                            iyy += FNY_GRID;
                        if (izz < 0)
                            izz += FNZ_GRID;
                        if (ixx >= FNX_GRID)
                            ixx -= FNX_GRID;
                        if (iyy >= FNY_GRID)
                            iyy -= FNY_GRID;
                        if (izz >= FNZ_GRID)
                            izz -= FNZ_GRID;

                        if (ixx >= ix0 && ixx < ix1
                            && iyy >= iy0 && iyy < iy1 && izz >= iz0 && izz < iz1)

                        {
                            idx = (ix - ixmin) * iydim * izdim + (iy - iymin) * izdim + iz - izmin;
                            idx1 =
                                (ixx - ix0) * FPY0_GRID * FPZ0_GRID + (iyy -
                                                                       iy0) * FPZ0_GRID + izz - iz0;

                            rho[idx1] += rho_out[idx];


                        }

                    }
                }
            }
        }
    }


    my_barrier();
    my_free(rho_tem);
    my_free(map);
    if (pct.thispe == 0)
        printf(" initial rho  done  \n");


}                               /* end init_rho_atom */
