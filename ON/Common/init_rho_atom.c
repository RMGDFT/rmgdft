/************************** SVN Revision Information **************************
 **    $Id$    **
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
#include "main.h"
#include "prototypes_on.h"


void init_rho_atom(double *rho)
{

    int idx1, idx, ix, iy, iz;
    int ix1, iy1, iz1;
    int ix0, iy0, iz0;
    int ixx, iyy, izz;
    char newname[MAX_PATH + 200];
    int ion, species, fhand;

    double hxgrid, hygrid, hzgrid;
    double hxgrid_new, hygrid_new, hzgrid_new;
    double *rho_tem, *rho_out, crds[3], *crds1;
    int ixdim, iydim, izdim;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int pex, pey, pez;
    int *map, mapx, mapy, mapz;
    int ibuf[6];
    double dbuf[6];

    void *RT = BeginRmgTimer("1-TOTAL: init: init_rho");
    my_malloc_init( map, ct.num_ions, int );
    hxgrid_new = get_hxxgrid() * get_xside();
    hygrid_new = get_hyygrid() * get_yside();
    hzgrid_new = get_hzzgrid() * get_zside();

    pe2xyz(pct.gridpe, &pex, &pey, &pez);

    ix0 = pex * get_FPX0_GRID();
    iy0 = pey * get_FPY0_GRID();
    iz0 = pez * get_FPZ0_GRID();
    ix1 = ix0 + get_FPX0_GRID();
    iy1 = iy0 + get_FPY0_GRID();
    iz1 = iz0 + get_FPZ0_GRID();



    if (pct.gridpe == 0)
        printf(" initial rho from atom \n");
    my_barrier();

    ixdim = 0;
    iydim = 0;
    izdim = 0;

    for (idx = 0; idx < get_FP0_BASIS(); idx++)
    {
        rho[idx] = 0.0;
    }


    void *RT1 = BeginRmgTimer("1-TOTAL: init: init_rho_map");
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        species = ct.ions[ion].species;

        sprintf(newname, "%s%s%s", pct.image_path[pct.thisimg], ct.file_atomic_orbit[species], ".rho_firstatom");
        if(pct.gridpe == 0)
        {
            fhand = open(newname, O_RDWR);
            if (fhand < 0)
            {
                printf("\n unable to open file: %s \n", newname);
                error_handler(" Unable to open file ");
            }

            read(fhand, ibuf, 6*sizeof(int));
            read(fhand, dbuf, 6*sizeof(double));
            close(fhand);
        }

        
        MPI_Bcast(ibuf, 6, MPI_INT, 0, pct.grid_comm);
        MPI_Bcast(dbuf, 6, MPI_DOUBLE, 0, pct.grid_comm);

        ixmin = ibuf[0];
        ixmax = ibuf[1];
        iymin = ibuf[2];
        iymax = ibuf[3];
        izmin = ibuf[4];
        izmax = ibuf[5];
        hxgrid = dbuf[0];
        hygrid = dbuf[1];
        hzgrid = dbuf[2];
        crds[0] = dbuf[3];
        crds[1] = dbuf[4];
        crds[2] = dbuf[5];
        ixdim = rmg_max(ixdim, ixmax - ixmin);
        iydim = rmg_max(iydim, iymax - iymin);
        izdim = rmg_max(izdim, izmax - izmin);


        crds1 = &ct.ions[ion].crds[0];
        ixmin += (int) (crds1[0] / hxgrid_new) - (int) (crds[0] / hxgrid);
        iymin += (int) (crds1[1] / hygrid_new) - (int) (crds[1] / hygrid);
        izmin += (int) (crds1[2] / hzgrid_new) - (int) (crds[2] / hzgrid);

        ixmax += (int) (crds1[0] / hxgrid_new) - (int) (crds[0] / hxgrid);
        iymax += (int) (crds1[1] / hygrid_new) - (int) (crds[1] / hygrid);
        izmax += (int) (crds1[2] / hzgrid_new) - (int) (crds[2] / hzgrid);

        map[ion] = 0;
        mapx = 0;
        for (ix = ixmin; ix < ixmax; ix++)
        {
            ixx = (ix + get_FNX_GRID()) % get_FNX_GRID();
            if (ixx >= ix0 && ixx < ix1)
                mapx = 1;
        }
        if(!mapx) continue;

        mapy = 0;
        for (iy = iymin; iy < iymax; iy++)
        {
            iyy = (iy + get_FNY_GRID()) % get_FNY_GRID();
            if(iyy >= iy0 && iyy < iy1)
                mapy = 1;
        }
        if(!mapy) continue;
        mapz = 0;
        for (iz = izmin; iz < izmax; iz++)
        {
            izz = (iz + get_FNZ_GRID()) % get_FNZ_GRID();
            if ( izz >= iz0 && izz < iz1)
                mapz = 1;
        }
        if(!mapz) continue;

        map[ion] = 1;
    }



    EndRmgTimer(RT1);
    my_malloc_init( rho_tem, ixdim * iydim * izdim, double );
    my_malloc_init( rho_out, ixdim * iydim * izdim, double );



    void *RT2 = BeginRmgTimer("1-TOTAL: init: init_rho_value");
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        if (map[ion] == 1)
        {
            species = ct.ions[ion].species;
            crds1 = &ct.ions[ion].crds[0];
            ION *iptr;
            iptr = &ct.ions[ion];
            
            sprintf(newname, "%s%s%s", pct.image_path[pct.thisimg], ct.file_atomic_orbit[species], ".rho_firstatom");
            fhand = open(newname, O_RDWR);
            if (fhand < 0)
            {
                printf("\n unable to open file: %s \n", newname);
                error_handler(" Unable to open file ");
            }

            read(fhand, ibuf, 6*sizeof(int));
            ixmin = ibuf[0];
            ixmax = ibuf[1];
            iymin = ibuf[2];
            iymax = ibuf[3];
            izmin = ibuf[4];
            izmax = ibuf[5];
            read(fhand, dbuf, 6*sizeof(double));
            hxgrid = dbuf[0];
            hygrid = dbuf[1];
            hzgrid = dbuf[2];
            crds[0] = dbuf[3];
            crds[1] = dbuf[4];
            crds[2] = dbuf[5];

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
                            ixx += get_FNX_GRID();
                        if (iyy < 0)
                            iyy += get_FNY_GRID();
                        if (izz < 0)
                            izz += get_FNZ_GRID();
                        if (ixx >= get_FNX_GRID())
                            ixx -= get_FNX_GRID();
                        if (iyy >= get_FNY_GRID())
                            iyy -= get_FNY_GRID();
                        if (izz >= get_FNZ_GRID())
                            izz -= get_FNZ_GRID();

                        if (ixx >= ix0 && ixx < ix1
                                && iyy >= iy0 && iyy < iy1 && izz >= iz0 && izz < iz1)

                        {
                            idx = (ix - ixmin) * iydim * izdim + (iy - iymin) * izdim + iz - izmin;
                            idx1 =
                                (ixx - ix0) * get_FPY0_GRID() * get_FPZ0_GRID() + (iyy -
                                        iy0) * get_FPZ0_GRID() + izz - iz0;

                            if (ct.spin_flag == 1) 
                            { 
                                if (pct.spinpe == 0)
                                    rho[idx1] += rho_out[idx] * (0.5 + iptr->init_spin_rho) ;
                                else
                                    rho[idx1] += rho_out[idx] * (0.5 - iptr->init_spin_rho) ;
                            }
                            else
                                rho[idx1] += rho_out[idx];


                        }

                    }
                }
            }
        }
    }


    my_barrier();
    EndRmgTimer(RT2);
    my_free(rho_tem);
    my_free(rho_out);
    my_free(map);
    if (pct.gridpe == 0)
        printf(" initial rho  done  \n");

    EndRmgTimer(RT);

}                               /* end init_rho_atom */


