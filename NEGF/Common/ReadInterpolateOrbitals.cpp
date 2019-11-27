#include "negf_prototypes.h"
 
/*

read_orbital.c

This routine reads the orbitals of all probes and the central part. 
Reading data for the left probe, right probe and the central part is 
(from ON cal) straightforward. However, for the down and up probes
the operation orbitals_NEGF(x,y,z) = orbitals_ON(y,x,z) is needed.

*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "init_var.h"
#include "LCR.h"
#include "LocalObject.h"
#include "transition.h"
#include "blas.h"

static void read_one_orbital(double *psi, int st1, int *ixmin, int *ixmax, int *iymin, int *iymax, int ipart);
static void interpolation_orbit (int st, double *psi_old, double *psi_new, 
        int ixmin_old, int iymin_old, int ipart);

void ReadInterpolateOrbitals ()
{

    int ixmin, ixmax, iymin, iymax;
    double *psi, *psi_rotate, *psi_new;

 
    /* Wait until everybody gets here */
    MPI_Barrier(pct.img_comm);

    psi = new double[ct.max_orbit_size];
    psi_new = new double[ct.max_orbit_size];
    psi_rotate = new double[ct.max_orbit_size];

    for(int st = 0; st < LocalOrbital->num_thispe; st++)
    {
        int st_glob = LocalOrbital->index_proj_to_global[st];

        int st_new = 0;
        for (int subsystem = 0; subsystem < cei.num_subsystem; subsystem++)
        {
            int ipart = cei.subsystem_idx[subsystem];

            for (int st1 = lcr[ipart].state_begin; st1 < lcr[ipart].state_end; st1++)
            {
                if(st_new == st_glob)
                {

                    read_one_orbital(psi, st1, &ixmin, &ixmax, &iymin, &iymax, ipart);

                    // for up and down probe, we swap x and y directions. 3 and/or 4 probe
                    if(ipart > 2 && ipart <= cei.num_probe) 
                    {

                        std::swap(ixmin, iymin);
                        std::swap(ixmax, iymax);

                        int incx = states[st_glob].orbit_nz * states[st_glob].orbit_ny;
                        int incy = states[st_glob].orbit_nz * states[st_glob].orbit_nx;
                        for(int ix = 0; ix < states[st_glob].orbit_nx; ix++)
                        {
                            for(int iy = 0; iy < states[st_glob].orbit_ny; iy++)
                            {
                                for(int iz = 0; iz < states[st_glob].orbit_nz; iz++)
                                {
                                    int idx = iz + iy * states[st].orbit_nz + ix * incx;
                                    int idx2= iz + ix * states[st].orbit_nz + iy * incy; 
                                    psi_rotate[idx] = psi[idx2];
                                }
                            }
                        }

                        int ione = 1;
                        dcopy(&states[st_glob].size, psi_rotate, &ione, psi, &ione);

                    }   
                    interpolation_orbit (st_glob, psi, psi_new, ixmin, iymin, ipart);

                    LocalOrbital->AssignOrbital(st, psi_new);

                }
                st_new++;
            }   /* st1 loop ends */
        }   /* subsystem loop ends */

    }

    delete [] psi;
    delete [] psi_new;
    delete [] psi_rotate;

    MPI_Barrier(pct.img_comm);

}                               /* end read_data */
static void read_one_orbital(double *psi, int st1, int *ixmin, int *ixmax, int *iymin, int *iymax, int ipart)
{
    char newname[500];
    int fhand;

    sprintf (newname, "%s_spin%d%s%d", lcr[ipart].name, pct.spinpe, ".orbit_", st1);

    fhand = open(newname, O_RDWR, S_IREAD | S_IWRITE );
    if (fhand < 0)
    {
        printf ("\n %s, st1 = %d", newname, st1);
        rmg_error_handler (__FILE__, __LINE__, " Unable to open file ");
    }

    int idx = states[st1].size * (int) sizeof (double);

    int nbytes = read (fhand, psi, idx);
    if (nbytes != idx)
    {
        printf ("\n read %d is different from %d for state %d", nbytes, idx, st1);
        rmg_error_handler (__FILE__, __LINE__, "Unexpected end of file orbit");
    }

    nbytes = read (fhand, ixmin, sizeof (int));
    nbytes = read (fhand, ixmax, sizeof (int));
    nbytes = read (fhand, iymin, sizeof (int));
    nbytes = read (fhand, iymax, sizeof (int));
    if (nbytes != sizeof (int))
    {
        printf ("\n read %d is different from %d for state %d", (int) nbytes,
                (int) sizeof (int), st1);
        rmg_error_handler (__FILE__, __LINE__, "Unexpected end of file orbit");
    }

}

static void interpolation_orbit (int st, double *psi_old, double *psi_new, 
        int ixmin_old, int iymin_old, int ipart)
{


    int idx1, incx;
    int NX, NY, NZ;
    int ixmin_new, iymin_new;
    int max_orbit_nx_ny;

    double *psi_1d_old, *psi_1d_new; 
    double hx_old, hx_new, hy_old, hy_new;
    double x1_old, x1_new, y1_old, y1_new;

    /*  when the transport calculation use different gird space hx because of the different length
     *  in lead and conductor, interpolate the readin orbitals, potential and rho to the current grid
     *  numbers of grids for both are same, the hx is only slightly different
     */

    hx_new = get_hxgrid() * get_xside();
    hy_new = get_hygrid() * get_yside();

    max_orbit_nx_ny = std::max(ct.max_orbit_nx, ct.max_orbit_ny);
    psi_1d_old = new double[max_orbit_nx_ny];
    psi_1d_new = new double[max_orbit_nx_ny];

#if 	LDEBUG
    printf ("\n PE: %d  xyside %f  %f   ", pct.gridpe, lcr[ipart].xside, lcr[ipart].yside);
    printf ("\n PE: %d  xy_shift %f  %f %f  ", pct.gridpe, lcr[ipart].x_shift, lcr[ipart].y_shift);
    printf ("\n PE: %d  get_NX_GRID() %d  %d ", pct.gridpe, lcr[ipart].NX_GRID, lcr[ipart].NY_GRID);
    printf ("\n PE: %d  num_states %d   ", pct.gridpe, lcr[ipart].num_states);
#endif

    NX = states[st].orbit_nx;
    NY = states[st].orbit_ny;
    NZ = states[st].orbit_nz;

    incx = NY * NZ;

    ixmin_new = states[st].ixmin;
    iymin_new = states[st].iymin;

    /* printf ("\n %d %d %d st ixmn\n", st, ixmin, ixmin_old);*/
    x1_new = ixmin_new * hx_new;
    y1_new = iymin_new * hy_new;


    if(cei.num_probe == 2)
    {
        hx_old = lcr[ipart].xside / lcr[ipart].NX_GRID;
        x1_old = lcr[ipart].x_shift + ixmin_old  * hx_old;
        int item = (x1_new - x1_old)/hx_old;
        if(item < -2)
            printf("\n inpterpolation_orb alon x error: ixmin_new old, shift %d %d %f", 
                    states[st].ixmin, ixmin_old, lcr[ipart].x_shift);

        /* ==== Interpolation along x-direction ==== */
        for (int iy = 0; iy < NY; iy++)
        {
            for (int iz = 0; iz < NZ; iz++)
            {
                for (int ix = 0; ix < NX; ix++)
                {
                    idx1 = iz + iy * NZ + ix * incx;
                    psi_1d_old[ix] = psi_old[idx1];
                }
                diff_hx_interpolation (st, psi_1d_new, psi_1d_old, NX, hx_new, hx_old, x1_new, x1_old);

                for (int ix = 0; ix < NX; ix++)
                {
                    idx1 = iz + iy * NZ + ix * incx;
                    psi_new[idx1] = psi_1d_new[ix];
                }
            }
        }

    }
    else if(cei.num_probe < 5)
    {
        hx_old = lcr[ipart].xside / lcr[ipart].NX_GRID;
        x1_old = lcr[ipart].x_shift + ixmin_old * hx_old;

        /* ==== Interpolation along x-direction ==== */
        for (int iy = 0; iy < NY; iy++)
        {
            for (int iz = 0; iz < NZ; iz++)
            {
                for (int ix = 0; ix < NX; ix++)
                {
                    idx1 = iz + iy * NZ + ix * incx;
                    psi_1d_old[ix] = psi_old[idx1];
                }
                diff_hx_interpolation (st, psi_1d_new, psi_1d_old, NX, hx_new, hx_old, x1_new, x1_old);

                for (int ix = 0; ix < NX; ix++)
                {
                    idx1 = iz + iy * NZ + ix * incx;
                    psi_new[idx1] = psi_1d_new[ix];
                }
            }
        }


        hy_old = lcr[ipart].yside / lcr[ipart].NY_GRID;
        y1_old = lcr[ipart].y_shift + iymin_old  * hy_old;
        /*
           y1_old = lcr[idx0].y_shift + iymin_old * hy_old;
           y1_new = (iymin_new + lcr[idx0].y0) * hy_new;
         */
        /*
           if(idx0 == 0 || idx0 == 3 || idx0 == 4)
           printf (" old_urgent %d %d %f %f %d %f %d %f\n", st, idx0, 
           hy_old, hy_new, iymin_old, y1_old, iymin_new, y1_new);
         */

        /* ==== Interpolation along y-direction ==== */

        for (int iz = 0; iz < NZ; iz++)
        {
            for (int ix = 0; ix < NX; ix++)
            {
                for (int iy = 0; iy < NY; iy++)
                {
                    idx1 = iz + iy * NZ + ix * incx;
                    psi_1d_old[iy] = psi_old[idx1];
                }
                /*
                   if(pct.gridpe ==0) 
                   printf (" urgent %d %d %f %f %f %f \n", idx0, NY, hy_new, hy_old, y1_new, y1_old);
                 */
                diff_hx_interpolation (st, psi_1d_new, psi_1d_old, NY, hy_new, hy_old, y1_new, y1_old); 

                for (int iy = 0; iy < NY; iy++)
                {
                    idx1 = iz + iy * NZ + ix * incx;
                    psi_new[idx1] = psi_1d_new[iy];
                    /*
                       if(pct.gridpe ==0) 
                       printf (" urgent %f %f \n", psi_old[iy], psi_new[iy]);
                     */
                }
            }
        }

    } 
    else
    {
        rmg_error_handler (__FILE__, __LINE__, " WARNING: Code may need to modify for num_probe > 4 ");
    }


    delete []psi_1d_old;
    delete []psi_1d_new;

#if 	LDEBUG
    printf ("\n %d %f %f state x x  %d ixmin_old %d  hx_old %f", st, x1, x1_old, ixmin,
            ixmin_old, hx_old);
    fflush (NULL);
#endif

}
