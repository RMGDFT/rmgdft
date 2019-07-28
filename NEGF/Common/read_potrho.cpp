#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
 *
 * read_potrho.c
 *
 * This routine works with the following steps:
 * 1. Patches the potentials (and rhos) of all subsystems and store them 
 * in global array(s).
 * 2. While patching the vh is adjusted to the Fermi energy of the subsystem.
 * 3. For the 3rd, 4th probes use: pot/rho_NEGF(x,y,z) = pot/rho_ON(y,x,z).
 * 4. Interpolate data in a new fine-grid used in NEGF calculation. 
 * 5. Move the data from the global-array to the distributed array.
 *
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
#include "init_var.h"
#include "LCR.h"




void read_potrho (double *vh, int iflag, char *file_ex)
{
    long fhand;
    long nbytes, nbytes_first, nbytes_last;
    char newname[MAX_PATH + 200];
    char msg[200];
    MPI_Offset position;

    long idx0, idx, idx1, idx_sub;
    double x0_old, hx_new, hx_old;
    double y0_old, hy_new, hy_old;
    int NX0, NY0, NZ0;

    int subsystem, ix, iy, iz;
    int x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
    int FNYZ, NYZ0, NXZ0, FPYZ0, FNYPZ;
    double *array_tmp, tem;
    double *array_tmp_tmp;
    double *vh_global, *vh_all;
    double *xold_global, *xold, xseed, sfit;
    double *yold_global, *yold, yseed;
    double *vh_old, *vh_new;
    int dis1, dis2, dis12;
    double fac1, fac2, V1, V2;
    int i, ii, jj, kk;
    int MaxNumBytes = 2000000000;  //maximum number of bytes one can read from ON calculation for each read
    long NumBytesLeft;  //number of bytes one need to read from ON calculation for each read

    
    int amode,  subsizes[1], starts[1];
    int *rcount, *rdisp;
    MPI_Info fileinfo;
    MPI_Datatype  filetype, newtype;
    MPI_Status status;
    MPI_Offset disp;


    /* Wait until everybody gets here */
    MPI_Barrier(pct.img_comm);

    rcount = (int *) malloc(pct.grid_npes * sizeof(int));
    rdisp = (int *) malloc(pct.grid_npes * sizeof(int));

    FPYZ0 = get_FPZ0_GRID() * get_FPY0_GRID();
    FNYPZ = get_FPZ0_GRID() * get_FNY_GRID(); 
    FNYZ = get_FNY_GRID() * get_FNZ_GRID(); 

    ii = get_FPX_OFFSET();
    jj = get_FPY_OFFSET();
    kk = get_FPZ_OFFSET();


    hx_new = get_xside() * get_hxxgrid();
    hy_new = get_yside() * get_hyygrid();
    if(pct.gridpe ==0) printf (" hx_new, hy_new  =  %f  %f \n", hx_new, hy_new); 


/* =========== Allocate memory & Reading ON2 data ================ */


    idx = get_FNX_GRID() * get_FNY_GRID();
    my_malloc_init(xold_global, idx, double);   
    my_malloc_init(yold_global, idx, double);   
    idx = get_FNX_GRID() * get_FNY_GRID() * get_FPZ0_GRID();
    my_malloc_init(vh_global, idx, double);   


 
    for(subsystem = 0; subsystem < cei.num_subsystem; subsystem++)
    {

        NX0 = lcr[subsystem].NX_GRID *get_FG_RATIO();
        NY0 = lcr[subsystem].NY_GRID *get_FG_RATIO();
        NZ0 = lcr[subsystem].NZ_GRID *get_FG_RATIO();

        idx = NX0 * NY0 * NZ0;

        for(i = 0; i < pct.grid_npes; i++)
        {
            rcount[i] = idx/pct.grid_npes;
            rdisp[i] = (idx/pct.grid_npes) * i;
        }
        rcount[pct.grid_npes-1] += idx%pct.grid_npes;
        subsizes[0] = rcount[pct.gridpe];
        starts[0] = rdisp[pct.gridpe];

        my_malloc_init(array_tmp, idx, double);

        int order = MPI_ORDER_C;

        MPI_Type_create_subarray(1, &idx, subsizes, starts, order, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);

        MPI_Info_create(&fileinfo);
        amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
        MPI_File mpi_fhand ;

        sprintf(newname, "%s.%s", lcr[subsystem].lead_name, file_ex);
        MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

        //position = ( idx * data_indicator) * sizeof(double);
        MPI_File_read(mpi_fhand, &array_tmp[starts[0]], subsizes[0],MPI_DOUBLE, &status);

        MPI_File_close(&mpi_fhand);

        fflush(NULL);
        MPI_Allgatherv(MPI_IN_PLACE, subsizes[0], MPI_DOUBLE, array_tmp, rcount, rdisp, MPI_DOUBLE, pct.grid_comm);
        /* ================ Patches the potentials and rhos ================ */



        x0 = lcr[subsystem].x0 * get_FG_RATIO();
        y0 = lcr[subsystem].y0 * get_FG_RATIO();
        z0 = lcr[subsystem].z0 * get_FG_RATIO();
        /*if(pct.gridpe ==0) printf (" x0, y0, z0 = %d %d %d %d \n", subsystem, x0, y0, z0 );*/


        x1 = lcr[subsystem].x1 * get_FG_RATIO();
        y1 = lcr[subsystem].y1 * get_FG_RATIO();
        /*z1 = lcr[subsystem].z1 * get_FG_RATIO();
          if(pct.gridpe ==0) printf (" x1, y1, z1 = %d %d %d %d \n", subsystem, x1, y1, z1 );*/


        x2 = lcr[subsystem].x2 * get_FG_RATIO();
        y2 = lcr[subsystem].y2 * get_FG_RATIO();
        z2 = lcr[subsystem].z2 * get_FG_RATIO();
        /*if(pct.gridpe ==0) printf (" x2, y2, z2 = %d %d %d %d \n", subsystem, x2, y2, z2 );*/

        x3 = x2 + x1 - x0;
        y3 = y2 + y1 - y0;
        /*z3 = z2 + z1 - z0;
          if(pct.gridpe ==0) printf (" x3, y3, z3 = %d %d %d %d \n", subsystem, x3, y3, z3 );*/


        NYZ0 = NY0 * NZ0;
        NXZ0 = NX0 * NZ0;

        hx_old = lcr[subsystem].xside/lcr[subsystem].NX_GRID/get_FG_RATIO();
        x0_old = lcr[subsystem].x_shift + lcr[subsystem].x0 * get_FG_RATIO() * hx_old ;
        /*if(pct.gridpe ==0) printf (" x0_old, hx_old = %d %f %f \n", subsystem, x0_old, hx_old);*/

        hy_old = lcr[subsystem].yside/lcr[subsystem].NY_GRID/get_FG_RATIO();
        y0_old = lcr[subsystem].y_shift + lcr[subsystem].y0 * get_FG_RATIO() * hy_old;
        /*if(pct.gridpe ==0) printf (" y0_old, hy_old = %d %f %f \n", subsystem, y0_old, hy_old);*/

        tem = (lcr[subsystem].EF_new - lcr[subsystem].EF_old) * eV_Ha; /* update */

        if(subsystem <= 2 | subsystem > cei.num_probe) /* Satisfies left probe, central parta & right probe */
        {

            for(ix = x2; ix < x3; ix++)
            {
                for(iy = y2; iy < y3; iy++)
                {
                    idx0 = iy + ix * get_FNY_GRID(); 
                    xold_global[idx0] = x0_old + hx_old * (ix - x2);
                    yold_global[idx0] = y0_old + hy_old * (iy - y2);


                    for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                    {
                        idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                        idx = (iz + kk) + (iy - y2 + y0) * NZ0 + (ix - x2 + x0) * NYZ0; 

                        vh_global[idx1]  = array_tmp[idx] + tem * iflag;
                    }
                }
            }
        }

        else /* Satisfies down and up probes */
            /* pot/rhos_NEGF(x,y,z) = pot/rho_ON(y,x,z) */
        {
            for(ix = x2; ix < x3; ix++)
            {
                for(iy = y2; iy < y3; iy++)
                {
                    idx0 = iy + ix * get_FNY_GRID(); 
                    xold_global[idx0] = x0_old + hx_old * (ix - x2);
                    yold_global[idx0] = y0_old + hy_old * (iy - y2);

                    for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                    {
                        idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                        idx = (iz + kk) + (ix - x2 + x0) * NZ0 + (iy - y2 + y0) * NXZ0; 

                        vh_global[idx1] = array_tmp[idx] + tem * iflag;
                    }
                }
            }
        }   /* if statement ends */



        my_free(array_tmp );
    } /*  subsystem loop ends here */


    /*if(pct.gridpe ==0) printf (" Potential patching is done...  \n" );*/

    /* ====== Fillin the vaccuam space: put a ramp ============== */


    if(cei.num_probe == 4)
    {


        /*  For corner: Left-bottom */
        /*  For corner: Right-bottom */
        /*  For corner: Left-top */
        /*  For corner: Right-top */


        x2 = 0;
        x3 = lcr[3].x2 * get_FG_RATIO();
        y2 = 0;
        y3 = lcr[1].y2 * get_FG_RATIO();


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = x3 - ix;
                dis2 = y3 - iy;
                dis12 = dis1 + dis2;

                idx0 = iy + ix * get_FNY_GRID(); 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                {

                    idx1 = iz + iy * get_FPZ0_GRID() + x3 * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + y3 * get_FPZ0_GRID() + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 
                }
            }
        }



        x2 = (lcr[3].x2 + lcr[3].x1 - lcr[3].x0) * get_FG_RATIO();
        x3 = get_FNX_GRID(); 
        y2 = 0;
        y3 = lcr[2].y2 * get_FG_RATIO();


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = ix - (x2 - 1);
                dis2 = y3 - iy;
                dis12 = dis1 + dis2;

                idx0 = iy + ix * get_FNY_GRID(); 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                {

                    idx1 = iz + iy * get_FPZ0_GRID() + (x2 - 1) * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + y3 * get_FPZ0_GRID() + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }



        x2 = 0;
        x3 = lcr[4].x2 * get_FG_RATIO();
        y2 = (lcr[1].y2 + lcr[1].y1 - lcr[1].y0) * get_FG_RATIO();
        y3 = get_FNY_GRID(); 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = x3 - ix;
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * get_FNY_GRID(); 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                {

                    idx1 = iz + iy * get_FPZ0_GRID() + x3 * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * get_FPZ0_GRID() + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }



        x2 = (lcr[4].x2 + lcr[4].x1 - lcr[4].x0) * get_FG_RATIO();
        x3 = get_FNX_GRID(); 
        y2 = (lcr[2].y2 + lcr[2].y1 - lcr[2].y0) * get_FG_RATIO();
        y3 = get_FNY_GRID(); 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = ix - (x2 - 1);
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * get_FNY_GRID(); 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                {

                    idx1 = iz + iy * get_FPZ0_GRID() + (x2 - 1) * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * get_FPZ0_GRID() + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }


    }  /* end of if statement */


    /* ====== Fillin the vaccuam space: put a ramp ============== */


    if(cei.num_probe == 3)
    {


        /*  For corner: Left-top */
        /*  For corner: Right-top */


        x2 = 0;
        x3 = lcr[3].x2 * get_FG_RATIO();
        y2 = (lcr[1].y2 + lcr[1].y1 - lcr[1].y0) * get_FG_RATIO();
        y3 = get_FNY_GRID(); 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = x3 - ix;
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * get_FNY_GRID(); 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                {

                    idx1 = iz + iy * get_FPZ0_GRID() + x3 * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * get_FPZ0_GRID() + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }



        x2 = (lcr[3].x2 + lcr[3].x1 - lcr[3].x0) * get_FG_RATIO();
        x3 = get_FNX_GRID(); 
        y2 = (lcr[2].y2 + lcr[2].y1 - lcr[2].y0) * get_FG_RATIO();
        y3 = get_FNY_GRID(); 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = ix - (x2 - 1);
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * get_FNY_GRID(); 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
                {

                    idx1 = iz + iy * get_FPZ0_GRID() + (x2 - 1) * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * get_FPZ0_GRID() + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }


    }  /* end of if statement */



    /* ============= Interpolate data along x-axis =================== */

    my_malloc_init(xold,    get_FNX_GRID(), double);   
    my_malloc_init(vh_old,  get_FNX_GRID(), double);   
    my_malloc_init(vh_new,  get_FNX_GRID(), double);   


    for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
    {
        for (iy = 0; iy < get_FNY_GRID(); iy++)
        {
            for (ix = 0; ix < get_FNX_GRID(); ix++)
            {
                idx = iy + ix * get_FNY_GRID(); 
                xold[ix] = xold_global[idx]; 
                idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                vh_old[ix] = vh_global[idx1];
            }


            /* 
             *           diff_hx_interpolation2 (vh_new,  vh_old,  get_FNX_GRID(), hx_new, xold, 0.0, 0.0); 
             */

            spline(xold, vh_old, get_FNX_GRID(), 0.0, 0.0, vh_new); 


            for (ix = 0; ix < get_FNX_GRID(); ix++)
            {
                idx = iy + ix * get_FNY_GRID(); 
                xseed = hx_new * ix;  

                splint(xold, vh_old, vh_new, get_FNX_GRID(), xseed, &sfit);
                idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                vh_global[idx1]  = sfit; 

            }
        }
    }

    /*if(pct.gridpe ==0) printf (" x-interpolation is done \n" );*/


    my_free(xold);
    my_free(vh_old);
    my_free(vh_new);


    /* ============= Interpolate data along y-axis =================== */

    if(cei.num_probe >= 3)
    {

        my_malloc_init(yold,    get_FNY_GRID(), double);   
        my_malloc_init(vh_old,  get_FNY_GRID(), double);   
        my_malloc_init(vh_new,  get_FNY_GRID(), double);   


        for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
        {
            for (ix = 0; ix < get_FNX_GRID(); ix++)
            {
                for (iy = 0; iy < get_FNY_GRID(); iy++)
                {
                    idx = iy + ix * get_FNY_GRID(); 
                    yold[iy] = yold_global[idx]; 
                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_old[iy] = vh_global[idx1];
                }

                /* 
                   diff_hx_interpolation2 (vh_new,  vh_old,  get_FNY_GRID(), hy_new, yold, 0.0, 0.0); 
                 */

                spline(yold, vh_old, get_FNY_GRID(), 0.0, 0.0, vh_new); 

                for (iy = 0; iy < get_FNY_GRID(); iy++)
                {
                    idx = iy + ix * get_FNY_GRID(); 
                    yseed = hy_new * iy;  

                    splint(yold, vh_old, vh_new, get_FNY_GRID(), yseed, &sfit);
                    idx1 = iz + iy * get_FPZ0_GRID() + ix * FNYPZ; 
                    vh_global[idx1]  = sfit; 

                }
            }
        }
        /*if(pct.gridpe ==0) printf (" y-interpolation is done \n" );*/


        my_free(yold);
        my_free(vh_old);
        my_free(vh_new);

    }

    /*======================================================================*/


    for (iz = 0; iz < get_FPZ0_GRID(); iz++)  
    {
        for (ix = 0; ix < get_FPX0_GRID(); ix++)
        {
            for (iy = 0; iy < get_FPY0_GRID(); iy++)  
            {

                idx1 = iz + iy * get_FPZ0_GRID() + ix * FPYZ0;
                idx  = iz + (iy + jj) * get_FPZ0_GRID() + (ix + ii) * FNYPZ; 

                vh[idx1]  = vh_global[idx];

            }
        }
    }



    my_free(vh_global);
    my_free(xold_global);
    my_free(yold_global);



    /* ======== Move data from Global to distributive array =========== */



    MPI_Barrier(pct.img_comm);

    fflush (NULL);

}                               /* end read_data */

/* ============================================================ */


