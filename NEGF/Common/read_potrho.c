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
#include "salloc.h"
#include "main.h"




void read_potrho (double *vh, int iflag, int data_indicator)
{
    long fhand;
    long nbytes, position, nbytes_first, nbytes_last;
    char newname[MAX_PATH + 200];
    char msg[200];

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
    int ii, jj, kk;
    int MaxNumBytes = 2000000000;  //maximum number of bytes one can read from ON calculation for each read
    long NumBytesLeft;  //number of bytes one need to read from ON calculation for each read

    /* Wait until everybody gets here */
    my_barrier ();

    FPYZ0 = FPZ0_GRID * FPY0_GRID;
    FNYPZ = FPZ0_GRID * FNY_GRID; 
    FNYZ = FNY_GRID * FNZ_GRID; 

    ii = pct.FPX_OFFSET;
    jj = pct.FPY_OFFSET;
    kk = pct.FPZ_OFFSET;


    hx_new = ct.xside * ct.hxxgrid;
    hy_new = ct.yside * ct.hyygrid;
    if(pct.gridpe ==0) printf (" hx_new, hy_new  =  %f  %f \n", hx_new, hy_new); 


/* =========== Allocate memory & Reading ON2 data ================ */


    idx = FNX_GRID * FNY_GRID;
    my_malloc_init(xold_global, idx, double);   
    my_malloc_init(yold_global, idx, double);   
    idx = FNX_GRID * FNY_GRID * FPZ0_GRID;
    my_malloc_init(vh_global, idx, double);   


 
    for(subsystem = 0; subsystem < cei.num_subsystem; subsystem++)
    {

        NX0 = lcr[subsystem].NX_GRID *FG_NX;
        NY0 = lcr[subsystem].NY_GRID *FG_NY;
        NZ0 = lcr[subsystem].NZ_GRID *FG_NZ;

        idx = NX0 * NY0 * NZ0;
        /*if(pct.gridpe ==0) printf (" idx +++++  =   %d \n", idx );*/

        my_malloc_init(array_tmp, idx, double);

        sprintf(newname, "%s%s", lcr[subsystem].lead_name, ".pot_rho");
        my_open( fhand, newname, O_RDWR, S_IREAD | S_IWRITE );

        position = ( idx * data_indicator) * sizeof(double);
        lseek(fhand, position, 0);
        
        NumBytesLeft = idx * sizeof(double); // this is the number of bytes we need to read at the very beginning
        
	printf ("\n  NumBytes to read at the beginning is %ld ", NumBytesLeft);
        
        array_tmp_tmp = array_tmp;  //initialize the two pointers point to the same location
        int time = 0;
        nbytes = 0;
	while (NumBytesLeft > MaxNumBytes)
	{

		nbytes_first = read (fhand, array_tmp_tmp, MaxNumBytes );
                array_tmp_tmp = array_tmp_tmp + MaxNumBytes / sizeof(double);  //shift the tmp_tmp pointer by num of doubles read already
		NumBytesLeft = NumBytesLeft - MaxNumBytes;// read the bytes left behind
		printf ("\n  NumBytes to read is %ld after the %d time", NumBytesLeft, time);
                nbytes = nbytes + nbytes_first;
		printf ("\n  NumBytes already read is %ld after the %d time", nbytes, time);
		time++;
	} 
	
	nbytes_last = read (fhand, array_tmp_tmp, NumBytesLeft ); 
	printf ("\n  NumBytes read is  %ld from the last time %d", nbytes_last, time);
	//out of the while loop, means NumBytesLeft is less than MaxNumBytes
	// read one last time to collect whatever is left from the while loop

	nbytes = nbytes + nbytes_last;


	if(nbytes != idx * sizeof(double)) 
	{
		dprintf ("\n read %ld is different from %ld ", nbytes, idx * sizeof(double));
		dprintf ("\n NX0 = %d NY0 = %d NZ0= %d subsystem = %d", NX0, NY0, NZ0, subsystem);
		error_handler ("\n Unexpected end of file vh");
	}


        close(fhand);


        /* ================ Patches the potentials and rhos ================ */



        x0 = lcr[subsystem].x0 * FG_NX;
        y0 = lcr[subsystem].y0 * FG_NY;
        z0 = lcr[subsystem].z0 * FG_NZ;
        /*if(pct.gridpe ==0) printf (" x0, y0, z0 = %d %d %d %d \n", subsystem, x0, y0, z0 );*/


        x1 = lcr[subsystem].x1 * FG_NX;
        y1 = lcr[subsystem].y1 * FG_NY;
        /*z1 = lcr[subsystem].z1 * FG_NZ;
          if(pct.gridpe ==0) printf (" x1, y1, z1 = %d %d %d %d \n", subsystem, x1, y1, z1 );*/


        x2 = lcr[subsystem].x2 * FG_NX;
        y2 = lcr[subsystem].y2 * FG_NY;
        z2 = lcr[subsystem].z2 * FG_NZ;
        /*if(pct.gridpe ==0) printf (" x2, y2, z2 = %d %d %d %d \n", subsystem, x2, y2, z2 );*/

        x3 = x2 + x1 - x0;
        y3 = y2 + y1 - y0;
        /*z3 = z2 + z1 - z0;
          if(pct.gridpe ==0) printf (" x3, y3, z3 = %d %d %d %d \n", subsystem, x3, y3, z3 );*/


        NYZ0 = NY0 * NZ0;
        NXZ0 = NX0 * NZ0;

        hx_old = lcr[subsystem].xside/lcr[subsystem].NX_GRID/FG_NX;
        x0_old = lcr[subsystem].x_shift + lcr[subsystem].x0 * FG_NX * hx_old ;
        /*if(pct.gridpe ==0) printf (" x0_old, hx_old = %d %f %f \n", subsystem, x0_old, hx_old);*/

        hy_old = lcr[subsystem].yside/lcr[subsystem].NY_GRID/FG_NY;
        y0_old = lcr[subsystem].y_shift + lcr[subsystem].y0 * FG_NY * hy_old;
        /*if(pct.gridpe ==0) printf (" y0_old, hy_old = %d %f %f \n", subsystem, y0_old, hy_old);*/

        tem = (lcr[subsystem].EF_new - lcr[subsystem].EF_old) * eV_Ha; /* update */

        if(subsystem <= 2 | subsystem > cei.num_probe) /* Satisfies left probe, central parta & right probe */
        {

            for(ix = x2; ix < x3; ix++)
            {
                for(iy = y2; iy < y3; iy++)
                {
                    idx0 = iy + ix * FNY_GRID; 
                    xold_global[idx0] = x0_old + hx_old * (ix - x2);
                    yold_global[idx0] = y0_old + hy_old * (iy - y2);


                    for (iz = 0; iz < FPZ0_GRID; iz++)  
                    {
                        idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
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
                    idx0 = iy + ix * FNY_GRID; 
                    xold_global[idx0] = x0_old + hx_old * (ix - x2);
                    yold_global[idx0] = y0_old + hy_old * (iy - y2);

                    for (iz = 0; iz < FPZ0_GRID; iz++)  
                    {
                        idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
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
        x3 = lcr[3].x2 * FG_NX;
        y2 = 0;
        y3 = lcr[1].y2 * FG_NY;


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = x3 - ix;
                dis2 = y3 - iy;
                dis12 = dis1 + dis2;

                idx0 = iy + ix * FNY_GRID; 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)  
                {

                    idx1 = iz + iy * FPZ0_GRID + x3 * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + y3 * FPZ0_GRID + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 
                }
            }
        }



        x2 = (lcr[3].x2 + lcr[3].x1 - lcr[3].x0) * FG_NX;
        x3 = FNX_GRID; 
        y2 = 0;
        y3 = lcr[2].y2 * FG_NY;


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = ix - (x2 - 1);
                dis2 = y3 - iy;
                dis12 = dis1 + dis2;

                idx0 = iy + ix * FNY_GRID; 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)  
                {

                    idx1 = iz + iy * FPZ0_GRID + (x2 - 1) * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + y3 * FPZ0_GRID + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }



        x2 = 0;
        x3 = lcr[4].x2 * FG_NX;
        y2 = (lcr[1].y2 + lcr[1].y1 - lcr[1].y0) * FG_NY;
        y3 = FNY_GRID; 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = x3 - ix;
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * FNY_GRID; 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)  
                {

                    idx1 = iz + iy * FPZ0_GRID + x3 * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * FPZ0_GRID + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }



        x2 = (lcr[4].x2 + lcr[4].x1 - lcr[4].x0) * FG_NX;
        x3 = FNX_GRID; 
        y2 = (lcr[2].y2 + lcr[2].y1 - lcr[2].y0) * FG_NY;
        y3 = FNY_GRID; 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = ix - (x2 - 1);
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * FNY_GRID; 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)  
                {

                    idx1 = iz + iy * FPZ0_GRID + (x2 - 1) * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * FPZ0_GRID + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
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
        x3 = lcr[3].x2 * FG_NX;
        y2 = (lcr[1].y2 + lcr[1].y1 - lcr[1].y0) * FG_NY;
        y3 = FNY_GRID; 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = x3 - ix;
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * FNY_GRID; 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)  
                {

                    idx1 = iz + iy * FPZ0_GRID + x3 * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * FPZ0_GRID + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }



        x2 = (lcr[3].x2 + lcr[3].x1 - lcr[3].x0) * FG_NX;
        x3 = FNX_GRID; 
        y2 = (lcr[2].y2 + lcr[2].y1 - lcr[2].y0) * FG_NY;
        y3 = FNY_GRID; 


        for (iy = y2; iy < y3; iy++)
        {
            for (ix = x2; ix < x3; ix++)
            {
                dis1 = ix - (x2 - 1);
                dis2 = iy - (y2 - 1);
                dis12 = dis1 + dis2;

                idx0 = iy + ix * FNY_GRID; 
                xold_global[idx0] = hx_new * ix;
                yold_global[idx0] = hy_new * iy;

                for (iz = 0; iz < FPZ0_GRID; iz++)  
                {

                    idx1 = iz + iy * FPZ0_GRID + (x2 - 1) * FNYPZ; 
                    V1 = vh_global[idx1];
                    idx1 = iz + (y2 - 1) * FPZ0_GRID + ix * FNYPZ; 
                    V2 = vh_global[idx1];

                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                    vh_global[idx1] = (V1 * dis2 + V2 * dis1) / dis12 * iflag; 

                }
            }
        }


    }  /* end of if statement */



    /* ============= Interpolate data along x-axis =================== */

    my_malloc_init(xold,    FNX_GRID, double);   
    my_malloc_init(vh_old,  FNX_GRID, double);   
    my_malloc_init(vh_new,  FNX_GRID, double);   


    for (iz = 0; iz < FPZ0_GRID; iz++)  
    {
        for (iy = 0; iy < FNY_GRID; iy++)
        {
            for (ix = 0; ix < FNX_GRID; ix++)
            {
                idx = iy + ix * FNY_GRID; 
                xold[ix] = xold_global[idx]; 
                idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                vh_old[ix] = vh_global[idx1];
            }


            /* 
             *           diff_hx_interpolation2 (vh_new,  vh_old,  FNX_GRID, hx_new, xold, 0.0, 0.0); 
             */

            spline(xold, vh_old, FNX_GRID, 0.0, 0.0, vh_new); 


            for (ix = 0; ix < FNX_GRID; ix++)
            {
                idx = iy + ix * FNY_GRID; 
                xseed = hx_new * ix;  

                splint(xold, vh_old, vh_new, FNX_GRID, xseed, &sfit);
                idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
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

        my_malloc_init(yold,    FNY_GRID, double);   
        my_malloc_init(vh_old,  FNY_GRID, double);   
        my_malloc_init(vh_new,  FNY_GRID, double);   


        for (iz = 0; iz < FPZ0_GRID; iz++)  
        {
            for (ix = 0; ix < FNX_GRID; ix++)
            {
                for (iy = 0; iy < FNY_GRID; iy++)
                {
                    idx = iy + ix * FNY_GRID; 
                    yold[iy] = yold_global[idx]; 
                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
                    vh_old[iy] = vh_global[idx1];
                }

                /* 
                   diff_hx_interpolation2 (vh_new,  vh_old,  FNY_GRID, hy_new, yold, 0.0, 0.0); 
                   */

                spline(yold, vh_old, FNY_GRID, 0.0, 0.0, vh_new); 

                for (iy = 0; iy < FNY_GRID; iy++)
                {
                    idx = iy + ix * FNY_GRID; 
                    yseed = hy_new * iy;  

                    splint(yold, vh_old, vh_new, FNY_GRID, yseed, &sfit);
                    idx1 = iz + iy * FPZ0_GRID + ix * FNYPZ; 
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


    for (iz = 0; iz < FPZ0_GRID; iz++)  
    {
        for (ix = 0; ix < FPX0_GRID; ix++)
        {
            for (iy = 0; iy < FPY0_GRID; iy++)  
            {

                idx1 = iz + iy * FPZ0_GRID + ix * FPYZ0;
                idx  = iz + (iy + jj) * FPZ0_GRID + (ix + ii) * FNYPZ; 

                vh[idx1]  = vh_global[idx];

            }
        }
    }



    my_free(vh_global);
    my_free(xold_global);
    my_free(yold_global);



    /* ======== Move data from Global to distributive array =========== */



    my_barrier ();

    fflush (NULL);

}                               /* end read_data */

/* ============================================================ */


