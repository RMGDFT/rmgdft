#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*

     apply_potential_drop.c


Apply linear potential drop  

*/




#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "method.h"
#include "twoParts.h"

#include "pmo.h"
#include "transition.h"


#define eps 1.0E-8
#define nmax 800000000

void apply_potential_drop (double *vbias)
{
    int idx;
    int i, j, xoff;
    double V1, V2, V3, V4;
    int x1, y1, x2, y2;

    MPI_Barrier(pct.img_comm);

    xoff = get_FPX_OFFSET();

    V1 = (lcr[1].bias) * eV_Ha;
    V2 = (lcr[2].bias) * eV_Ha;

    x1 = potentialCompass.box1.x1;   
    x2 = potentialCompass.box1.x2;   
    y1 = potentialCompass.box1.y1;   
    y2 = potentialCompass.box1.y2;   


/*    my_malloc_init( vbias, get_FPX0_GRID() * get_FPY0_GRID(), double ); */

 
    /* apply a linear potential drop through conductor */

    if (cei.num_probe == 2)
    {
        for (j = 0; j < get_FPY0_GRID(); j++)
        {
            for (i = 0; i < get_FPX0_GRID(); i++)
            {
                idx = j + i * get_FPY0_GRID();
               
                if (i + xoff < x1) vbias[idx] += V1;
                else if (i + xoff > x2) vbias[idx] += V2;
                else vbias[idx] += V1 + (V2 - V1) * (i + xoff - x1)/(x2 - x1);
            }
        }
    }



    if (cei.num_probe == 4)
    {

        int iter, idx2, idx3;
        double tresh, absval; 
        double *vtemp, *vtempold;
        int FNXY_GRID;
        int nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4;


        ny1 = cei.probe_window_start[0] * get_FG_RATIO();
        ny2 = cei.probe_window_end[0] * get_FG_RATIO();
        ny3 = cei.probe_window_start[1] * get_FG_RATIO();
        ny4 = cei.probe_window_end[1] * get_FG_RATIO();
        nx1 = cei.probe_window_start[2] * get_FG_RATIO();
        nx2 = cei.probe_window_end[2] * get_FG_RATIO();
        nx3 = cei.probe_window_start[3] * get_FG_RATIO();
        nx4 = cei.probe_window_end[3] * get_FG_RATIO();


        V3 = (lcr[3].bias) * eV_Ha;
        V4 = (lcr[4].bias) * eV_Ha;

 
        if(pct.gridpe ==0)
        {
            rmg_printf (" hello0 %d %d %d %d \n", x1, x2, y1, y2);
            rmg_printf (" hello2 %d %d %d %d \n", ny1, ny2, ny3, ny4);
            rmg_printf (" hello3 %d %d %d %d \n", nx1, nx2, nx3, nx4);
            rmg_printf (" hello4 %d %d %d \n", get_FNX_GRID(), get_FNY_GRID(), get_FNZ_GRID());
        }



        FNXY_GRID = get_FNX_GRID() * get_FNY_GRID();
        my_malloc_init( vtemp, FNXY_GRID, double );
        my_malloc_init( vtempold, FNXY_GRID, double );





        /* boundary conditions */


        for (i = 0; i < get_FNX_GRID(); i++)
        {
            for (j = 0; j < get_FNY_GRID(); j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx] = 0.0; 
            }
        }


        for (i = 0; i < x1; i++)
        {
            for (j = ny1-1; j < ny2; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V1;
            }
        }


        for (i = 0; i < x1; i++)
        {
            for (j = 0; j < ny1-1; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V1*j/(ny1-2);
            }
        }


        for (i = 0; i < x1; i++)
        {
            for (j = ny2; j < get_FNY_GRID(); j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V1 - V1*(j-ny2+1)/(get_FNY_GRID()-ny2);
            }
        }


        for (i = x2-1; i < get_FNX_GRID(); i++)
        {
            for (j = ny3-1; j < ny4; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V2;
            }
        }


        for (i = x2-1; i < get_FNX_GRID(); i++)
        {
            for (j = 0; j < ny3-1; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V2*j/(ny3-2);
            }
        }


        for (i = x2-1; i < get_FNX_GRID(); i++)
        {
            for (j = ny4; j < get_FNY_GRID(); j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V2 - V2*(j-ny4+1)/(get_FNY_GRID()-ny4);
            }
        }



        for (i = nx1-1; i < nx2; i++)
        {
            for (j = 0; j < y1; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V3;
            }
        }



        for (i = 0; i < nx1-1; i++)
        {
            for (j = 0; j < y1; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V3*i/(nx1-2);
            }
        }


        for (i = nx2; i < get_FNX_GRID(); i++)
        {
            for (j = 0; j < y1; j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V3 - V3*(i-nx2+1)/(get_FNX_GRID()-nx2);
            }
        }


        for (i = nx3-1; i < nx4; i++)
        {
            for (j = y2-1; j < get_FNY_GRID(); j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V4;
            }
        }
    

        for (i = 0; i < nx3-1; i++)
        {
            for (j = y2-1; j < get_FNY_GRID(); j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V4*i/(nx3-2);
            }
        }
  


        for (i = nx4; i < get_FNX_GRID(); i++)
        {
            for (j = y2-1; j < get_FNY_GRID(); j++)
            {
                idx = j + i * get_FNY_GRID();
                vtemp[idx]+=V4 - V4*(i-nx4+1)/(get_FNX_GRID()-nx4);
            }
        }

/*
                for (i = 0; i < get_FNX_GRID() ; i++)
                {
                    for (j = 0; j < get_FNY_GRID(); j++)
                    {
                        idx = j + i * get_FNY_GRID();
                        rmg_printf (" %d %d %f \n", i, j, vtemp[idx] );
                    }
                        rmg_printf (" \n");
                }

*/

  /*Solve Laplace equation: once forward */

        for (iter = 0; iter < nmax ; iter++)
        {


            for (i = 0; i < get_FNX_GRID(); i++)
            {
                for (j = 0; j < get_FNY_GRID(); j++)
                {
                    idx = j + i * get_FNY_GRID();
                    vtempold[idx] = vtemp[idx]; 
                }
            }

/*
            for (i = x1; i < x2+1 ; i++)
            {
                for (j = y1; j < y2+1; j++)
                {
                    idx = j + i * get_FNY_GRID();
                    idx2 = j + (i-1) * get_FNY_GRID();
                    idx3 = j + (i+1) * get_FNY_GRID();
                    vtemp[idx]=(vtemp[idx2]+vtemp[idx3])/2.0;
                }
            }


            for (i = x1; i < x2+1 ; i++)
            {
                for (j = y1; j < y2+1; j++)
                {
                    idx = j + i * get_FNY_GRID();
                    idx2 = (j-1) + i * get_FNY_GRID();
                    idx3 = (j+1) + i * get_FNY_GRID();
                    vtemp[idx]=(vtemp[idx2]+vtemp[idx3])/2.0;
                }
            }

*/

            for (i = 1; i < get_FNX_GRID()-1; i++)
            {
                for (j = 1; j < get_FNY_GRID()-1; j++)
                {
                    idx = j + i * get_FNY_GRID();
                    idx2 = j + (i-1) * get_FNY_GRID();
                    idx3 = j + (i+1) * get_FNY_GRID();
                    vtemp[idx]=(vtemp[idx2]+vtemp[idx3])/2.0;
                    if((i >= nx1-1 && i < nx2) && (j < y1 || j >= y2-1)) vtemp[idx] = vtempold[idx];
                    if((j >= ny1-1 && j < ny2) && (i < x1 || i >= x2-1)) vtemp[idx] = vtempold[idx];
                }
            }
                                                                                                       
                                                                                                       
            for (i = 1; i < get_FNX_GRID()-1; i++)
            {
                for (j = 1; j < get_FNY_GRID()-1; j++)
                {
                    idx = j + i * get_FNY_GRID();
                    idx2 = (j-1) + i * get_FNY_GRID();
                    idx3 = (j+1) + i * get_FNY_GRID();
                    vtemp[idx]=(vtemp[idx2]+vtemp[idx3])/2.0;
                    if((i >= nx1-1 && i < nx2) && (j < y1 || j >= y2-1)) vtemp[idx] = vtempold[idx];
                    if((j >= ny1-1 && j < ny2) && (i < x1 || i >= x2-1)) vtemp[idx] = vtempold[idx];
                }
            }



            tresh = -1000.0;

            for (i = 0; i < get_FNX_GRID(); i++)
            {
                for (j = 0; j < get_FNY_GRID(); j++)
                {
                    idx = j + i * get_FNY_GRID();
                    absval = fabs(vtemp[idx]-vtempold[idx]);
                    if( absval > tresh) tresh = absval; 
                }
            }
            /*rmg_printf (" hello %d %f %f \n", iter, absval, tresh );*/



            if(tresh < eps)
            {



              //  if(pct.gridpe ==0)
              //  {
              //      for (i = 0; i < get_FNX_GRID(); i++)
              //      {
              //          for (j = 0; j < get_FNY_GRID(); j++)
              //          {
              //              idx = j + i * get_FNY_GRID();
              //              rmg_printf (" hello  %f \n", vtemp[idx] );
              //          }
              //      }
              // }



                global_to_distribute3 (vtemp, vbias);

                break;
            }



        } /* iter statement ends */


        if(tresh > eps)
        rmg_printf (" Did not quite converge, tresh, iter = %f %d \n", tresh, iter ); 




    my_free(vtemp);
    my_free(vtempold);

    }  /* if statement ends */




    MPI_Barrier(pct.img_comm);


}                              

