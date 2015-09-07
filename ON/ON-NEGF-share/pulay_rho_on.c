/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/* 
 *     Pulay  method subroutine 
 *     modified by Wenchang Lu 10-06-2005
 *      
 *
 *  Input:  step-- iteration step for Pulay,  MUST START from 0.0 
 *      N  ---  Dinmention of array  xm and fm  
 *      xm ---  current sollution 
 *      fm ---  current residual 
 *      NsavedSteps --- Number of steps saved, 
 *           including the current step, Must greater than 2
 *      Nrefresh:  refresh the pulay mixing each Nrefresh steps
 *
 *  Output: xm ----- updated 
 * 
 *      Note:  No D here, since preconditioning is done outside of this subroutine  
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "prototypes_on.h"
#include "BaseGrid.h"
#include "transition.h"

#define     MAX_STEPS   300


void pulay_rho_on (int step0, int N, double *xm, double *fm, int NsavedSteps,
                int Nrefresh, double scale, int preconditioning)
{
    static double *x;
    static double *f;
    static int step_effect = 0;
    double A[MAX_STEPS * MAX_STEPS];
    double b[MAX_STEPS];
    double b0[MAX_STEPS], bmin, bmax;
    int ipvt[MAX_STEPS];
    int i, j;
    int ione = 1;
    int info;
    int size, step;
    int s2;
    double one = 1.0;
    double alpha;
    double t1;
    double *x1, *x2, *f1, *f2;
    double *fi, *fj;
    int A_size;

    int incx, incy, dimx, dimy, dimz;
    int ix, iy ,iz;
    int ixs, ixms, ixps, iys, iyps, iyms;
    double *res_s, cc, fc, ec, crn;


    cc = 1.0/27.0;
    fc =  1.0/27.0;
    ec =  1.0/27.0;
    crn =  1.0/27.0;
 //   dimx = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
 //   dimy = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
 //   dimz = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);
    dimx = get_FPX0_GRID();
    dimy = get_FPY0_GRID();
    dimz = get_FPZ0_GRID();
    incx = (dimy+2) * (dimz+2);
    incy = dimz+2;

    res_s = (double *) malloc((dimx+2)*(dimy+2)*(dimz+2)*sizeof(double));

    pack_ptos(res_s, fm, dimx, dimy, dimz);
    trade_images(res_s, dimx, dimy, dimz, FULL_TRADE);
    
    for(ix = 1; ix <=dimx; ix++)
    {
        ixs = ix * incx;
        ixms = (ix - 1) * incx;
        ixps = (ix + 1) * incx;

        for (iy = 1; iy <= dimy; iy++)
        {

            iys = iy * incy;
            iyms = (iy - 1) * incy;
            iyps = (iy + 1) * incy;

            for (iz = 1; iz <= dimz; iz++)
            {

                fm[(ix-1) * dimy * dimz + (iy-1)* dimz + iz-1]
                   = cc * res_s[ixs + iys + iz] +
                    fc * (res_s[ixms + iys + iz] +
                          res_s[ixps + iys + iz] +
                          res_s[ixs + iyms + iz] +
                          res_s[ixs + iyps + iz] +
                          res_s[ixs + iys + (iz - 1)] +
                          res_s[ixs + iys + (iz + 1)]) +
                    ec * (res_s[ixms + iys + (iz + 1)] +
                          res_s[ixps + iys + (iz + 1)] +
                          res_s[ixms + iys + (iz - 1)] +
                          res_s[ixps + iys + (iz - 1)] +
                          res_s[ixs + iyms + (iz + 1)] +
                          res_s[ixs + iyps + (iz + 1)] +
                          res_s[ixs + iyms + (iz - 1)] +
                          res_s[ixs + iyps + (iz - 1)] +
                          res_s[ixms + iyms +  iz    ] +
                          res_s[ixps + iyms +  iz    ] +
                          res_s[ixms + iyps +  iz    ] +
                          res_s[ixps + iyps +  iz    ]) +
                    crn *(res_s[ixms + iyms + (iz - 1)] +
                          res_s[ixms + iyms + (iz + 1)] +
                          res_s[ixms + iyps + (iz - 1)] +
                          res_s[ixms + iyps + (iz + 1)] +
                          res_s[ixps + iyms + (iz - 1)] +
                          res_s[ixps + iyms + (iz + 1)] +
                          res_s[ixps + iyps + (iz - 1)] +
                          res_s[ixps + iyps + (iz + 1)]);
            }
        }
    }


     
    pack_stop(res_s, fm, dimx, dimy, dimz);

    step = step_effect % Nrefresh;

    printf("\n  step_effect  %d  %d ", step_effect, step);
    if (NsavedSteps == 1)
    {
        alpha = -ct.mix;
        daxpy (&N, &alpha, fm, &ione, xm, &ione);
        return;
    }

/*  allocate memory for previous steps */
    if (step == 0)
    {
        if(x == NULL) my_malloc( x, N * (NsavedSteps - 1), double );
        if(f == NULL)my_malloc( f, N * (NsavedSteps - 1), double );
        if ((x == NULL) || (f == NULL))
        {
            printf ("pulay.c ---Could not allocate memory for x or f \n");
            fflush (NULL);
            exit (-1);
        }
    }

    if (step == 0)
    {
        dcopy (&N, xm, &ione, x, &ione);
        dcopy (&N, fm, &ione, f, &ione);

        alpha = -ct.mix;
        daxpy (&N, &alpha, fm, &ione, xm, &ione);


    }
    else
    {
        size = step + 1;
            if (size > NsavedSteps)
                size = NsavedSteps;
            A_size = size + 1;
            s2 = A_size * A_size;

            for (i = 0; i < MAX_STEPS * MAX_STEPS; i++)
                A[i] = 0.0;
            for (i = 0; i < size; i++)
            {
                for (j = i; j < size; j++)
                {
                    fi = f + i * N;
                    fj = f + j * N;
                    if (i == (size - 1))
                        fi = fm;
                    if (j == (size - 1))
                        fj = fm;

                    /*  get A(i,j)  */
                    A[j * (size + 1) + i] = ddot (&N, fi, &ione, fj, &ione);
                    A[i * (size + 1) + j] = A[j * (size + 1) + i];
                }

                b[i] = 0.0;
            }
            /*  Real_sum_all ( A, b )  if mutil-processing */
            global_sums (A, &s2, pct.grid_comm);
            global_sums (b, &A_size, pct.grid_comm);

            b[size] = 1.0;
            for (i = 0; i < size; i++)
            {
                A[i * (size + 1) + size] = 1.0;
                A[size * (size + 1) + i] = 1.0;
            }



            /*   b = A^(-1) * b     */
            dgesv (&A_size, &ione, A, &A_size, ipvt, b, &A_size, &info);

        bmin = 1000.0;
        bmax = -1000.0;
        for(i = 0; i < size; i++) 
        {
            bmin = rmg_min(bmin, b[i]);
            bmax = rmg_max(bmax, b[i]);
        }
            


#if 0
        bmin = -1.5;
        for (i = 0; i < size; i++) 
        {
            if(b[i] < bmin)
            {
                for(j = 0; j <size; j++)
                {
                    b[j] = b[j] * (1.0 + (b0[i]-bmin)/(1.0-b0[i]));
                }
                b[i] = bmin;
                break;
            }
        }
#endif


        if (pct.gridpe == 0)
        {
            printf ("\n");
            for (i = 0; i < size; i++)
                printf ("   pulay_b[%d]: %10.6f  \n", i, b[i]);
        }


        if (step <= (NsavedSteps - 2))
        {
            x1 = x + (size - 1) * N;
            f1 = f + (size - 1) * N;
            dcopy (&N, xm, &ione, x1, &ione);
            dcopy (&N, fm, &ione, f1, &ione);

            t1 = b[size - 1];
            dscal (&N, &t1, xm, &ione);
            for (i = 0; i < size - 1; i++)
            {
                x1 = x + i * N;
                daxpy (&N, &b[i], x1, &ione, xm, &ione);
            }

            t1 = 0.0;
            dscal (&N, &t1, fm, &ione);
            for (i = 0; i < size; i++)
            {
                t1 = -1.0 * b[i];
                f1 = f + i * N;
                daxpy (&N, &t1, f1, &ione, fm, &ione);
            }


            t1 = scale;
            daxpy (&N, &t1, fm, &ione, xm, &ione);

        }
        else
        {
            t1 = b[0];
            dscal (&N, &t1, x, &ione);
            for (i = 1; i < size - 1; i++)
            {
                x1 = x + i * N;
                daxpy (&N, &b[i], x1, &ione, x, &ione);
            }
            daxpy (&N, &b[size - 1], xm, &ione, x, &ione);

            t1 = -1.0 * b[0];
            dscal (&N, &t1, f, &ione);
            for (i = 1; i < size - 1; i++)
            {
                t1 = -1.0 * b[i];
                f1 = f + i * N;
                daxpy (&N, &t1, f1, &ione, f, &ione);
            }
            t1 = -1.0 * b[size - 1];
            daxpy (&N, &t1, fm, &ione, f, &ione);


            t1 = scale;
            dscal (&N, &t1, f, &ione);
            daxpy (&N, &one, x, &ione, f, &ione);


            for (i = 0; i < size - 2; i++)
            {
                x1 = x + i * N;
                x2 = x + i * N + N;
                dcopy (&N, x2, &ione, x1, &ione);
            }
            x1 = x + (size - 2) * N;
            dcopy (&N, xm, &ione, x1, &ione);

            dcopy (&N, f, &ione, xm, &ione);

            for (i = 0; i < size - 2; i++)
            {
                f1 = f + i * N;
                f2 = f + i * N + N;
                dcopy (&N, f2, &ione, f1, &ione);
            }
            f1 = f + (size - 2) * N;
            dcopy (&N, fm, &ione, f1, &ione);
        }

    }
    step_effect++;
}
