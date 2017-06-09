/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/
 
/* 
 *     Pulay  method subroutine 
 *
 *  	
 *
 * 	Input: 	step-- iteration step for Pulay,  MUST START from ZERO 
 * 		N  ---	Dinmention of arrays rho and rho_new  
 * 		rho_new -- new charge density calculated from current wavefunctions
 * 		rho_old ---    charge density from previous step
 * 		NsavedSteps ---	Number of steps saved, including the current step, Must be greater than 2
 *              hist, rhist --- Pointers to arrays to keep history for charge density and residual
 *                              These should be setup as static double pointers:
 *                              static double **rhohist=NULL, **residhist=NULL;
 *                              Memory allocation takes place inside this function when the pointers are NULL at step 0
 * 	
 * 	Output: rho_old ----- updated 
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "common_prototypes.h"

#define 	MAX_STEPS 	160

double special_dot_product(double *a, double *b, double *b_trade, double weight, int length_x, int length_y, int length_z);

void pulay_rho(int step, int N, int N_x, int N_y, int N_z, double *rho_new, double *rho_old, int NsavedSteps, double ***hist, double ***rhist, int special_metric, double weight)
{
    double  *residual, *trade_space=NULL;
    double A[MAX_STEPS * MAX_STEPS];
    double b[MAX_STEPS];
    int ipvt[MAX_STEPS];
    int i, j;
    int ione = 1;
    int info;
    int size;
    int s2;
    double alpha;
    double t1;
    double *fi, *fj, *tptr1, *tptr2;
    int A_size;
    double scale, first_mix;

    scale = -1.0 * ct.charge_pulay_scale;

    my_malloc (residual, N, double);
    
    /*Calculate residual*/
    for (i=0; i<N; i++)
	residual[i] = rho_old [i] - rho_new [i]; 

    if (step == 0)
    {

	if (*hist == NULL)
	{
	    my_malloc(*hist, NsavedSteps - 1, double *);
	    my_malloc((*hist)[0], N * (NsavedSteps - 1), double );

	    for (i=1; i<NsavedSteps - 1; i++)
		(*hist)[i] = (*hist)[0] + i*N;
	}
	
	if (*rhist == NULL)
	{
	    my_malloc(*rhist, NsavedSteps - 1, double *);
	    my_malloc((*rhist)[0], N * (NsavedSteps - 1), double );

	    for (i=1; i<NsavedSteps - 1; i++)
		(*rhist)[i] = (*rhist)[0] + i*N;
	}
    

        /* For the first step, mixing quite a bit of new density is helpful, 
         * proabably becasue initial guess can be pretty bad*/
        first_mix = ct.mix;
        if (first_mix < 0.5)
            first_mix = 0.5;
    
	/* For first step do linear mixing only*/
        alpha = - first_mix;
	for (i=0; i<N; i++)
	    rho_old [i]  = rho_old[i] + alpha * residual [i];
        
	/*Add current residual and charge density into history*/
	QMD_dcopy(N, rho_new, ione, (*hist)[0], ione);
        QMD_dcopy(N, residual, ione, (*rhist)[0], ione);
    }
    else
    {
        size = step + 1;
        if (size > NsavedSteps)
            size = NsavedSteps;
        A_size = size + 1;
        s2 = A_size * A_size;

	if (s2 >= MAX_STEPS) 
	    error_handler("Pulay order is too high");

	/*Open memory space for trade_image, which will only be used for special metrics*/
	if (special_metric) my_malloc(trade_space, (N_x+2) * (N_y+2) * (N_z+2), double);

        for (i = 0; i < MAX_STEPS * MAX_STEPS; i++)
            A[i] = 0.0;
        for (i = 0; i < size; i++)
        {
            for (j = i; j < size; j++)
            {
                if (i == (size - 1))
                    fi = residual;
		else
		    fi = (*rhist)[i];
                
                if (j == (size - 1))
                    fj = residual;
		else
		    fj = (*rhist)[j];

                /*  get A(i,j)  */
                A[j * (size + 1) + i] = QMD_ddot(N, fi, ione, fj, ione);

		/* Evaluate dot product with special metric as used by gpaw, if requested*/
		if (special_metric) 
		    A[j * (size + 1) + i] = (1.0 + weight/8.0) *  A[j * (size + 1) + i] + 
			special_dot_product(fi, fj, trade_space, weight, N_x, N_y, N_z); 
                
		A[i * (size + 1) + j] = A[j * (size + 1) + i];
            }

            b[i] = 0.0;
        }
	
	if (special_metric) my_free(trade_space);
        
	/*  Sum A over all processors*/
        global_sums(A, &s2, pct.grid_comm);

        b[size] = 1.0;
        for (i = 0; i < size; i++)
        {
            A[i * (size + 1) + size] = 1.0;
            A[size * (size + 1) + i] = 1.0;
        }


        /*   b = A^(-1) * b     */
        info = 0;
        dgesv(&A_size, &ione, A, &A_size, ipvt, b, &A_size, &info);
        //double work[MAX_STEPS*20];
        //int lwork = MAX_STEPS*20;
        //dgels("N", &A_size, &A_size, &ione, A, &A_size, b, &A_size, work, &lwork, &info);
	/*Print pulay constants*/
#if 0
        printf ("\n\nPulay Constants:");
	for (i = 0; i < size; i++)
	{
	    if (i % 5 == 0)
		printf("\n");
	    
	    printf("b[%d]:%.4f  ", i, b[i]);
	}
#endif

	/*Now that the vector b[i] is determined we determine new charge density as
	 * a linear combination of previous and current density and residuals
	 * b[i] gives weights for the linear combination */

	/*Scale rho history by b[i], store result in rho_old*/
	for (i=0; i<N; i++)
	    rho_old[i] = b[0] * (*hist)[0][i];

	for (i = 1; i < size - 1; i++)
	    QMD_daxpy(N, b[i], (*hist)[i], ione, rho_old, ione);

	QMD_daxpy(N, b[size - 1], rho_new, ione, rho_old, ione);
	    
	/* Scale residual history by -1 * b[i] * scale and add to the result stored in rho_old*/
	for (i = 0; i < size - 1; i++)
	{
	    t1 = -1.0 * scale * b[i];
	    QMD_daxpy(N, t1, (*rhist)[i], ione, rho_old, ione);
	}

	t1 = -1.0 * scale * b[size - 1];
	QMD_daxpy(N, t1, residual, ione, rho_old, ione);


	if (step > (NsavedSteps - 2))
	{
	    /*Rotate history for charge density and residual*/
	    tptr1 = (*hist)[0];
	    tptr2 = (*rhist)[0];
	    for (i = 0; i < size - 2; i++)
	    {
		(*hist)[i] = (*hist)[i+1];
		(*rhist)[i] = (*rhist)[i+1];
	    }

	    (*hist)[size - 2] = tptr1;
	    (*rhist)[size - 2] = tptr2;

	    /*Add new charge density and residual into history*/
	    QMD_dcopy(N, rho_new, ione, (*hist)[size - 2], ione);
	    QMD_dcopy(N, residual, ione, (*rhist)[size - 2], ione);
	}

	/*If we do not have all history saved, add current values there
	 * without rotation and overwriting*/
	else
	{
	    QMD_dcopy(N, rho_new,  ione, (*hist)[(size - 1)], ione);
	    QMD_dcopy(N, residual, ione, (*rhist)[(size - 1)], ione);
	}

    }

    my_free(residual);

}


double special_dot_product(double *a, double *b, double *b_trade, double weight, int length_x, int length_y, int length_z)
{
    int ix, iy, iz;
    int idx, idy, idz;
    int idxp1, idyp1, idzp1;
    int idxm1, idym1, idzm1;
    double rvalue = 0.0, a_value;
    int incx, incy;
   
    trade_imagesx (b, b_trade, length_x, length_y, length_z, 1, FULL_TRADE);

    incy = length_z + 2;
    incx = incy * (length_y + 2);

    for (ix=1; ix < length_x+1; ix++)
    {
	idx = ix*incx;
	
	idxp1 = ix + 1;
	idxp1 *= incx;

	idxm1 = ix - 1;
	idxm1 *= incx;
	
	for (iy=1; iy < length_y+1; iy++)
	{
	    idy = iy*incy;
		
	    idyp1 = iy + 1;
	    idyp1 *= incy;
		
	    idym1 = iy - 1;
	    idym1 *= incy;
	    
	    for (iz=1; iz < length_z+1; iz++)
	    {

		idz = iz; 

		idzp1 = iz + 1;
		idzm1 = iz - 1;

		a_value = a[(ix -1) * length_y * length_z + (iy -1) * length_z + iz -1];

		/*First NN*/
		rvalue += a_value * (weight / 16.0) *
		        (b_trade[idxp1 + idy + idz]
		       + b_trade[idxm1 + idy + idz]
		       + b_trade[idx + idyp1 + idz]
		       + b_trade[idx + idym1 + idz]
		       + b_trade[idx + idy + idzp1]
		       + b_trade[idx + idy + idzm1]);

		rvalue += a_value * (weight / 32.0) *
		         (b_trade[idxp1 + idyp1 + idz]
		        + b_trade[idxp1 + idym1 + idz]
		        + b_trade[idxp1 + idy + idzp1]
		        + b_trade[idxp1 + idy + idzm1]
		        + b_trade[idxm1 + idyp1 + idz]
		        + b_trade[idxm1 + idym1 + idz]
		        + b_trade[idxm1 + idy + idzp1]
		        + b_trade[idxm1 + idy + idzm1]
			+ b_trade[idx + idyp1 + idzp1]
			+ b_trade[idx + idyp1 + idzm1]
			+ b_trade[idx + idym1 + idzp1]
			+ b_trade[idx + idym1 + idzm1]);
		
		rvalue += a_value * (weight / 64.0) *
		        (b_trade[idxp1 + idyp1 + idzp1]
		        + b_trade[idxp1 + idyp1 + idzm1]
		        + b_trade[idxp1 + idym1 + idzp1]
		        + b_trade[idxp1 + idym1 + idzm1]
		        + b_trade[idxm1 + idyp1 + idzp1]
		        + b_trade[idxm1 + idyp1 + idzm1]
		        + b_trade[idxm1 + idym1 + idzp1]
		        + b_trade[idxm1 + idym1 + idzm1]);

	    }
	}
    }

		
    return rvalue;

}


