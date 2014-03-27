/*
 **    $Id$    **
 ******************************************************************************/


/*
 *   Copyright (C) 2014   Wenchang Lu, Jerzy Bernholc
 *   All rights reserved.
 */
//  
//  
//  find the new energy points with 5-point Simpson rule
//  input:   
//        tot_energy_point: current number of energy points
//        cond:  value of integrand, dimension of tot_energy_point*2
//        ener1: original energy points , dimension of
//        tot_energy_point*2
//        simpson_tol:  tol control
//  output:
//        energy_insert_index:  index of the inserted energy points in
//                              new array, dimension of tot_energy_point
//        ener1_temp:  the new energy points  dimension of
//        tot_energy_point
//        EP_final:  number of new points     


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>

#include "main.h"
#include "init_var.h"
#include "LCR.h"
#include "pmo.h"

double find_new_energy_point(double *cond, double *ener1, int tot_energy_point, double simpson_tol, int *EP_final, 
                           int *energy_insert_index, double *ener1_temp)

{


    int i, EP, ione =1;
    double tol, max_tol, de, sum_left, sum_right, sum_coarse;

    double *cond_work, *ener_work;


    cond_work = (double *)malloc(tot_energy_point* sizeof(double));
    ener_work = (double *)malloc(tot_energy_point* sizeof(double));


    if(tot_energy_point % 4 != 1) 
    {   
        error_handler("\n tot_energy_point is not equal n*4+1 %d", tot_energy_point);
    }

    //  copy the original conductance and energy point to a work array
    //  On output, the cond and ener1 will shift to the new positions
    //  with inserted energy points.

    dcopy (&tot_energy_point, cond, &ione, cond_work, &ione);
    dcopy (&tot_energy_point, ener1, &ione, ener_work, &ione);

    // for each five energy points, check the convergence of integration by simpson rule 
    // if not converged, insert 4 energy points
    //
    
    EP = 0;
    max_tol = 0.0;
    for(i = 0; i < tot_energy_point-1; i+=4)
    {
       
        de = ener_work[i+1] - ener_work[i]; 
        sum_left = de/3.0 * ( cond_work[i] + 4.0 * cond_work[i+1] + cond_work[i+2]);
        sum_right = de/3.0 * ( cond_work[i+2] + 4.0 * cond_work[i+3] + cond_work[i+4]);
        sum_coarse =  2.0 *de /3.0 *( cond_work[i] + 4.0 * cond_work[i+2] + cond_work[i+4]);

        tol = fabs(sum_coarse - sum_left - sum_right);
        if (tol > max_tol ) max_tol = tol;

        if (tol > simpson_tol )
        {

            ener1_temp[EP+0] = ener_work[i] + de/2.0;
            ener1_temp[EP+1] = ener_work[i+1] + de/2.0;
            ener1_temp[EP+2] = ener_work[i+2] + de/2.0;
            ener1_temp[EP+3] = ener_work[i+3] + de/2.0;
            
            energy_insert_index[EP+0] = i +EP +1;
            energy_insert_index[EP+1] = i +EP +3;
            energy_insert_index[EP+2] = i +EP +5;
            energy_insert_index[EP+3] = i +EP +7;

            cond[i+EP + 0] = cond_work[i];
            cond[i+EP + 2] = cond_work[i+1];
            cond[i+EP + 4] = cond_work[i+2];
            cond[i+EP + 6] = cond_work[i+3];

            ener1[i+EP + 0] = ener_work[i];
            ener1[i+EP + 2] = ener_work[i+1];
            ener1[i+EP + 4] = ener_work[i+2];
            ener1[i+EP + 6] = ener_work[i+3];
            EP += 4;
        }   
        
        else
        {
            
            cond[i+EP + 0] = cond_work[i];
            cond[i+EP + 1] = cond_work[i+1];
            cond[i+EP + 2] = cond_work[i+2];
            cond[i+EP + 3] = cond_work[i+3];

            ener1[i+EP + 0] = ener_work[i];
            ener1[i+EP + 1] = ener_work[i+1];
            ener1[i+EP + 2] = ener_work[i+2];
            ener1[i+EP + 3] = ener_work[i+3];
        }

    }
    
   // the last point is again the last point in new set of energy point
   // the new energy points at which conductance need to be calculated
   // are returned  in energ1_temp
   //
    ener1[tot_energy_point+ EP-1] = ener_work[tot_energy_point-1];
    cond[tot_energy_point+ EP-1] = cond_work[tot_energy_point-1];

    *EP_final = EP;
    free(cond_work);
    free(ener_work);
    return max_tol;
}
