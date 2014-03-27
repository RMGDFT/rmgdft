/*
 **    $Id$    **
 ******************************************************************************/


/*
 *   Copyright (C) 2014   Wenchang Lu, Jerzy Bernholc
 *   All rights reserved.
 */

//  insert two energy points near a peak until the value difference is
//  smaller than critical_val
//  input:   
//        tot_energy_point: current number of energy points
//        cond:  value of integrand, dimension of tot_energy_point*2
//        ener1: original energy points , dimension of tot_energy_point*2
//        critical_val:  tol control
//  output:
//        energy_insert_index:  index of the inserted energy points in
//                              new array, dimension of tot_energy_point
//        ener1_temp:  the new energy points  dimension of tot_energy_point
//        EP_final:  number of new points     

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>
#include "main.h"

void find_new_energy_point_sharp_peak(double *cond, double *ener1, int
tot_energy_point, double critical_val, int *EP_final, 
                           int *energy_insert_index, double *ener1_temp)

{


    int i, EP, ione =1;

    double value;
    double *cond_work, *ener_work;

    cond_work = (double *)malloc(tot_energy_point* sizeof(double));
    ener_work = (double *)malloc(tot_energy_point* sizeof(double));


    //  copy the original value and energy point to a work array
    //  On output, the cond and ener1 will shift to the new positions
    //  with inserted energy points.

    dcopy (&tot_energy_point, cond, &ione, cond_work, &ione);
    dcopy (&tot_energy_point, ener1, &ione, ener_work, &ione);

    
    EP = 0;
   
    for(i = 1; i < tot_energy_point-1; i++)
    {
       
        value = -cond_work[i-1] + 2* cond_work[i] - cond_work[i+1];

        if( (fabs(value) > critical_val) 
            && (cond_work[i] > cond_work[i-1] )  
            && (cond_work[i] > cond_work[i+1]) ) 
        {      
        
            
            ener1_temp[EP+0] = (ener_work[i-1] + ener_work[i]) /2.0;
            ener1_temp[EP+1] = (ener_work[i] + ener_work[i+1]) /2.0;
            
            energy_insert_index[EP+0] = i + EP ;
            energy_insert_index[EP+1] = i + EP + 2;

            cond[i+EP + 1] = cond_work[i];

            ener1[i+EP + 1] = ener_work[i];

            EP += 2;
        }   
        
        else
        {
            
            cond[i+EP] = cond_work[i];

            ener1[i+EP] = ener_work[i];
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
}
