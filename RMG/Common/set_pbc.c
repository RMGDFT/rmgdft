/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

// 
//
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "main.h"

// set periodic boundary conditions on the vector position
//
void set_pbc(double *position, int num_ions, int num_images)
{

    double xtal[3];
    int i;

    
    for(i =0; i < num_ions * num_images; i++)
    {

        to_crystal(xtal, &position[i*3]);

        if(xtal[0] >0.5 ) xtal[0] -=1.0;
        if(xtal[0] <-0.5 ) xtal[0] +=1.0;

        if(xtal[1] >0.5 ) xtal[1] -=1.0;
        if(xtal[1] <-0.5 ) xtal[1] +=1.0;

        if(xtal[2] >0.5 ) xtal[2] -=1.0;
        if(xtal[2] <-0.5 ) xtal[2] +=1.0;

        to_cartesian( xtal, &position[i*3]);
    }

}
