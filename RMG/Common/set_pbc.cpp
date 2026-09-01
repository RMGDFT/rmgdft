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
