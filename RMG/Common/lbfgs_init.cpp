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

#include "lbfgs.h"

void lbfgs_init(int num_coeff)
{


    int item;

    memory_steps = 20;
    finite_step = 0.005;
    maxmove = 0.2;

    lineopt = 0;

    invcurv = 0.01;

    //my_malloc_init(ro, memory_steps, double);
    //my_malloc_init(alpha_lbfgs, memory_steps, double);
    ro = new double[memory_steps]();
    alpha_lbfgs = new double[memory_steps]();

    item = num_coeff * memory_steps;

    //my_malloc_init(change_in_G, item, double);
    //my_malloc_init(change_in_R,  item,double);
    change_in_G = new double[item]();
    change_in_R = new double[item]();

    item = num_coeff;
    //my_malloc_init(Rold, item, double);
    //my_malloc_init(Fold, item, double);
    //my_malloc_init(direction, item, double);
    Rold = new double[item]();
    Fold = new double[item]();
    direction = new double[item]();

    reset_flag = 1;
    fdstep = 1;
}

