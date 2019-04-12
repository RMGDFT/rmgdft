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

#include <stdbool.h>

#if __cplusplus
extern "C" {
#endif

void nlforce_par_gamma (double * par_gamma, int ion, int nh, double *force);
void nlforce_par_omega (double * par_omega, int ion, int nh, double *force);
void nlforce_par_Q (double *veff, double *, double *, double *gamma, int ion, ION *iptr, int nh,
                     double *forces);
#if __cplusplus
}
#endif

