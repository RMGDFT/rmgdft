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


#include <unistd.h>
#include <errno.h>
#include "main.h"
#include "common_prototypes.h"

void output_wave (STATE * states, int ik, int fhand)
{
    int status, size, is;

    /* Wait until everyone gets here */
    my_barrier ();

    status = write (fhand, &ik, sizeof (int));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kpt[0], sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kpt[1], sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kpt[2], sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    status = write (fhand, &ct.kp[ik].kweight, sizeof (rmg_double_t));
    if (status == -1)
        error_handler ("Unable to write wave state. ERRNO is %d.", errno);


    size = (GAMMA_PT) ? get_P0_BASIS() : 2 * get_P0_BASIS();

    for (is = 0; is < ct.num_states; is++)
    {
        status = write (fhand, &states[is].eig, sizeof (rmg_double_t));
        if (status == -1)
            error_handler ("Unable to write wave state. ERRNO is %d.", errno);
        status = write (fhand, states[is].psiR, size * sizeof (rmg_double_t));
        if (status == -1)
            error_handler ("Unable to write wave state. ERRNO is %d.", errno);
    }

}                               /* end write_data */

/******/
