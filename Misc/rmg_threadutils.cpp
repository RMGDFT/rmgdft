/*
 *
 * Copyright 2025 The RMG Project Developers. See the COPYRIGHT file 
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

#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "transition.h"


int rmg_get_active_threads(void)
{
    int active_threads = ct.MG_THREADS_PER_NODE;
    if(ct.mpi_queue_mode && (active_threads > 1)) active_threads--;
    return active_threads;
}

int rmg_get_max_threads(void)
{
    return ct.MG_THREADS_PER_NODE;
}
