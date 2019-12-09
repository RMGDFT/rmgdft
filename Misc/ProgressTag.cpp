
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


#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "RmgException.h"
#include "Kpoint.h"
#include "transition.h"
#include "Plots.h"
#include "Functional.h"
#include "Exxbase.h"
#include "../Headers/macros.h"


void ProgressTag(double step_time, double elapsed_time)
{
    if (pct.imgpe == 0) {
        rmg_printf (" quench: [md: %3d/%-d  scf: %3d/%-d  step time: %6.2f  scf time: %8.2f secs  RMS[dV]: %8.2e ]\n\n\n",
                ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps, step_time, elapsed_time, ct.rms);

        /*Also print to stdout*/
        if(pct.images == 1)
            fprintf (stdout,"\n quench: [md: %3d/%-d  scf: %3d/%-d  step time: %6.2f  scf time: %8.2f secs  RMS[dV]: %8.2e ]",
                    ct.md_steps, ct.max_md_steps, ct.scf_steps, ct.max_scf_steps, step_time, elapsed_time, ct.rms);
    }
}

void ExxProgressTag(double step_time, double elapsed_time)
{
    if (pct.imgpe == 0) {
        rmg_printf (" quench: [md: %3d/%-d  exx: %3d/%-d  step time: %6.2f  scf time: %8.2f secs  EXX[dV]: %8.2e ]\n\n\n",
                ct.md_steps, ct.max_md_steps, ct.exx_steps, ct.max_exx_steps, step_time, elapsed_time, ct.exx_delta);

        /*Also print to stdout*/
        if(pct.images == 1)
            fprintf (stdout,"\n quench: [md: %3d/%-d  exx: %3d/%-d  step time: %6.2f  scf time: %8.2f secs  EXX[dV]: %8.2e ]",
                    ct.md_steps, ct.max_md_steps, ct.exx_steps, ct.max_exx_steps, step_time, elapsed_time, ct.exx_delta);
    }
}
