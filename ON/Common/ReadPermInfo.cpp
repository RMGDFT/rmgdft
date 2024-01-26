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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#if !(defined(_WIN32) || defined(_WIN64))
    #include <unistd.h>
#else
    #include <io.h>
#endif
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "transition.h"


/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void ReadPermInfo(char *name, unsigned int *perm_index)
{
    int fhand;
    char newname[MAX_PATH + 20];

    sprintf (newname, "%s%s", name, ".perm");

    fhand = open(newname, O_RDWR);
    if (fhand < 0)
    {
        rmg_printf("Unable to open file %s", newname);
        exit(0);
    }

    size_t nbytes = read(fhand, perm_index, ct.num_ions * sizeof(unsigned int));

    if(nbytes != (size_t)(ct.num_ions * sizeof(unsigned int)))
    {
        rmg_printf("read perminfo failed: read %zu != %zu", nbytes, ct.num_ions*sizeof(unsigned int));
        exit(0);
    }

    close (fhand);

}

