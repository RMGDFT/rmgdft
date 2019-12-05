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
#include <boost/filesystem.hpp>
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
#include "Kpoint.h"
#include "transition.h"

void WriteXyz (char *name)
{
    char newname[MAX_PATH + 20];
    int amode;
    FILE *fhandle;

    if (pct.imgpe == 0)
    {
        std::string new_file(name);
        new_file = new_file + ".xyz";

        FILE *fhandle = fopen (new_file.c_str(), "w");
        if (!fhandle)
        {
             rmg_error_handler(__FILE__, __LINE__, "Unable to write atomic coordinate xyz file. Terminating.");
        }

        fprintf(fhandle,"%lu\n\n", Atoms.size());

        for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        {   
            ION &Atom = Atoms[ion];
            SPECIES &AtomType = Species[Atom.species];
            fprintf(fhandle,"%s %#15.12g %#15.12g %#15.12g\n", AtomType.atomic_symbol, 
                    a0_A*Atom.crds[0], a0_A*Atom.crds[1], a0_A*Atom.crds[2]);
        }


        fclose (fhandle);
        fflush(NULL);

    }

}                               /* end write_data */

