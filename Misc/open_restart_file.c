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

#include "portability.h"
#include <sys/types.h>
#include <sys/stat.h>
#if !(defined(_WIN32) || defined(_WIN64))
    #include <libgen.h>
#endif
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"



/*This opens s file for writing, returns a file handle
 * If opening the file fails, PE 0 tries to create a directory, since it is possible that
 * the reason for failure is that the directory does not exist*/

FILE *open_restart_file (char *filename)
{

    char newname[MAX_PATH + 20], tmpname[MAX_PATH];
    FILE *fhand;


    /* Make the new output file name */
    sprintf (newname, "%s.restart", filename);



    fhand = fopen (newname, "w");

    /*Previous call may have failed because directory did not exist
     * Let us try to to create it*/
    if (fhand == NULL)
    {

	/*Make a copy of output filename, dirname overwrites it */
	strcpy (tmpname, filename);

	printf ("\n write_data: Opening output file '%s' failed\n"
		"  Trying to create subdirectory in case it does not exist\n", newname);

#if !(defined(_WIN32) || defined(_WIN64))
	if (mkdir (dirname (tmpname), S_IRWXU))
        {
	    printf ("\n Creating directory %s FAILED\n\n", tmpname);
            rmg_error_handler("Terminating.");
        }
#else
        char dirname[_MAX_DIR];
        _splitpath(tmpname, NULL, dirname, NULL, NULL);
        if (!_mkdir(dirname))
        {
	    printf ("\n Creating directory %s FAILED\n\n", tmpname);
            rmg_error_handler("Terminating.");
        }
#endif

	printf ("\n Creating directory %s successfully\n\n", tmpname);
	fflush (NULL);

	/*try opening file again */
	my_fopen (fhand, newname, "w");

    }


    return (fhand);

}
