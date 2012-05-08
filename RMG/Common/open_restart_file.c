/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "main.h"



/*This opens s file for writing, returns a file handle
 * If opening the file fails, PE 0 tries to create a directory, since it is possible that
 * the reason for failure is that the directory does not exist*/

FILE *open_restart_file (char *filename)
{

    char newname[MAX_PATH + 20], tmpname[MAX_PATH];
    int amode;
    FILE *fhand;


    /* Make the new output file name */
    sprintf (newname, "%s.restart", filename);

    amode = S_IREAD | S_IWRITE;


    fhand = fopen (newname, "w");

    /*Previous call may have failed because directory did not exist
     * Let us try to to create it*/
    if (fhand == NULL)
    {

	/*Make a copy of output filename, dirname overwrites it */
	strcpy (tmpname, filename);

	printf ("\n write_data: Opening output file '%s' failed\n"
		"  Trying to create subdirectory in case it does not exist\n", newname);


	if (!mkdir (dirname (tmpname), S_IRWXU))
	    printf ("\n Creating directory %s succesfull\n\n", tmpname);
	else
	    printf ("\n Creating directory %s FAILED\n\n", tmpname);

	fflush (NULL);

	/*try opening file again */
	my_fopen (fhand, newname, "w");

    }


    return (fhand);

}
