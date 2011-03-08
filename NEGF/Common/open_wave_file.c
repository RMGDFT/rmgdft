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
#include "md.h"



/*This opens s file for writing, returns a file handle
 * If opening the file fails, PE 0 tries to create a directory, since it is possible that
 * the reason for failure is that the directory does not exist*/

int open_wave_file (char *filename)
{

    char newname[MAX_PATH + 20], tmpname[MAX_PATH];
    int amode;
    int fhand;

    amode = S_IREAD | S_IWRITE;

    fhand = open (filename, O_CREAT | O_TRUNC | O_RDWR, amode);

    /*Previous call may have failed because directory did not exist
     * Let us try to to create it*/
    if (fhand < 0)
    {

	/*Make a copy of output filename, dirname overwrites it*/
	strcpy(tmpname, filename);

	printf( "\n write_data: Opening output file '%s' failed\n" 
		"  Trying to create subdirectory in case it does not exist\n", 
		filename );


	if (!mkdir(dirname(tmpname),S_IRWXU))
	    printf ("\n Creating directory '%s' succesful\n\n", dirname(tmpname));
	else
	    printf ("\n Creating directory '%s' FAILED\n\n", dirname(tmpname));

	fflush (NULL);

	/*try opening file again */
	fhand = open (filename, O_CREAT | O_TRUNC | O_RDWR, amode);

    }

    return fhand;

}
