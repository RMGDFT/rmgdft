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

int open_wave_file (char *filename)
{

    char newname[MAX_PATH + 20], tmpname[MAX_PATH];
    int amode;
    int fhand;


    /* Make the new output file name */

    if (pct.spin_flag)
    {   
	if (pct.thisspin==0)
    		sprintf (newname, "%s.up%d", filename, pct.gridpe);
	else if(pct.thisspin==1) 
    		sprintf (newname, "%s.dw%d", filename, pct.gridpe);
		
    }
    else
    sprintf (newname, "%s%d", filename, pct.gridpe);

    amode = S_IREAD | S_IWRITE;


    /*PE 0 will be first one to try to open output file
     * and create directory if needed*/
    if (pct.imgpe == 0)
    {
        fhand = open (newname, O_CREAT | O_TRUNC | O_RDWR, amode);

        /*Previous call may have failed because directory did not exist
         * Let us try to to create it*/
        if (fhand < 0)
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
            my_open (fhand, newname, O_CREAT | O_TRUNC | O_RDWR, amode);

        }

    }                           /*end if (pct.gridpe == 0) */


    /*All processors should wait until 0 is done */
    //my_barrier ();
    MPI_Barrier(pct.img_comm);

    /*Other processors can now try to open */
    if (pct.imgpe)
        my_open (fhand, newname, O_CREAT | O_TRUNC | O_RDWR, amode);

    return fhand;

}
