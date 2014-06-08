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
#include "macros.h"

int filename_increment(char *filename)
{

    int lognum = 0, status;
    char workdir[MAX_PATH], basename[MAX_PATH];
    struct stat buffer;


    strcpy(basename, filename);


    if (status = stat ( dirname(basename), &buffer ) == 0)
    {
        if ( !S_ISDIR( buffer.st_mode))
            error_handler 
                ("Found %s, that is not a directory ", basename);
    }
    else
    {
        if ((status = mkdir (basename, S_IRWXU)) == -1)
            error_handler
                ("Unable to create directory \"%s\", mkdir returned %d, check file system for permissions/diskfull.",
                 basename, status);
    }



    snprintf (basename, MAX_PATH, "%s.%02d", filename, lognum );
    while ((status = stat (basename, &buffer)) != -1)
    {
        if (++lognum > 99)
            error_handler
                ("You have over 100 logfiles, you need to think of a better job scenario!\n");
        snprintf (basename, MAX_PATH, "%s.%02d", filename, lognum );
    }

    /* open and save logfile handle, printf is stdout before here */
    return lognum;
}
