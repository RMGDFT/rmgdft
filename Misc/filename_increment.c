/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/
#include "portability.h"
#include <sys/types.h>
#include <sys/stat.h>
#if !(defined(_WIN32) || defined(_WIN64))
    #include <libgen.h>
#else
    #define S_ISDIR(ST_MODE) (((ST_MODE) & _S_IFMT) == _S_IFDIR)
#endif
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "main.h"
#include "macros.h"

int filename_increment(char *filename)
{

    int lognum = 0, status;
    int pathlength;
#if (defined(_WIN32) || defined(_WIN64))
    char workdir[_MAX_PATH], basename[_MAX_PATH], dirname[_MAX_DIR];
    _splitpath(filename, NULL, dirname, NULL, NULL);
    pathlength = _MAX_PATH;
#else
    char workdir[MAX_PATH], basename[MAX_PATH];
    pathlength = MAX_PATH;
#endif
    struct stat buffer;


    strncpy(basename, filename, sizeof(basename));


#if (defined(_WIN32) || defined(_WIN64))
    if (status = stat ( basename, &buffer ) == 0)
    {
#else
    if (status = stat ( dirname(basename), &buffer ) == 0)
    {
#endif
        if ( !S_ISDIR( buffer.st_mode))
            error_handler 
                ("Found %s, that is not a directory ", basename);
    }
    else
    {
#if (defined(_WIN32) || defined(_WIN64))
        if (_mkdir (basename))
#else
        if ((status = mkdir (basename, S_IRWXU)) == -1)
#endif
            error_handler
                ("Unable to create directory \"%s\", mkdir returned %d, check file system for permissions/diskfull.",
                 basename, status);
    }


    snprintf (basename, pathlength, "%s.%02d", filename, lognum );
    while ((status = stat (basename, &buffer)) != -1)
    {
        if (++lognum > 99)
            error_handler
                ("You have over 100 logfiles, you need to think of a better job scenario!\n");
        snprintf (basename, pathlength, "%s.%02d", filename, lognum );
    }

    /* open and save logfile handle, printf is stdout before here */
    return lognum;

}
