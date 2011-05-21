/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

/***** ftostr.c *****
 * NAME
 *      ftostr.c is part of the input functionality for RMG
 * COPYRIGHT
 *   Copyright (C) 2008  Frisco Rose
 *   
 * FUNCTION
 * 		Call with path/filename in char * fname. 
 * 		Returns stringified file contents (including '\0' terminator).
 */
#include <sys/stat.h>
#include "main.h"

char *filetostr (char *fname)
{
    int size, tsize;
    char *dstr;
    struct stat statbuf;
    FILE *fp;

    if (stat (fname, &statbuf) == 0)
    {
        size = (int) statbuf.st_size;

        my_malloc (dstr, size, char);
        fp = fopen (fname, "r");
        tsize = fread (dstr, 1, size - 1, fp);
        fclose (fp);
        if (tsize == size - 1)
        {
            dstr[size - 1] = '\0';
        }
    }
    else
    {
        error_handler ("Fatal Error: Unable to find essential file \"%s\"\n", fname);
    }
    return dstr;
}
