/************************** SVN Revision Information **************************
 *  **    $Id: sstripcmnt.c 1212 2011-01-25 18:52:03Z froze $    **
 *  ******************************************************************************/


#include <stdio.h>
#include "input.h"

/* Remove comments from a string. The comment is defined as
delim...\n and all characters after the comment char are
moved up in the string. The newline is not removed. e.g.

"All foos are bars\n but not all bars are foos. #foos, not fus.\n\0"
becomes
"All foos are bars\n but not all bars are foos. \n\0..............\0"
... indicates garbage chars in new string.

sstripcmnt returns number of chars remaining.*/

int sstripcmnt (char *string, const char delim)
{
    int i, j;
    int region = 0;

    for (i = j = 0; string[i] != '\0'; i++)
    {                           /* Clean comments. */
        if (region)
        {                       /* Are we in a comment region? */
            if (string[i] == '\n')
            {                   /* Found end of comment. */
                string[j] = '\n';
                j++;
                region = 0;
            }
        }
        else
        {
            if (string[i] == delim)
            {                   /* Found beginning of comment. */
                region = 1;
            }
            else
            {                   /* Direct copy other characters. */
                string[j] = string[i];
                j++;
            }
        }
    }                           /* Scrubed comments from input string. */
    string[j] = '\0';
    return ++j;
}
