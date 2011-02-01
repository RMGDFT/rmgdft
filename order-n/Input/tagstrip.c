/************************** SVN Revision Information **************************
 *  **    $Id: tagstrip.c 1212 2011-01-25 18:52:03Z froze $    **
 *  ******************************************************************************/


#include <string.h>
#include <ctype.h>
#include "main.h"

int tagstrip (void)
{
    int i, j;
    int region = 0;
    char snip[8];
    char *string = this->is->the.string;

    /* remove comments from input string. */
    sstripcmnt (string, '#');   /* here comment is #...\n */

    if (isspace (string[0]))
        string[0] = ' ';
    for (i = j = 0; string[i] != '\0'; i++)
    {                           /* Single white space outside data, ws={"="|isspace()}. */
        if (string[i] == '"')
        {                       /* are we in a data region? */
            region++;
            region %= 2;
        }
        if (region)
        {                       /* direct copy data regions. */
            string[j] = string[i];
            j++;
            if (string[i - 1] == '"')
            {
                if (isspace (string[j - 1]))
                {
                    if (string[j - 1] == '\n')
                    {
                        j--;
                    }
                    else
                    {
                        snprintf (snip, 8, "%s", &string[i]);
                        printf
                            ("Illegal white space near ...%s...\nHope you know what your doing!\n",
                             snip);
                    }
                }
            }
        }
        else
        {                       /* Single space meta-structure region. */
            if (isspace (string[i]) || string[i] == '=')
            {
                string[j] = ' ';
                if (string[j - 1] != ' ')
                    j++;
            }
            else
            {
                string[j] = string[i];
                j++;
            }
        }
    }
    string[j] = '\0';

    for (i = j = 0; string[i] != '\0'; i++)
    {                           /* Strip spaces in meta-structure. */
        if (string[i] == '"')
        {                       /* are we in a data region? */
            region++;
            region %= 2;
        }
        if (region)
        {                       /* direct copy data regions. */
            string[j] = string[i];
            j++;
        }
        else
        {                       /* Strip space meta-structure region. */
            if (isspace (string[i]))
            {
                if (j > 0 && string[j - 1] != '"' && string[i + 1] != '"')
                {               /* Found space were none should be! */
                    snprintf (snip, 8, "%s", &string[i]);
                    printf ("Malformed input near ...%s...\nAttemtping recovery!\n", snip);
                    while (string[j - 1] != '"' && j > 0)
                        j--;
                }
            }
            else
            {
                string[j] = string[i];
                j++;
            }
        }
    }
    string[j] = '\0';

    /*
     * Free up unused memory at end. 
     */
    my_malloc (string, ++j, char);
    memcpy (string, this->is->the.string, j);
    my_free (this->is->the.string);
    this->is->the.string = string;
    return j;
}
