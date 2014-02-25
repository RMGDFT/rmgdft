/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"



/* Reads a line from the input file into tbuf skips blank lines and
   comment lines. Calls error and terminates if eof encountered. */
char *get_line (char *buf, FILE * fh)
{

    char *tptr;
    int len;

    do
    {

        if (NULL == (tptr = fgets (buf, 200, fh)))
        {

            error_handler ("Unexpected end of file");

        }                       /* end if */
        len = strlen (buf);
        if (!len)
            error_handler ("Unexpected end of file");

        /* strip the new line character at the end */
        buf[strlen (buf) - 1] = 0;

        /* skip leading whitespace */
        while ((*tptr == ' ') || (*tptr == 9))
            tptr++;

    }
    while ((*tptr == 0) || (*tptr == '#'));

    return tptr;

}                               /* end get_line */

/* Returns a pointer to the second int or real number in a
   string. Returns NULL if only one number is found in the string.  */
char *get_num (char *str)
{

    /* Skip leading white space */
    while ((*str == ' ') || (*str == 9))
        str++;
    if (!*str)
        return NULL;

    /* Find next whitespace */
    while (*str)
    {

        if ((*str == ' ') || (*str == 9))
            break;
        str++;

    }                           /* end while */

    if (!*str)
        return NULL;

    return str;

}                               /* end get_num */
