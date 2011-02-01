/************************** SVN Revision Information **************************
  **    $Id: listlen.c 809 2007-05-29 23:15:00Z miro $    **
******************************************************************************/

/*****QMD-MGDFT/listlen.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2006  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 * 	int listlen(FILE *fh, const char *id);
 * 		computes the number of list entries for a given file hndle and list name.		
 *    
 * INPUTS
 *   file handle and list name.
 * OUTPUT
 *   integer, number of lines
 * PARENTS
 *   md.c
 * CHILDREN
 * SOURCE */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "md.h"

#define del_space() del_space(fh, &tchr, isdata)

/* retrieve data from an input file.  */
int listlen (FILE * fh, char *id)
{

    int ischar = 0, isdata = 0, count = 0, tchr;

    rewind (fh);
    fputc (' ', fh);

    while ((tchr = fgetc (fh)) != EOF)
    {

        if (del_space ())
        {
            /* The first char of a tag must follow a space, see the
             * fputc after rewind */
            if ((ischar == 0) && (id[ischar] == tchr))
            {
                tchr = fgetc (fh);
                ischar = 1;
            }
            if (tchr == '#')
                while ((tchr != '\n') && (tchr != EOF))
                    tchr = fgetc (fh);
        }

        if (isdata)
        {
            tchr = fgetc (fh);

            do {
                if (tchr == '"' || tchr == '\n' || tchr == '#')
                {
                    if (tchr == '"')
                        goto exit_routine;
                    if (tchr == '#')
                        while ((tchr != '\n') && (tchr != EOF))
                            tchr = fgetc (fh);
                    count++;
					del_space();
					fseek (fh, -1, SEEK_CUR);
                }
            } while ((tchr = fgetc (fh)) != EOF);

            printf ("While getting data for tag %s,", id);
            error_handler ("EOF before closing quote.");
        }

        if (ischar && id[ischar] == '\0' && isdata == 0)
        {
            /* found word == id[] */
            if (tchr == '=')
            {
                /* found sequentially id[] and equal sign, almost tag match. */
                if ((tchr = fgetc (fh)) == EOF)
                {
                    printf ("While getting data for tag %s,", id);
                    error_handler ("EOF reached before opening quote.");
                }

                del_space ();

                if (tchr == '"')        /* found sequentially id[] and equal sign and quote, tag match! */
                    isdata = 1;
                else
                    ischar = 0;

            }
            else
                ischar = 0;

        }

        /* only continue finding tag if id has matched at least first ischar */
        if (ischar && id[ischar] == tchr)
            /* found matching character */
            ischar++;
        else
            /* nonmatch, reset */
            ischar = 0;

    }
  exit_routine:

    rewind (fh);
    return count;
    /* how much data was obtained in number of lines */
}
