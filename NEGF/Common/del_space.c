/************************** SVN Revision Information **************************
  **    $Id$    **
******************************************************************************/

/*****QMD-MGDFT/del_space.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2006  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 *	int del_space (FILE * fh, int *tchr, int isdata)
 * subroutine for eliminating redundant space in input 
 *    
 * INPUTS
 *   file handle , current character and isdata status (ie, quoted data).
 * OUTPUT
 *   int, true or false based on if space chars where passed over.
 * PARENTS
 *   get_input.c, listlen.c
 * CHILDREN
 * SOURCE */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "md.h"

/* subroutine for eliminating redundant space in input */
int del_space (FILE * fh, int *tchr, int isdata)
{

    int status = FALSE;

    while (isspace (*tchr) || (isdata == 0 && *tchr == '#'))
    {

        status = TRUE;

        /* never del comments in quoted data,
         * isdata only nonzero in quoted assignment */
        if (isdata == 0 && *tchr == '#')
            while ((*tchr = fgetc (fh)) != EOF)
                if (*tchr == '\n')
                    break;

        do
            if (!isspace (*tchr))
                break;
        while (((*tchr = fgetc (fh)) != EOF) && isspace (*tchr));

        if (*tchr == EOF)
            status = FALSE;

    }

    return status;
}


