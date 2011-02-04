
/************************** SVN Revision Information **************************
  **    $Id: verify.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/

/*****QMD-MGDFT/verify.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2006  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 *	int verify (char *tagname, char *optname)
 *  check if input agrees with provided option value 
 *    
 * INPUTS
 *   name of the tag to check and the option value to check against.
 * OUTPUT
 *   int, true or false based on if option was set in input.
 * PARENTS
 *   
 * CHILDREN
 * SOURCE */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "md.h"

extern struct t_list tag_list[];
extern struct b_list bool_list[];
extern int ntags;
extern int nbools;

/* check if input agrees with provided option value */
int verify (char *tagname, char *optname)
{

    int i;

    for (i = 0; i < ntags; i++)
    {
        if (strcmp (tag_list[i].TagName, tagname) == 0)
        {

            if (strcmp (optname, tag_list[i].OptName) == 0)
                return TRUE;
            else
                return FALSE;

        }
    }

    for (i = 0; i < nbools; i++)
    {
        if (strcmp (bool_list[i].TagName, tagname) == 0)
        {

            if ((strcmp (optname, "true") == 0) && bool_list[i].Flag)
                return TRUE;
            else
                return FALSE;

        }
    }

    printf ("get_input.c: Tag %s not found in inputs.h", tagname);
    error_handler ("Tag Not Found");
    /* This never happens, but it suppresses warnings. */
    return FALSE;
}
