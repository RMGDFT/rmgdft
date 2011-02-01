/************************** SVN Revision Information **************************
 **    $Id$    **
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
 *   integer, true or false based on if option was set in input.
 * PARENTS
 *   
 * CHILDREN
 * SOURCE */

#include <string.h>
#include "input.h"
#include "main.h"

/* check if input agrees with provided option value */
bool verify (char *tagname, const void *optvalue)
{
    if (findNode (tagname))
    {
        Dprintf ("Verifying %s is ", tagname);
        if (optvalue == NULL)
        {
            Dprintf ("in existence, it is");
            return true;
        }
        if (this->is->as & BOOL)
        {
            if (*(bool *) optvalue == this->is->the.boolean)
            {
                Dprintf ("a BOOL that was true, as expected");
                return true;
            }
            else
            {
                Dprintf ("a BOOL that was false, not as expected");
                return false;
            }
        }
        else if (this->is->as & INT)
        {
            if (*(int *) optvalue == this->is->the.integer)
            {
                Dprintf ("an INT that was %d, as expected", this->is->the.integer);
                return true;
            }
            else
            {
                Dprintf ("an INT that was %d, not %d", this->is->the.integer, *(int *) optvalue);
                return false;
            }
        }
        else if (this->is->as & DBL)
        {
            if (*(double *) optvalue == this->is->the.rational)
            {
                Dprintf ("a DBL that was %f, as expected", this->is->the.rational);
                return true;
            }
            else
            {
                Dprintf ("a DBL that was %f, not %f", this->is->the.rational, *(double *) optvalue);
                return false;
            }
        }
        else if (this->is->as & STR)
        {
            if (strcmp (this->is->the.string, optvalue) == 0)
            {
                Dprintf ("a STR that was %s, as expected", this->is->the.string);
                return true;
            }
            else
            {
                Dprintf ("a STR that is %s, not query value of %s", this->is->the.string, (char *) optvalue);
                return false;
            }
        }
        else
        {
            error_handler  ("Tag %s has Invalid Type!", tagname);
        }
    }
    Dprintf ("not found in input");
    if (optvalue != NULL)
    /* Die only if not testing for tagname existence. */
        error_handler ("Tag %s not found in input", tagname);
    return false;
}
