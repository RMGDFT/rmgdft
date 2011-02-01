/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include <string.h>
#include "main.h"


bool cast (flags_t flags)
{
    int i;
    char *dstr;
    item_t *iptr;
    item_t dest;
    if (this == NULL)
    {
        Dprintf ("Trying to cast in non-existent node");
        return false;
    }
    else if (this->is == NULL)
    {
        Dprintf ("Trying to cast in non-existent item");
        return false;
    }
    if (~this->is->as & STR)
    {
        Dprintf ("Trying to cast from type %d", this->is->as);
        error_handler ("Can only cast from strings.");
    }
    else
    {
        dstr = this->is->the.string;
    }
    if (flags & STR )
    {
        Dprintf ("Cast %s=\"%s\" as STR", this->name, this->is->the.string);
        this->is->as |= STR;
        return true;
    }
    else
    {
        this->is->as &= ~STR;
        dest.as = this->is->as;
    }
    if (flags & BOOL)
    {
        if (strcmp (this->is->the.string, "true") == 0)
        {
            this->is->the.boolean = true;
        }
        else if (strcmp (this->is->the.string, "false") == 0)
        {
            this->is->the.boolean = false;
        }
        else
        {
            error_handler ("failure to match boolean input string");
        }
        this->is->as |= BOOL;
        //catNode( this );
        Dprintf ("Cast %s=\"%d\" as BOOL", this->name, this->is->the.boolean);
        return true;
    }
    else if (flags & INT)
    {
        dest.the.integer = (int)strtol (this->is->the.string, &dstr, 10);
        if (this->is->the.string == dstr)
        {
            error_handler ("In data of integer tag %s, argument :%s: is not a number!", this->name, this->is->the.string);
        }
        else
        {
            this->is->the.integer = dest.the.integer;
            this->is->as |= INT;
        }
        Dprintf ("Cast %s=\"%d\" as INT", this->name, this->is->the.integer);
        return true;
    }
    else if (flags & DBL)
    {
        if (flags & LIST)
        {                       /* Here there be a VEC */
            for (i = 0; i < 3; i++)
            {
                dest.the.rational = strtod (this->is->the.string, &dstr);
                if (this->is->the.string == dstr)
                {
                    error_handler ("In this->is->the.string of vector tag %s, argument %d is not a number!", this->name, i);
                }
                else
                {
                    this->is->the.string = dstr;
                }
                append (newItem (DBL | this->is->as, &dest.the.rational));
                Dprintf ("Cast %s=\"%f\" as VEC item %d", this->name, this->is->the.rational, i);
            }
            iptr = pop ();
            my_free (iptr);
            return true;
        }
        else
        {                       /* Here there be a scalar */
            dest.the.rational = strtod (this->is->the.string, &dstr);
            if (this->is->the.string == dstr)
            {
                error_handler ("In data of floating point tag %s, argument is not a number!", this->name);
            }
            else
            {
                this->is->as |= DBL;
                this->is->the.rational = dest.the.rational;
            }
            Dprintf ("Cast %s=\"%f\" as DBL", this->name, this->is->the.rational);
            return true;
        }
    }
    else
    {
        error_handler ("While casting to %d in tag %s, unhandled conversion case.", flags, this->name);
    }
    Dprintf ("Intended logic should never print this statement");
    return false;
}
