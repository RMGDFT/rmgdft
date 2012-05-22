/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

bool assign (flags_t flags, void *dest)
{

    Dprintf ("Reached assign for %s of type %d requested as %d", this->name, this->is->as, flags);
    if (this->is->as & END)
    {
        Dprintf(" here is where to reset END on all items so that reread of list can be performed.");
        item_t *item = this->is;
        while (item != NULL)
        {
            item->as &= ~END;
            item = item->next;
        }

        Dprintf ("Reached the end of unread items in %s, returning false", this->name);
        dest = NULL;
        return false;
    }
    else
    {
        if (this->is->as & TAGS)
        {
            *(int *) dest = this->is->next->the.integer;
            Dprintf ("Granted TAGS second item INT assignment %d from %s",
                     this->is->next->the.integer, this->name);
            return true;
        }
        if (this->is->as & BOOL)
        {
            *(bool *) dest = this->is->the.boolean;
            Dprintf ("Granted BOOL assignment %d from %s", this->is->the.boolean, this->name);
        }
        else if (this->is->as & INT)
        {
            *(int *) dest = this->is->the.integer;
            Dprintf ("Granted INT assignment %d from %s", this->is->the.integer, this->name);
        }
        else if (this->is->as & DBL)
        {
            *(double *) dest = this->is->the.rational;
            Dprintf ("Granted DBL assignment %f from %s", this->is->the.rational, this->name);
        }
        else if (flags & (STR | ITEM) && this->is->as & (STR | ITEM))
        {
            Dprintf ("Copying %d characters into %p", strlen (this->is->the.string), dest);
            strcpy ((char *) dest, this->is->the.string);
            Dprintf ("Granted STR assignment %s from %s", this->is->the.string, this->name);
        }
    }
    if (this->is->as & ITEM)
    {
        /* The items of a list dont go away they just
         * move to the bottom of the stack. */
        Dprintf ("Top item being demoted to bottom item");
        return append (pop ());
    }
    return true;
}
