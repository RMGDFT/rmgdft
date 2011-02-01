/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

item_t *newItem (flags_t type, void *value)
{
    item_t *new;
    my_malloc (new, 1, item_t);
    new->next = NULL;
    new->as = type;
    Dprintf ("Declaring new item of type %d", type);
    if (type & BOOL)
    {
        new->the.boolean = *(bool *) value;
    }
    else if (type & INT)
    {
        new->the.integer = *(int *) value;
    }
    else if (type & DBL)
    {
        new->the.rational = *(double *) value;
    }
    else if (type & STR)
    {
        new->the.string = (char *) value;
    }
    else if (type & LIST)
    {
        new->the.list = (node_t *) value;
    }
    else
    {
        error_handler ("Unkown item type");
    }
//      catItem( new ); printf( "\n" );
    return new;
}
