/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

int countItems ()
{
    int i = 0;
    item_t *item;
    if (this == NULL)
    {
        error_handler ("No nodes");
    }

    item = this->is;

    while (item != NULL)
    {
        i++;
        item = item->next;
    }
    Dprintf ("Have found %d items", i);
    return i;
}
