/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

bool appendItem (item_t * new, node_t * here)
{
    item_t *tmpItem;
    if (here == NULL)
    {
        error_handler ("Appending item to NULL node doesn't make sense.");
    }
    new->as |= END;
    if ((tmpItem = here->is) == NULL)
    {
        here->is = new;
    }
    else
    {
        while (tmpItem->next != NULL)
        {
            tmpItem = tmpItem->next;
        }
        tmpItem->next = new;
    }
    Dprintf ("Appended new item to node %s", this->name);
    DcatItem (new);
    fflush (NULL);
    return true;
}
