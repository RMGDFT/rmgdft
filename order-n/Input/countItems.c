/************************** SVN Revision Information **************************
 *  **    $Id: countItems.c 1066 2009-08-31 18:41:09Z froze $    **
 *  ******************************************************************************/


#include "main.h"

int countItems ()
{
    int i = 0;
    item_t *item;
    if (this == NULL)
    {
        Dprintf ("No nodes");
    }
    else
    {
        item = this->is;
    }
    while (item != NULL)
    {
        i++;
        item = item->next;
    }
    Dprintf ("Have found %d items", i);
    return i;
}
