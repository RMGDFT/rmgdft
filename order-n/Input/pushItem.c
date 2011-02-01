/************************** SVN Revision Information **************************
 *  **    $Id: pushItem.c 1074 2009-09-17 18:02:05Z froze $    **
 *  ******************************************************************************/


#include "main.h"

bool pushItem (item_t * new, node_t * here)
{
    if (here != NULL)
    {
        new->next = here->is;
        here->is = new;
        Dprintf ("Pushing new item to node %s", here->name);
        DcatItem (new);
        return true;
    }
    return false;
}
