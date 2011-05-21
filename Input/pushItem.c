/************************** SVN Revision Information **************************
 *  **    $Id$    **
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
