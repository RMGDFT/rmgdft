/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

bool killNode (node_t * node)
{
    item_t *tmpItem;

    if (node == NULL)
    {
        return false;
    }
    Dprintf ("Removing node %s from memory", node->name);
    while (node->is != NULL)
    {
        tmpItem = popItem (node);
        my_free (tmpItem);
    }
    my_free (node);
    return true;
}
