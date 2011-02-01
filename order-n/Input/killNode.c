/************************** SVN Revision Information **************************
 *  **    $Id: killNode.c 1066 2009-08-31 18:41:09Z froze $    **
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
