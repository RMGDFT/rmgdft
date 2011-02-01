/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include <string.h>
#include "main.h"

int tagsload (void)
{
    int i;
    char delim = '"';
    item_t *tmpItem, *varItem;
    node_t *here = this;

    i = itemize (delim);
    if (odd (i))
    {
        error_handler ("Tag/Value paring error\n");
    }
    while ((tmpItem = popItem (here)) != NULL)
    {
        Dprintf ("Creating a new node with tag/value pair from input file");
        varItem = popItem (here);
        varItem->as &= ~ITEM;
        varItem->as |= INFO;
        newNode (tmpItem->the.string, varItem);
        my_free (tmpItem);
    }
    this = here;
    here = unlinkNode ();
    my_free (here);
    return i;
}
