/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

void pdb_clean (void)
{
    item_t *item;
    node_t *node;

    if (findNode ("pdb_atoms"))
    {
        node = this;
        newNode ("#tmp", NULL);
        while ((item = popItem (node)) != NULL)
        {
            if (strncmp (item->the.string, "ATOM", 4) == 0
                || strncmp (item->the.string, "HETATM", 5) == 0)
            {
                push (item);
            }
            else
            {
                my_free (item);
            }
        }
        while ((item = pop ()) != NULL)
        {
            pushItem (item, node);
        }
        killNode( unlinkNode() );
        this = node;
    }
    return;
}
