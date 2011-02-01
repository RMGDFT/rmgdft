/************************** SVN Revision Information **************************
 *  **    $Id: pdb_clean.c 1066 2009-08-31 18:41:09Z froze $    **
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
        kill ();
        this = node;
    }
    return;
}
