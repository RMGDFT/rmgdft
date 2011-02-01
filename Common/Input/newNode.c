/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

node_t *newNode (char *name, item_t * item)
{
    node_t tmpNode;

    Dprintf ("Make sure no node of the name %s exists", name);

    if (findNode (name))
    {
        Dprintf ("Must remove crufty node of the same name as %s", name);
        kill ();
    }

    Dprintf ("Making new keyword(TAG): %s with item", name);
    DcatItem (item);
    Dprintf ("\n");

    node_t *new;
    my_malloc (new, 1, node_t);
    new->name = name;
    new->next = new->last = new;
    new->is = item;
    return linkNode (new);
}
