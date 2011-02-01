/************************** SVN Revision Information **************************
 *  **    $Id: popItem.c 1074 2009-09-17 18:02:05Z froze $    **
 *  ******************************************************************************/


#include "main.h"

item_t *popItem (node_t * here)
{
    if (here != NULL)
    {
        item_t *tmp = here->is;
        if (here->is != NULL)
        {
            here->is = here->is->next;
            tmp->next = NULL;
            tmp->as &= ~END;
            Dprintf ("Popping item from node %s", here->name);
            DcatItem (tmp);
            Dprintf ("\n");
        }
        return tmp;
    }
    return NULL;
}
