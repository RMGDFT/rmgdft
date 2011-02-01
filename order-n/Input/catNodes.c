/************************** SVN Revision Information **************************
 *  **    $Id: catNodes.c 1066 2009-08-31 18:41:09Z froze $    **
 *  ******************************************************************************/


#include "main.h"

void catNodes (void)
{
    int i = 0;
    node_t *here = this;
    if (here == NULL)
    {
        error_handler (" /this->next/ is uninitialized! ");
    }
    do
    {
        Dprintf ("# Node count is: %d, ", ++i);
        catNode (here);
        fflush (NULL);
        here = here->next;
    }
    while (here != this);
    return;
}
