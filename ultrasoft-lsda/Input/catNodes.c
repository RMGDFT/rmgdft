/************************** SVN Revision Information **************************
 *  **    $Id$    **
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
