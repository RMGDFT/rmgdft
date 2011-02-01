/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

node_t *linkNode (node_t * new)
{
    if (this != NULL)
    {
        this->next->last = new;
        new->next = this->next;
        this->next = new;
        new->last = this;
    }
    return this = new;
}
