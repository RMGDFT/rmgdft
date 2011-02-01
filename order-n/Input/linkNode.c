/************************** SVN Revision Information **************************
 *  **    $Id: linkNode.c 1066 2009-08-31 18:41:09Z froze $    **
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
