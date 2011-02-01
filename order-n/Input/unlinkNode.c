/************************** SVN Revision Information **************************
 *  **    $Id: unlinkNode.c 1066 2009-08-31 18:41:09Z froze $    **
 *  ******************************************************************************/


#include "main.h"

node_t *unlinkNode (void)
{
    node_t *tmp = this;
    if (this != NULL)
    {
        if (this->next == this)
        {
            Dprintf ("Unlinking nascence node\n");
            this = NULL;
        }
        else
        {
            this->next->last = this->last;
            this->last->next = this->next;
            this = tmp->last;
        }
        tmp->next = tmp->last = tmp;
    }
    return tmp;
}
