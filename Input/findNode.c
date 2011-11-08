/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

bool findNode (char *name)
{
    node_t *here = this;
    if (this == NULL)
    {
        return false;
    }
    Dprintf("this is currently named %s", here->name);
    do
    {
        //Dprintf( "#Looking at %s for %s\n", here->name, name );
        if (strcmp (name, here->name))
        {
            here = here->next;
        }
        else
        {
            Dprintf ("Found %s of type %d", here->name, here->is->as);
            this = here;
            return true;
        }
    }
    while (this != here);
    Dprintf ("Unable to find %s", name);
    return false;
}
