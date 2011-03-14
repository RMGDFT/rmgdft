/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

void fcatNode( FILE *stream, node_t *here )
{
    int i = 1;
    item_t *item;
    if (pct.gridpe == 0)
    {
        if (here == NULL)
        {
            Dprintf ("Non-existent node\n");
            return;
        }
        else
        {
            if (here->is == NULL)
            {
                return;
            }
            else
            {
                item = here->is;
            }
        }
        if (item->as & INFO)
            fprintf(stream, "#WARNING!!! The following tag was not validated and may be misspelled.\n");
        fprintf (stream, "%s = \"", here->name);
        if (item->as & ITEM)
        {
            while (item != NULL)
            {
                if (item->as & STR)
                {
                    if (item->the.string[0] != '#')
                    {
                        fprintf (stream, "\n");
                    }
                }
                else
                {
                    fprintf (stream, "\n");
                }
                catItem (item);
                item = item->next;
            }
            fprintf (stream, "\n\"\n");
        }
        else
        {
            catItem (item);
            fprintf (stream, "\"");
            while (item->next != NULL)
            {
                item = item->next;
                fprintf (stream, "\t");
                catItem (item);
            }
            fprintf (stream, "\n");
        }
        fflush (NULL);
    }
    return;
}
