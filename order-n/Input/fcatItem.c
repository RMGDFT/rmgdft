/************************** SVN Revision Information **************************
 *  **    $Id: fcatItem.c 1118 2010-04-19 18:14:43Z froze $    **
 *  ******************************************************************************/


#include "main.h"

void fcatItem( FILE *stream, item_t *item )
{
    if (item == NULL)
    {
        return;
    }
    if (item->as & BOOL)
    {
        if (item->the.boolean == true)
        {
            fprintf (stream, "true");
            Dprintf ("Last line displayed as BOOL");
        }
        else
        {
            fprintf (stream, "false");
            Dprintf ("Last line displayed as BOOL");
        }
    }
    else if (item->as & INT)
    {
        if (item->as & TAGS)
        {
            fprintf (stream, "#%d", item->the.integer);
            Dprintf ("Last line displayed as INT");
        }
        else
        {
            fprintf (stream, "%d", item->the.integer);
            Dprintf ("Last line displayed as INT");
        }
    }
    else if (item->as & DBL)
    {
        fprintf (stream, "%g", item->the.rational);
        Dprintf ("Last line displayed as DBL");
    }
    else if (item->as & STR)
    {
        fprintf (stream, "%s", item->the.string);
        Dprintf ("Last line displayed as STR @%p", item->the.string);
    }
    else if (item->as & LIST)
    {
        fprintf (stream, "%s", Root.name);
        Dprintf ("Last line displayed as LIST - listing not yet supported");
    }
    else
    {
        fprintf (stream, "\n\t# Unknown item type is: %d", item->as);
    }
    Dprintf ("Item was internally typed as %d", item->as);
    return;
}
