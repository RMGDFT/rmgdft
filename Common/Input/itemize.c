/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

/* split a string in to 'delimiter' seprated items in list and
set the ITEM flag for counting */
int itemize (char delimiter)
{
    int i,k, j = 0;
    item_t *tmp = pop ();
    char *string = tmp->the.string;

    for (i = strlen (string); i >= 0; i--)
    {
        if (string[i] == delimiter)
        {
            string[i] = '\0';
            /* don't itemize for empty strings */
            if (string[i + 1] != '\0')
            {
                for(k =0; k<strlen (&string[i+1]); k++)
                {
                    if(string[i+1+k] != ' ') break;
                }
                Dprintf ("Found bounded string \n%s", &string[i + 1 + k]);
                push (newItem (STR | ITEM, &string[i + 1 + k]));
                j++;
            }
        }
    }
    if (strlen (string))
    {
        for(k =0; k<strlen (string); k++)
        {
            if(string[k] != ' ') break;
        }
        Dprintf ("Found bounded string \n%s", &string[k]);
        tmp->as |= ITEM;
        push (tmp);
    }
    else
    {
        my_free (tmp);
    }
    Dprintf ("Found %d items by splitting string on delimiter %c", j + 1, delimiter);
    return ++j;
}
