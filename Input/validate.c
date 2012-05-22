/************************** SVN Revision Information **************************
 *  **    $Id$    **
 *  ******************************************************************************/


#include "main.h"

bool validate (char *optlist)
{
    int zero = 0;
    item_t *validation_count = newItem (INT | TAGS, &zero);
    item_t *tmpItem, *varItem;
    if (this->is == NULL || optlist == NULL)
    {
        printf ("Validation of NULL is nonsense\n");
        return false;
    }

    Dprintf("flagged as %d\n", this->is->as);
    if (this->is->as & TAGS)
    {
        Dprintf ("Already INIT validated tag, must be a duplicate call");
        return true;
    }

    if (countItems () == 1)
    {
        Dprintf ("Validating from input generated node against validation list in optlist");
        node_t *here = this;
        newNode ("#tmp", newItem (STR, optlist));
        itemize ('\n');
        while (this->is != NULL)
        {
            if (strcmp (this->is->the.string, here->is->the.string))
            {
                validation_count->the.integer++;
                tmpItem = pop ();
                my_free (tmpItem);
            }
            else
            {
                Dprintf ("Validated %s against input %s", this->is->the.string,
                         here->is->the.string);
                kill ();
                this = here;
                this->is->as |= TAGS;
                append (validation_count);
                return true;
            }
        }
        error_handler("Fatal unmatched input %s=\"%s\",\nCheck for (in)valid %s options in read_control.c", here->name, here->is->the.string, here->name);
        return false;
    }
    else
    {
        Dprintf
            ("Validating default in optlist against validation list in node made from INIT|OPT");
        while (this->is != NULL)
        {
            tmpItem = pop ();
            //              Dprintf("Validating %s against default %s", tmpItem->the.string, optlist );
            if (strcmp (tmpItem->the.string, optlist))
            {
                validation_count->the.integer++;
                my_free (tmpItem);
            }
            else
            {
                while (this->is != NULL)
                {
                    varItem = pop ();
                    my_free (varItem);
                }
                Dprintf ("Validated %s against default %s", tmpItem->the.string, optlist);
                tmpItem->the.string = optlist;
                tmpItem->as &= ~ITEM;
                tmpItem->as |= TAGS;
                push (tmpItem);
                append (validation_count);
                append (newItem (STR, "#Was set to default!"));
                return true;
            }
        }
    }
    Dprintf ("Validation failure with %s trying to match %s", this->name, optlist);
    error_handler ("Not able to validate option.\n");
    return false;               /* should never get here */
}
