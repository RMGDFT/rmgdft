#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "main.h"




int assign_species (CONTROL * c, char *buf)
{
    int symbol_nonexistant = TRUE, i, sp_number = 0;


    /*Loop over species and see if symbol is known */
    for (i = 0; i < c->num_species; i++)
    {
        if (!strcmp (c->sp[i].pseudo_symbol, buf))
        {
            sp_number = i;
            symbol_nonexistant = FALSE;
            break;
        }
    }                           /*end for (i = 0; i < c->num_species; i++) */

    if (symbol_nonexistant)
    {
        printf ("\n\n PE:%d Specified symbol \"%s\" does not match any of known species symbol",
                pct.gridpe, buf);
        for (i = 0; i < c->num_species; i++)
            printf ("\n \"%s\"", c->sp[i].pseudo_symbol);

        error_handler ("Wrong specification of pseudo_symbol");
    }

    return sp_number;

}
