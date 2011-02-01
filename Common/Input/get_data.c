/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/***** RMG/get_data.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2008  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 * 	char * get_data(char *meta, void *dest, int flags, char *data);
 *
 * 		get_data is a parsing routine that handles file input with special
 * 		attention to parsing the input settings for RMG. Returns boolean representing
 * 		success of action. get_data is data set reentrant. Current data sets are either tags
 * 		or line lists, TAGS and LINE flag respectively. Data set remains active
 * 		until called with new TAGS or LINE, INIT. Previous data sets
 * 		may be returned to if called with TAGS or LINE and same string in 
 * 		char * meta as original set INITed with.
 *
 * 		For RMG tagged input the first call should set flags as TAGS|INIT and
 * 		meta as path/input_filename. 
 * 		For line parsed input the first call should set flags as LINE|INIT and
 * 		meta as path/input_filename.
 *
 * 		char *meta;
 * 			In general meta holds descriptive information regarding data,
 * 			if flags are TAGS or LINES then meta = path/filename, else meta = tag name.
 *
 *		void *dest;
 *			Destination for type converted data according to flag type.
 *
 *		int flags;
 *			Specifies action to be performed according to the definitions below.
 *			A set is defined as the collection of keyword/value pairs
 *			A list is defined as a LIFO stack of items addressed by a keyword
 *			An item is a datum of fundamental type BOOL, INT, DBL, STR
 *
 *			INIT	==	initialize or prepare an action, if input action then
 *						"meta" contains path/filename, if "meta" previously used
 *						old data set is lost, with TAGS unduplicated tags may remain, but
 *						this should be OK since they are most likely to be default values.
 *                      Where relevent returns an integer count in dest.
 *
 *						INIT|TAGS: parse "meta" input as dictionary (keyword/value) set
 *						INIT|LINES: parse "meta" input file as newline delimited list 
 *						INIT|LIST: convert keyword's value to newline delimited list
 *						INIT|OPT:  validate specific keyword's value
 *						INIT|INFO: append message to the "meta" named data node
 *						INIT|END:  Delete the "meta" named data node 
 *						INIT: (Re)Create "meta" named data node with contents of data
 *
 *
 *			TAGS	==	Used primarily in conjuntion with INIT,
 *						also used as an adjunct flag in the definition of OPT
 *
 * 			LINES	==  Like TAGS, but uses newline parsed path/file. Declared as LIST+STR
 *
 * 			LIST	==	Used in conjunction with INIT finds keyword in meta and
 * 						splits the value into a LIFO stack with newline as delimiter.
 * 						Also used as an adjunct flag in definition of VEC and LINES
 *
 *			OPT		==	Used in conjuction with INIT - causes validation of meta's value
 *						against newline delimited list of options passed in data.
 *						Used alone retrieves validated string, most useful in verify function call.
 *						Declared as TAGS+STR
 *
 *			INFO	==	Return number of remaining item(s) in dest as INT for keyword in meta.
 *
 *			END		==	Print all keyword/values that currently exist (sortof misnomer).
 *
 * 				flag conversion types, required for tag to be considered read,
 *				read means purged from available input data unless FIX.
 * 			BOOL	==	"false" or "true" (always lower case).
 * 			INT		==	int, i.e. signed integer
 * 			DBL		==	double floating point value
 * 			VEC		==	list of 3 double floating point values
 * 			STR		==	dest is string copied from value, does not iterate.
 * 			ITEM	==	iterate through typed items, returns false when last item has been read.
 *
 *		char *data;
 *			Specifies either comment or default string to be used in absence of input
 *			and will be converted according to flag. data must be persistent (ie. string constant)
 *          if used as a default that is cast to string. If data == NULL and input not
 *			found/available then get_data returns false.
 *		
 * 			
 *		get_data supersedes and extends get_input, get_line, and get_num
 *    
 * SOURCE */

#include <string.h>
#include "main.h"

/* retrieve data from an input file.  */
bool get_data (char *meta, void *dest, flags_t flags, char *data)
{
    char *dstr;
    int size;
    if (flags & INIT)
    {
        if (flags == (INIT | TAGS))
        {
            Dprintf ("Entering TAGS|INIT for %s\n", meta);
            if (pct.thispe == 0)
            {
                newNode (meta, newItem (STR, filetostr (meta)));
                size = tagstrip ();
                dstr = this->is->the.string;
            }
            MPI_Bcast (&size, 1, MPI_INT, 0, pct.thisgrp_comm);
            if (pct.thispe != 0)
            {
                my_malloc (dstr, size, char);
            }
            MPI_Bcast (dstr, size, MPI_CHAR, 0, pct.thisgrp_comm);
            if (pct.thispe != 0)
            {
                newNode (meta, newItem (STR, dstr));
            }
			size = tagsload();
			if( dest != NULL )
                *(int *) dest = size;
            return (size > 0)? true: false;

        }
        else if (flags == (INIT | LINES))
        {
            Dprintf ("Entering LINES|INIT for %s\n", meta);
            if (pct.thispe == 0)
            {
                newNode (meta, newItem (STR, filetostr (meta)));
                size = sstripcmnt (this->is->the.string, '#' );
                dstr = this->is->the.string;
            }
            MPI_Bcast (&size, 1, MPI_INT, 0, pct.thisgrp_comm);
            if (pct.thispe != 0)
            {
                my_malloc (dstr, size, char);
            }
            MPI_Bcast (dstr, size, MPI_CHAR, 0, pct.thisgrp_comm);
            if (pct.thispe != 0)
            {
                newNode (meta, newItem (STR, dstr));
            }
			size = itemize ('\n');
			if( dest != NULL )
                *(int *) dest = size;
            return (size > 0)? true: false;

        }
        else if (flags == (INIT | LIST))
        {
            Dprintf ("Entering LIST|INIT for %s\n", meta);
            if (findNode (meta))
            {
                Dprintf ("Itemizing node %s", this->name);
				size = itemize ('\n');
				if( dest != NULL )
					*(int *) dest = size;
				return (size > 0)? true: false;
            }
            else
            {
                Dprintf ("Attempted to itemize non-existent named list.");
                return false;
            }

        }
        else if (flags == (INIT | OPT))
        {
            Dprintf ("Entering OPT|INIT for %s", meta);
            if (findNode (meta))
            {
                Dprintf ("validate existing tag against string in data");
                return validate (data);
            }
            else
            {
                newNode (meta, newItem (STR, data));
                Dprintf ("Setting future validation list.");
				size = itemize ('\n');
				if( dest != NULL )
					*(int *) dest = size;
				return (size > 0)? true: false;
            }

        }
        else if (flags == (INIT | END))
        {
            Dprintf ("Entering INIT | END for %s", meta);
            if (findNode (meta))
            {
                Dprintf ("Remove existing data node");
                return kill();
            }
            else
            {
                Dprintf ("Unable to remove non-existent data node %s.", meta);
                return false;
            }

        }
        else if (flags == (INIT | INFO))
            {
            Dprintf ("Appending (creating?) message for %s", meta);

            if (findNode (meta))
            {  
                return append (newItem (STR, data));
            }
            else
            {
                newNode (meta, newItem (STR, data));
                return append (newItem (STR, "#Important message for non-existent item!"));
            }

			Dprintf ("INIT | INFO failed for some reason.");
			return false;
        }
        else if (flags == INIT)
        {                       /* perform if and only if INIT flag set */
            Dprintf ("(Re)Creating node:%s, with contents %s\n", meta, data);
            return newNode (meta, newItem (STR, data));

        }
        Dprintf ("INIT failed for some reason.");
        return false;
    }
    else if (flags & INFO)
    {
        Dprintf ("Determining remaining items in input tag %s", meta);
        if (findNode (meta))
        {
			size = countItems ();
			if( dest != NULL )
				*(int *) dest = size;
			return (size > 0)? true: false;
        }
        else
        {
            Dprintf ("input tag %s does not seem to exist", meta);
            return false;
        }
    }
    else if (flags & END)
    {
        catNodes ();
        return true;
    }
    else
    {
        Dprintf ("Preparing action according to flags=%d for %s", flags, meta);
        if (flags & (ITEM | BOOL | INT | DBL | STR))
        {
            Dprintf ("Passed flags value matches fundamentally allowed types");
            if (findNode (meta) == false)
            {
                Dprintf ("No input of name %s found.", meta);
                if (data == NULL)
                {
                    Dprintf ("No input and no default specified!");
                    return false;
                }
                else
                {
                    Dprintf ("Setting default to: %s", data);
                    newNode (meta, newItem (STR, data));
                    append (newItem (STR, "#Was set to default!"));
                }
            }
            else
            {
                this->is->as &= ~INFO;
            }
            if (flags == OPT)
            {
                if (this->is->as & TAGS)
                {
                    Dprintf ("Option from input already validated");
                }
                else
                {
                    Dprintf ("We have to validate the default value");
                    validate (data);
                }
            }
            if (this->is->as & STR)
            {
                Dprintf ("Considering cast according to flags=%d for %s=%s\n", flags, this->name, this->is->the.string);
                cast (flags);
            }
            if (dest != NULL)
            {
                Dprintf ("Considering assignment into pointer %p as type (flag_t)=%d for %s", dest,
                         flags, this->name);
                return assign (flags, dest);
            }
        }
    }
    Dprintf ("Failed for some reason we don't know");
    return false;
}
