/************************** SVN Revision Information **************************
 **    $Id: get_input.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/

/*****QMD-MGDFT/get_input.c *****
 * NAME
 *   Ab initio real space code with multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2005  Frisco Rose
 *                       Jerzy Bernholc
 *   
 * FUNCTION
 * 	int get_input(FILE *fh, const char *id, void *dest, int flag, const char *def);
 * 		get_input rewinds the file refered to by *fh and searches for a
 *		tag that matches char *id, it returns the number of characters of data
 *		in the tags argument on success. The data returned in void *dest 
 *		is whitespace substituted with a single space ' ' and converted into
 *		void *dest depending on what int flag is
 *		
 *			SEQ		== to rewind or not to rewind?, this is the answer.
 * 			LIST	== reads in sequence of lines from data of tag
 * 			
 * 			The following are variations of numeric types
 * 			BOOL	== "false" or "true" (allways lower case) listed in inputs.h
 * 			INT		== int, i.e. signed int
 * 			DBL		== double floating point value
 * 			VEC		== 3 double floating point values
 * 			
 * 			The following are variations of char* (of length not to exceed MAX_CHARS)
 * 			STR		== generic string
 * 			OPT		== set tag option (string) according to definitions in inputs.h
 * 			LINE	== reads in a white space compressed line
 * 			
 *		char *def represents the default value to load into dest if no element
 *		matching char *id is found and will be converted according to flag.
 *		if int flag is  "SEQ" then the file is not rewound, 
 *		instead it continues to read from the last position.
 *		get_input superceeds and extends get_line, and get_num
 *    
 * INPUTS
 *   main control structure
 * OUTPUT
 *   variables in structure CONTROL c are updated
 *   in most of other file, the name is ct.... 
 *   see md.h for structure CONTROL
 * PARENTS
 *   md.c
 * CHILDREN
 * SOURCE */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "inputs.h"

/* get_input subroutine for setting Flag of option[j] in tag[i] */
static void set_opt (char *tagname, int *dest, char *optname)
{

    int i, j;

    for (i = 0; i < ntags; i++)
    {
        if (strcmp (tag_list[i].TagName, tagname) == 0)
        {

            for (j = 0; j < tag_list[i].OptCount; j++)
            {

                if (strcmp (tag_list[i].Opt[j], optname) == 0)
                {
                    *dest = j;
                    tag_list[i].OptName = tag_list[i].Opt[j];
                    return;
                }

            }
            printf
                ("get_input.c: Attempting to set option %s of tag %s, option not found in inputs.h\n",
                 optname, tagname);
            if (pct.thispe == 0)
            {
                printf ("Possible tags of option '%s':\n", optname);
                for (j = 0; j < tag_list[i].OptCount; j++)
                    printf ("  '%s'\n", tag_list[i].Opt[j]);
            }
            error_handler ("Option Not Found");

        }
    }

    for (i = 0; i < nbools; i++)
    {

        if (strcmp (bool_list[i].TagName, tagname) == 0)
        {
            if (strcmp (optname, "true") == 0)
                *dest = bool_list[i].Flag = TRUE;
            else if (strcmp (optname, "false") == 0)
                *dest = bool_list[i].Flag = FALSE;
            else
            {
                printf ("get_input.c: Setting boolean %s of tag %s, it's neither true or false",
                        optname, tagname);
                error_handler ("Boolean neither true or false");
            }
            return;
        }

    }
    printf ("get_input.c: Tag %s not found in inputs.h", tagname);
    error_handler ("Tag Not Found");
}

#define del_space() del_space(fh, &tchr, isdata)
#define isnumber(arg) ((arg) == '+' || (arg) == '-' || (arg) == '.' || isdigit((arg))) ? TRUE : FALSE

/* retrieve data from an input file.  */
int get_input (FILE * fh, char *id, void *dest, unsigned int flag, char *def)
{

    int ischar = 0, isdata = 0, tchr;
    char buf[MAX_CHAR], *tdata, *tptr, *eptr;
    static long islist;
    tptr = buf;
    tdata = tptr;

    /* using LIST as first get_input from file probably broken */
    if (!(flag & LIST))
        islist = 0;

    for (tchr = 0; tchr < MAX_CHAR; tchr++)
        tdata[tchr] = '\0';

    if (flag & LINE)
    {
        /* LINE should almost always be called as flag = LINE|SEQ
         * (think about it) */
        del_space ();

        do
        {
            while (isspace (tchr) && (tchr != '\n') && (tchr != EOF) && (flag & RAW))
            {
                tdata[isdata] = ' ';
                tchr = fgetc (fh);
            }

            if (tchr == '\n' || tchr == '#' || tchr == EOF)
                break;

            tdata[isdata] = tchr;

            if (++isdata > MAX_CHAR)
                error_handler ("LINE data exceeds MAX_CHAR length.");
        }
        while ((tchr = fgetc (fh)) != EOF);

    }
    else
    {

        while ((tchr = fgetc (fh)) != EOF)
        {

            if (del_space ())
                /* The first char of a tag must follow a space, see the
                 * fputc after rewind */
                if ((ischar == 0) && (id[ischar] == tchr))
                {
                    tchr = fgetc (fh);
                    ischar = 1;
                }

            if (isdata)
            {

                if (flag & LIST)
                {
                    if (islist == 0)
                    {
                        if (pct.thispe == 0)
                            printf ("%s=\n\"", id);
                    }
                    else
                    {
                        fseek (fh, islist, SEEK_SET);
                        tchr = fgetc (fh);
                        if (!(flag & RAW))
                            del_space ();
                    }
                    do
                    {
                        if (tchr == '"' || tchr == '\n' || tchr == '#')
                        {
                            if (tchr == '#')
                                while ((tchr != '\n') && (tchr != EOF))
                                    tchr = fgetc (fh);
                            if (--isdata)
                                islist = ftell (fh) - 1;
                            else if (pct.thispe == 0)
                                printf ("\n\"\n");
                            goto exit_routine;
                        }
                        tdata[isdata - 1] = tchr;

                        if (++isdata >= MAX_CHAR)
                        {
                            /* This incorrectly fires if next fgetc is '\"' */
                            printf ("While getting data for tag %s,", id);
                            error_handler ("Data exceeds buffer size of MAX_CHAR length.");
                        }
                    }
                    while ((tchr = fgetc (fh)) != EOF);
                    printf ("While getting data for tag %s,", id);
                    error_handler ("EOF before closing quote.");

                }
                else
                {
                    do
                    {
                        if (!(flag & RAW))
                        {
                            if (del_space ())
                            {
                                tdata[isdata - 1] = ' ';
                                if (++isdata >= MAX_CHAR)
                                {
                                    /* This incorrectly fires if next fgetc is '\"' */
                                    printf ("While getting data for tag %s,", id);
                                    error_handler ("Data exceeds buffer size of MAX_CHAR length.");
                                }
                            }
                        }
                        if (tchr == '"')        /* found all data, this was closing quote. */
                            goto exit_routine;  /* break doubly nested while */


                        tdata[isdata - 1] = tchr;

                        if (++isdata >= MAX_CHAR)
                        {
                            /* This incorrectly fires if next fgetc is '\"' */
                            printf ("While getting data for tag %s,", id);
                            error_handler ("Data exceeds buffer size of MAX_CHAR length.");
                        }
                    }
                    while ((tchr = fgetc (fh)) != EOF);
                }
                printf ("While getting data for tag %s,", id);
                error_handler ("EOF reached before closing quote.");
            }

            if (ischar && id[ischar] == '\0' && isdata == 0)
            {
                /* found word == id[] */
                if (tchr == '=')
                {
                    /* found sequentially id[] and equal sign, almost tag match. */
                    if ((tchr = fgetc (fh)) == EOF)
                    {
                        printf ("While getting data for tag %s,", id);
                        error_handler ("EOF reached before opening quote.");
                    }

                    del_space ();

                    if (tchr == '"')    /* found sequentially id[] and equal sign and quote, tag match! */
                        isdata = 1;
                    else
					{
						if ( pct.thispe == 0 ) {
							printf ("While reading options for tag %s, leading quote on option missing.", id);
						}
						error_handler ("Encountered malformed option.");
					}
						

                }
                else
                    ischar = 0;

            }

            /* only continue finding tag if id has matched at least first ischar */
            if (ischar && id[ischar] == tchr)
                /* found matching character */
                ischar++;
            else
                /* nonmatch, reset */
                ischar = 0;

        }
    }

  exit_routine:
    tptr = tdata;

    if (isdata == 0)
    {
        /* always rewind on zero data for next call */
        rewind (fh);
        fputc (' ', fh);
        islist = 0;

        if (def != NULL)
        {
            tdata = def;
            if (pct.thispe == 0)
                printf ("#DEFAULT\t%s=\"", id);
        }
        else if (!(flag & SEQ) && !(flag & LIST))
        {
            /* If doing sequential or list reads just report if zero data */
            printf ("get_input.c: Required tag %s not found in input file", id);
            error_handler ("EOF before required tag");
        }

    }

    if (dest != NULL)
    {

        if (isdata && !islist)
            if (pct.thispe == 0)
                printf ("%s=\"", id);

        if (!(flag & SEQ) || isdata != 0 || def != NULL)
        {

            if ((flag & BOOL) || (flag & OPT))
            {

                set_opt (id, dest, tdata);
                if (pct.thispe == 0)
                    printf ("%s\"\n", tdata);

            }
            else if (flag & INT)
            {

                if (isnumber (tdata[0]))
                    *(int *) dest = strtol (tdata, NULL, 10);
                else
                {
                    printf ("In data of int tag %s,", id);
                    error_handler ("Argument is not a number!");
                }
                if (pct.thispe == 0)
                    printf ("%i\"\n", *(int *) dest);

            }
            else if (flag & DBL)
            {

                if (isnumber (tdata[0]))
                    *(REAL *) dest = atof (tdata);
                else
                {
                    printf ("In data of floating point tag %s,", id);
                    error_handler ("Argument is not a number!");
                }
                if (pct.thispe == 0)
                    printf ("%g\"\n", *(REAL *) dest);

            }
            else if (flag & VEC)
            {

                ((REAL *) dest)[0] = strtod (tdata, &eptr);
                if (tdata == eptr)
                {
                    printf ("In data of vector tag %s,", id);
                    error_handler ("Argument 1 is not a number!");
                }
                else
                    tdata = eptr;
                if (pct.thispe == 0)
                    printf ("%10.6f", ((REAL *) dest)[0]);

                ((REAL *) dest)[1] = strtod (tdata, &eptr);
                if (tdata == eptr)
                {
                    printf ("In data of vector tag %s,", id);
                    error_handler ("Argument 2 is not a number!");
                }
                else
                    tdata = eptr;
                if (pct.thispe == 0)
                    printf ("  %10.6f", ((REAL *) dest)[1]);

                ((REAL *) dest)[2] = strtod (tdata, &eptr);
                if (tdata == eptr)
                {
                    printf ("In data of vector tag %s,", id);
                    error_handler ("Argument 3 is not a number!");
                }
                else
                    tdata = eptr;
                if (pct.thispe == 0)
                    printf ("  %10.6f\"\n", ((REAL *) dest)[2]);

            }
            else if ((flag & STR) || (flag & LINE))
            {
                strcpy ((char *) dest, tdata);
                if (pct.thispe == 0)
                    printf ("%s\"\n", tdata);

            }
            else if (flag & LIST)
            {
                strcpy ((char *) dest, tdata);
                if (pct.thispe == 0 && isdata)
                {
                    if (flag & RAW)
                        printf ("\n%s", tdata);
                    else
                        printf ("\n\t%s", tdata);
                }
            }
        }

    }
    else
    {
        printf ("While converting data for tag %s,", id);
        error_handler ("Cannot load data in to Null dest pointer.");
    }

    if (!(flag & SEQ))
    {
        /* rewind for next call unless SEQ */
        rewind (fh);
        fputc (' ', fh);
    }

    fflush (NULL);
    return isdata;
    /* how much data was obtained in number of chars */
}
