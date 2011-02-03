/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/***** RMG: Common/init_IO.c *****
 * NAME
 *   Ab initio real space multigrid acceleration
 *   Quantum molecular dynamics package.
 * COPYRIGHT
 *   Copyright (C) 2009  Frisco Rose, Jerzy Bernholc
 * FUNCTION
 *   void init_IO( int argc, char **argv )
 *   Initializes settings and creates directory structures for ouput logging
 *   Make each run image manage its own directory of input/output
 * INPUTS
 *   argc and argv from main
 * OUTPUT
 *   none
 * PARENTS
 *   main.c
 * CHILDREN
 *   init_pe.c read_pseudo.c
 * SOURCE
 */



#include "main.h"
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>


void init_IO (int argc, char **argv)
{

    int status, lognum = 0;
    char workdir[MAX_PATH], logname[MAX_PATH], basename[MAX_PATH], *quantity, *extension, *endptr;
    struct stat buffer;
    time_t timer;

    /* Set start of program time */
    timer = time (NULL);
    /* Initialize MPI, we need it for error_handler, amongst others */
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &pct.thispe);
    pct.thisgrp = pct.thispe / NPES;

    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;



    /* 1st argument must exist as a file, use it as the control file or die */
    if ((status = stat (argv[1], &buffer)) == -1)
        error_handler ("Command line argument \"%s\" is malformed, file D.N.E.", argv[1]);
    strncpy (ct.cfile, argv[1], MAX_PATH);

    /* test to see if argument 1 determines a multi-image run, fail gracefully if not */
    pct.images = 1;
    strncpy (basename, argv[1], MAX_PATH);
    if ((extension = rindex (basename, '.')) != NULL)
        *extension++ = '\0';
    if ((quantity = rindex (basename, '.')) != NULL)
    {
        *quantity++ = '\0';
        pct.images = (int) strtol (quantity, &endptr, 10);

        /*Check if this is *just* a number, not some numerical naming convention */
        if (quantity == endptr)
            pct.images = 1;
    }

    Dprintf ("RMG will run with %d images", pct.images);

    /* PE(MPI) initialization, need mpi groups defined before logfile is initialized */
    init_pe ();

    if (pct.images > 1)
    {
        if (strcmp (extension, "rmg") != 0)
            error_handler
                ("Multi-image input file %s does not end with proper \"rmg\" extension.", argv[1]);

        /* logfile name is based on input file and this images group number */
        /* if second command line argument exists, use it as a basename for output */
        if (argc == 3)
        {
            snprintf (logname, MAX_PATH, "%s.%d.log", argv[2], pct.thisimg + 1);
        }
        else
        {
            snprintf (logname, MAX_PATH, "%s.%d.log", basename, pct.thisimg + 1);
        }

        /* After chdir, this groups control file will be in parent directory */
        snprintf (ct.cfile, MAX_PATH, "../%s.%d.rmg", basename, pct.thisimg + 1);

        /* every image has it own output/working directory */
        sprintf (workdir, "image.%d", pct.thisimg + 1);
        if (pct.imgpe == 0)
        {
            if (status = stat (workdir, &buffer) == 0)
            {
                if (!S_ISDIR (buffer.st_mode))
                    error_handler
                        ("Found %s, that is not a directory as required for multi-image run!",
                         workdir);
            }
            else
            {
                if ((status = mkdir (workdir, S_IRWXU)) == -1)
                    error_handler
                        ("Unable to create image directory \"%s\", mkdir returned %d, check file system for permissions/diskfull.",
                         workdir, status);
            }
        }
        /* Make sure that process 0 has successfully verified/created the image directory */
        my_barrier ();

        if ((status = chdir (workdir)) == -1)
            error_handler
                ("Unable to cd to image directory \"%s\", check file system for permissions/diskfull.",
                 workdir);
    }
    else
    {
        /* Not multi-image run, so set output file accordingly */
        /* if second command line argument exists, use it as a basename for output */
        if (argc == 3)
        {
            snprintf (logname, MAX_PATH, "%s.log", argv[2]);
        }
        else
        {
            /* Reset basename to be input file name (except rmg extension) */
            strncpy (basename, argv[1], MAX_PATH);
            if ((extension = rindex (basename, '.')) != NULL)
            {
                *extension++ = '\0';
                if (strcmp (extension, "rmg") == 0)
                    snprintf (logname, MAX_PATH, "%s.log", basename);
                else
                    snprintf (logname, MAX_PATH, "%s.%s.log", basename, extension);
            }
            else
            {
                snprintf (logname, MAX_PATH, "%s.log", basename);
            }
        }
    }

    /* if logname exists, increment until unique filename found */
    if (pct.imgpe == 0)
    {
        while ((status = stat (logname, &buffer)) != -1)
        {
            strncpy (basename, logname, MAX_PATH);
            if (++lognum > 99)
                error_handler
                    ("You have over 100 logfiles, you need to think of a better job scenario!\n");
            if ((extension = rindex (basename, '.')) != NULL)
                *extension = '\0';
            if ((quantity = rindex (basename, '.')) != NULL)
                if (quantity[1] == 'r')
                    *quantity = '\0';
            snprintf (logname, MAX_PATH, "%s.r%02d.log", basename, lognum);
        }

        /* open and save logfile handle, printf is stdout before here */
        my_fopen (ct.logfile, logname, "w");
    }

    MPI_Comm_size (pct.grid_comm, &status);
    printf ("\nRMG run started at GMT %s", asctime (gmtime (&timer)));
    printf ("\nRMG running in message passing mode with %d procs.\n", status);

    /* Read in our pseudopotential information */
    read_pseudo ();

    return;
}
