/************************** SVN Revision Information **************************
 **    $Id: read_states_info.c 1242 2011-02-02 18:55:23Z luw $    **
******************************************************************************/
 
/*


*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "md.h"



void read_states_info (char *name, STATE * states)
{
    int fhand;
    int st;
    char newname[MAX_PATH + 200];
    unsigned nbytes;

    my_barrier ();

    sprintf (newname, "%s%s", name, ".states_info");
    fhand = open (newname, O_RDWR);
    if (fhand < 0)
        error_handler (" Unable to open file ");

    for (st = 0; st < ct.num_states; st++)
    {
        nbytes = read (fhand, &states[st].pe, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].crds[0], sizeof (double));
        if (nbytes != sizeof (double))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].crds[1], sizeof (double));
        if (nbytes != sizeof (double))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].crds[2], sizeof (double));
        if (nbytes != sizeof (double))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].radius, sizeof (double));
        if (nbytes != sizeof (double))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].movable, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].frozen, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].index, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].ixmin, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].iymin, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].izmin, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].ixmax, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].iymax, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].izmax, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].xfold, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].yfold, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].zfold, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].ixstart, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].iystart, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].izstart, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].ixend, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].iyend, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].izend, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].orbit_nx, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].orbit_ny, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].orbit_nz, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

        nbytes = read (fhand, &states[st].size, sizeof (int));
        if (nbytes != sizeof (int))
            error_handler ("Unexpected end of file");

    }

}
