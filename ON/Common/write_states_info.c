/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*


	write_states_info.c

    Functions to write data to files.


*/



#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"



void write_states_info(char *name, STATE * states)
{
    int amode;
    char newname[MAX_PATH + 20];
    int st;
    int fhand;


    my_barrier();

    if (pct.gridpe == 0)
    {
        sprintf(newname, "%s%s", name, ".states_info");

        amode = S_IREAD | S_IWRITE;

        fhand = open(newname, O_CREAT | O_TRUNC | O_RDWR, amode);
        if (fhand < 0)
            error_handler(" Unable to write file ");

        for (st = 0; st < ct.num_states; st++)
        {
            write(fhand, &states[st].pe, sizeof(int));
            write(fhand, &states[st].crds[0], sizeof(double));
            write(fhand, &states[st].crds[1], sizeof(double));
            write(fhand, &states[st].crds[2], sizeof(double));
            write(fhand, &states[st].radius, sizeof(double));
            write(fhand, &states[st].movable, sizeof(int));
            write(fhand, &states[st].frozen, sizeof(int));
            write(fhand, &states[st].index, sizeof(int));
            write(fhand, &states[st].ixmin, sizeof(int));
            write(fhand, &states[st].iymin, sizeof(int));
            write(fhand, &states[st].izmin, sizeof(int));
            write(fhand, &states[st].ixmax, sizeof(int));
            write(fhand, &states[st].iymax, sizeof(int));
            write(fhand, &states[st].izmax, sizeof(int));
            write(fhand, &states[st].xfold, sizeof(int));
            write(fhand, &states[st].yfold, sizeof(int));
            write(fhand, &states[st].zfold, sizeof(int));
            write(fhand, &states[st].ixstart, sizeof(int));
            write(fhand, &states[st].iystart, sizeof(int));
            write(fhand, &states[st].izstart, sizeof(int));
            write(fhand, &states[st].ixend, sizeof(int));
            write(fhand, &states[st].iyend, sizeof(int));
            write(fhand, &states[st].izend, sizeof(int));
            write(fhand, &states[st].orbit_nx, sizeof(int));
            write(fhand, &states[st].orbit_ny, sizeof(int));
            write(fhand, &states[st].orbit_nz, sizeof(int));
            write(fhand, &states[st].size, sizeof(int));
        }
    }


}
