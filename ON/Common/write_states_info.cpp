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
#include "prototypes_on.h"
#include "write_wrapper.h"



void write_states_info(char *name, STATE * states)
{
    int amode;
    int st;
    int fhand;


    MPI_Barrier(pct.img_comm);

    std::string newname;
    if (pct.gridpe == 0)
    {
        newname = std::format("{}{}", name, ".states_info");

        amode = S_IREAD | S_IWRITE;

        fhand = open(newname.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
        if (fhand < 0)
            rmg::error(" Unable to write file ");

        for (st = 0; st < ct.num_states; st++)
        {
            write_int(fhand, &states[st].pe, 1);
            write_double(fhand, &states[st].crds[0], 1);
            write_double(fhand, &states[st].crds[1], 1);
            write_double(fhand, &states[st].crds[2], 1);
            write_double(fhand, &states[st].radius, 1);
            write_int(fhand, &states[st].movable, 1);
            write_int(fhand, &states[st].frozen, 1);
            write_int(fhand, &states[st].index, 1);
            write_int(fhand, &states[st].ixmin, 1);
            write_int(fhand, &states[st].iymin, 1);
            write_int(fhand, &states[st].izmin, 1);
            write_int(fhand, &states[st].ixmax, 1);
            write_int(fhand, &states[st].iymax, 1);
            write_int(fhand, &states[st].izmax, 1);
            write_int(fhand, &states[st].xfold, 1);
            write_int(fhand, &states[st].yfold, 1);
            write_int(fhand, &states[st].zfold, 1);
            write_int(fhand, &states[st].ixstart, 1);
            write_int(fhand, &states[st].iystart, 1);
            write_int(fhand, &states[st].izstart, 1);
            write_int(fhand, &states[st].ixend, 1);
            write_int(fhand, &states[st].iyend, 1);
            write_int(fhand, &states[st].izend, 1);
            write_int(fhand, &states[st].orbit_nx, 1);
            write_int(fhand, &states[st].orbit_ny, 1);
            write_int(fhand, &states[st].orbit_nz, 1);
            write_int(fhand, &states[st].size, 1);
        }
    }


}
