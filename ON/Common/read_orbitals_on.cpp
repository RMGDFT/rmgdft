/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   rdwrpg.c


   Functions to read data from files.


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
#include "init_var.h"



/* read orbitals to states.psiR
 */

void read_orbitals_on(char *name, STATE *sts)
{
    int fhand;
    int state;
    unsigned nbytes;
    char newname[MAX_PATH + 200];
    int idx;
    int pex, pey, pez;

    /* Wait until everybody gets here */
    MPI_Barrier(pct.img_comm);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        sprintf(newname, "%s%s%d", name, ".orbit_", state);
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
        {
            dprintf("\n  unable to open file %s", newname);
            exit(0);
        }

        nbytes = read(fhand, sts[state].psiR, sts[state].size * sizeof(double));
        idx = sts[state].size * sizeof(double);
        if (nbytes != idx)
        {
            printf("\n read %d is different from %d for state %d", nbytes, idx, state);
            error_handler("Unexpected end of file orbit");
        }


        close(fhand);
    }

    MPI_Barrier(pct.img_comm);


}                               /* end read_data */
