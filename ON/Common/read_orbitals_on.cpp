/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


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
#include "transition.h"



/* read orbitals to states.psiR
 */

void read_orbitals_on(char *name, STATE *sts)
{
    int fhand;
    int state;
    unsigned nbytes;
    char newname[MAX_PATH + 200];
    unsigned idx;

    /* Wait until everybody gets here */
    MPI_Barrier(pct.img_comm);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        sprintf(newname, "%s_spin%d%s%d", name, pct.spinpe, ".res_", state);
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
        {
            printf("\n  unable to open file %s", newname);
            exit(0);
        }

        nbytes = read(fhand, sts[state].psiR, sts[state].size * sizeof(double));
        idx = sts[state].size * sizeof(double);
        if (nbytes != idx)
        {
            rmg_printf("\n read %d is different from %d for state %d", nbytes, idx, state);
            rmg_error_handler(__FILE__, __LINE__, "Unexpected end of file orbit");
        }


        close(fhand);
    }

    MPI_Barrier(pct.img_comm);


}                               /* end read_data */
