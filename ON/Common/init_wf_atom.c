/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*
     	Just generates a random start.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "main.h"
#include "prototypes_on.h"
#include "init_var.h"

void init_wf_atom(STATE * states)
{

    int idx, state, ix, iy, iz;
    STATE *sp;
    long idum;
    int ix1, iy1, iz1;
    int ixx, iyy, izz;
    rmg_double_t temp;
    char newname[MAX_PATH + 200];
    int i, ii;
    int ion, species, ist, fhand, nbytes;


    if (pct.gridpe == 0)
        printf(" readin initial wavefunction \n");
    my_barrier();

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        ion = state_to_ion[state];
        species = ct.ions[ion].species;
        ist = states[state].atomic_orbital_index;



        sprintf(newname, "%s%s%d", ct.file_atomic_orbit[species], ".orbit_", ist);
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
        {
            printf("\n unable to open file: %s \n", newname);
            error_handler(" Unable to open file ");
        }

        idx = states[state].size * sizeof(double);
        nbytes = read(fhand, states[state].psiR, idx);
        if (nbytes != idx)
        {
            printf("\n read %d is different from %d ", nbytes, idx);
            printf("\n file name: %s\n", newname);

            error_handler("Unexpected end of file orbit");
        }


    }



    if (pct.gridpe == 0)
        printf(" readin initial orbitals  done  \n");


}                               /* end init_wf_atom */
