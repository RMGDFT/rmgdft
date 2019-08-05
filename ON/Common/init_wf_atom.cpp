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

    int idx, state;
    char newname[MAX_PATH + 200];
    int ion, species, ist, fhand, nbytes;

    
    RmgTimer *RT1 = new RmgTimer("1-TOTAL: init: aromic orbita");

    if (pct.gridpe == 0)
        printf(" readin initial wavefunction \n");
    MPI_Barrier(pct.img_comm);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        ion = states[state].atom_index;
        species = Atoms[ion].species;
        ist = states[state].atomic_orbital_index;



        sprintf(newname, "%s%s%s%d", pct.image_path[pct.thisimg], ct.file_atomic_orbit[species].c_str(), "_spin0.orbit_", ist);
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
        {
            printf("\n ddd %d %d %d", state, ion, species);
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
    delete RT1;


}                               /* end init_wf_atom */


