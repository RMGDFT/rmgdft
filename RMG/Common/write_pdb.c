/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include "main.h"



void write_pdb (void)
{
    int ion;
    ION *iptr;



    //printf("\n\n Updated PDB file");

    for (ion = 0; ion < ct.num_ions; ion++)
    {

        iptr = &ct.ions[ion];


        printf ("\n");
        printf ("%-6s", iptr->pdb.record_name);
        printf ("%5d ", iptr->pdb.serial_num);
        printf ("%4s", iptr->pdb.name);
        printf ("%1s", iptr->pdb.altLoc);
        printf ("%3s ", iptr->pdb.resName);
        printf ("%1s", iptr->pdb.chainID);
        printf ("%4d", iptr->pdb.resSeq);
#if 1
        printf ("%1s   ", iptr->pdb.iCode);
        printf ("%8.3f", a0_A * iptr->crds[0]);
        printf ("%8.3f", a0_A * iptr->crds[1]);
        printf ("%8.3f", a0_A * iptr->crds[2]);
        //printf("%6s", iptr->pdb_occupancy);
        printf ("%6.2f", iptr->pdb.occupancy);
        //printf("%6s      ", iptr->pdb_tempFactor);
        printf ("%6.2f          ", iptr->pdb.tempFactor);
        //printf("%4s", iptr->pdb_segID);
        printf ("%2s", ct.sp[iptr->species].pseudo_symbol);
        printf ("%2s", iptr->pdb.charge);
#endif

    }                           /*end for(ion = 0;ion < c->num_ions;ion++) */



    printf ("\nEND");
    printf ("\n\n");


}                               /* EOF */
