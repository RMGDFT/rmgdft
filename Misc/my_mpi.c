/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include "main.h"


void my_barrier ()
{
    MPI_Barrier (pct.img_comm);
    /*MPI_Barrier (MPI_COMM_WORLD);*/

}

