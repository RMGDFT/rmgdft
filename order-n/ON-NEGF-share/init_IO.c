/************************** SVN Revision Information **************************
 **    $Id: init_pe.c 1176 2011-01-11 21:34:44Z luw $    **
******************************************************************************/
 

#include "main.h"
#include <stdlib.h>
#include <stdio.h>




void init_IO(void)
{

    MPI_Comm_dup ( MPI_COMM_WORLD, &pct.thisgrp_comm);
    MPI_Comm_dup ( MPI_COMM_WORLD, &pct.img_comm);
    MPI_Comm_dup ( MPI_COMM_WORLD, &pct.grid_comm);
    ct.logfile = stdout;
    pct.instances = 1;
    pct.imgpe = pct.thispe;
    pct.images = pct.thispe;


}
