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
#include "md.h"


void read_rho_and_pot (char *name, double *vh, double *vxc, 
        double *vh_old, double *vxc_old, double *rho)
{
    int fhand;
    char newname[MAX_PATH + 200];
    int idx, pex, pey, pez;

    double tem1, tem2;

    tem2 = my_crtc();

    sprintf(newname, "%s%s", name, ".pot_rho");

    pe2xyz(pct.thispe, &pex, &pey, &pez);

    int sizes[3], subsizes[3], starts[3];
    MPI_Info fileinfo;
    MPI_Datatype  filetype; 
    MPI_Status status;
    MPI_Offset disp, offset;


    /* this datatype describes the mapping of the local array
     * to the global array (file)
     * */

    sizes[0] = FNX_GRID;
    sizes[1] = FNY_GRID;
    sizes[2] = FNZ_GRID;

    subsizes[0] = FPX0_GRID;
    subsizes[1] = FPY0_GRID;
    subsizes[2] = FPZ0_GRID;

    starts[0] = pex * FPX0_GRID;
    starts[1] = pey * FPY0_GRID;
    starts[2] = pez * FPZ0_GRID;

    /*int order = MPI_ORDER_FORTRAN;*/
    int order = MPI_ORDER_C;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, order, MPI_DOUBLE, &filetype);

    MPI_Type_commit(&filetype);

    MPI_Info_create(&fileinfo);


    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;
    MPI_File_open(MPI_COMM_WORLD, newname, amode, fileinfo, &mpi_fhand);


    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

    MPI_File_read(mpi_fhand, vh, FP0_BASIS,MPI_DOUBLE, &status);
    MPI_File_read(mpi_fhand, vxc, FP0_BASIS,MPI_DOUBLE, &status);
    MPI_File_read(mpi_fhand, rho, FP0_BASIS,MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);

    my_barrier();
    tem1 = my_crtc();
    if(pct.thispe == 0) printf("\n aaaa read pot_rho %f", tem2 - tem1);
    fflush(NULL);

}


