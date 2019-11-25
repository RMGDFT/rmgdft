#include "negf_prototypes.h"
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
#include "init_var.h"
#include "LCR.h"


void read_rho_and_pot (char *name, double *vh, double *vxc, 
        double *vh_old, double *vxc_old, double *rho)
{
    char newname[MAX_PATH + 200];

    double tem1, tem2;

    tem2 = my_crtc();



    int sizes[3], subsizes[3], starts[3];
    MPI_Info fileinfo;
    MPI_Datatype  filetype; 
    MPI_Status status;
    MPI_Offset disp;


    /* this datatype describes the mapping of the local array
     * to the global array (file)
     * */

    sizes[0] = get_FNX_GRID();
    sizes[1] = get_FNY_GRID();
    sizes[2] = get_FNZ_GRID();

    subsizes[0] = get_FPX0_GRID();
    subsizes[1] = get_FPY0_GRID();
    subsizes[2] = get_FPZ0_GRID();

    starts[0] = get_FPX_OFFSET();
    starts[1] = get_FPY_OFFSET();
    starts[2] = get_FPZ_OFFSET();

    /*int order = MPI_ORDER_FORTRAN;*/
    int order = MPI_ORDER_C;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, order, MPI_DOUBLE, &filetype);

    MPI_Type_commit(&filetype);

    MPI_Info_create(&fileinfo);


    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    sprintf(newname, "%s%s", name, ".vh");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read(mpi_fhand, vh, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);

    sprintf(newname, "%s%s", name, ".vxc");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read(mpi_fhand, vxc, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);

    sprintf(newname, "%s%s", name, ".rho");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read(mpi_fhand, rho, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);

    MPI_Barrier(pct.img_comm);
    tem1 = my_crtc();
    if(pct.gridpe == 0) printf("\n time for read vh, vxc, rho %f", tem2 - tem1);
    fflush(NULL);

}


