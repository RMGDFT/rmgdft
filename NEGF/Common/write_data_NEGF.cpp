/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   write_data.c


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
#include <time.h>
#include "main.h"
#include "prototypes_on.h"
#include <math.h>
#include "init_var.h"
#include <mpi.h>
#include "transition.h"


/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void write_data_NEGF(char *name, double *vh, double *vxc, double *rho)
{
    int amode;
    int state;
    char newname[MAX_PATH + 20];
    int idx, idx1;
    double hxgrid, hygrid, hzgrid;
    double *rho_tem;
    int ix, iy, iz, ixdim, iydim, izdim;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int ixx, iyy, izz;
    int PNX0, PNY0, PNZ0;
    int pex, pey, pez;


    /* Wait until everyone gets here */
    MPI_Barrier(pct.img_comm);

    /* Make the new output file name */

	if(pct.gridpe == 0) fflush(NULL);

     if (pct.imgpe == 0)
    {
        FILE *fhandle;
        fhandle = open_restart_file (name);
        fclose(fhandle);
    }


	pe2xyz (pct.gridpe, &pex, &pey, &pez);

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

	starts[0] = pex * get_FPX0_GRID();
	starts[1] = pey * get_FPY0_GRID();
	starts[2] = pez * get_FPZ0_GRID();

	/*int order = MPI_ORDER_FORTRAN;*/
	int order = MPI_ORDER_C;
	MPI_Type_create_subarray(3, sizes, subsizes, starts, order, MPI_DOUBLE, &filetype);

	MPI_Type_commit(&filetype);

	MPI_Info_create(&fileinfo);


	amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
	MPI_File mpi_fhand ;

    MPI_Barrier(pct.grid_comm);
	sprintf(newname, "%s%s", name, ".vh");
	MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
	disp=0;
	MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpi_fhand, vh, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

	sprintf(newname, "%s_spin%d%s", name, pct.spinpe, ".vxc");
	MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
	disp=0;
	MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpi_fhand, vxc, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

	sprintf(newname, "%s_spin%d%s", name, pct.spinpe, ".rho");
	MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
	disp=0;
	MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpi_fhand, rho, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

	MPI_Barrier(pct.img_comm);


}
