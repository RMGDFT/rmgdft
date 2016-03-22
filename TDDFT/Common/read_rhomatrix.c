/************************** SVN Revision Information **************************
 **    $Id: write_data.c 3152 2015-08-14 22:25:39Z luw $    **
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


/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void read_rhomatrix(char *name, double *rho_matrix_row)
{
    int amode;
    char newname[MAX_PATH + 20];
    int idx, idx1;
    int fhand;
    /* Wait until everyone gets here */
    my_barrier();

    /* Make the new output file name */

	int sizes[3], subsizes[3], starts[3];
	MPI_Info fileinfo;
	MPI_Datatype  filetype; 
	MPI_Status status;
	MPI_Offset disp;



	MPI_File mpi_fhand ;

	sprintf(newname, "%s%s", name, ".rho_matrix");

	/* this datatype describes the mapping of the local array
	 * to the global array (file)
	 * */

	sizes[0] = ct.num_states;
	sizes[1] = ct.num_states;

	subsizes[0] = ct.state_end - ct.state_begin;
	subsizes[1] = ct.num_states;

	starts[0] = pct.gridpe * (ct.state_end - ct.state_begin);
	starts[1] = 0;

    int order = MPI_ORDER_C;

	/*int order = MPI_ORDER_FORTRAN;*/
	MPI_Type_create_subarray(2, sizes, subsizes, starts, order, MPI_DOUBLE, &filetype);

	MPI_Type_commit(&filetype);

	MPI_Info_create(&fileinfo);


	amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
	MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);


	disp=0;
	MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

    
    idx = (ct.state_end - ct.state_begin) * ct.num_states;
	MPI_File_read(mpi_fhand, rho_matrix_row, idx,MPI_DOUBLE, &status);
	MPI_File_close(&mpi_fhand);

	my_barrier();
}
