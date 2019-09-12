/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   rdwrpg.c


   Functions to read data from files.


 */



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
#include "LdaU_on.h"



/* Reads the potential, the wavefunctions and various other things
   from a file. 
 */

void read_data(char *name, double *vh, double *vxc, double *vh_old,
        double *vxc_old, double *rho, double *vh_corr, STATE * states)
{
    int fhand;
    int state;
    size_t nbytes;
    char newname[MAX_PATH + 200];
    int idx;
    int pex, pey, pez;

    /* Wait until everybody gets here */
    MPI_Barrier(pct.img_comm);


    pe2xyz(pct.gridpe, &pex, &pey, &pez);

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


    int amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
    MPI_File mpi_fhand ;

    sprintf(newname, "%s%s", name, ".vh");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(mpi_fhand, vh, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

    sprintf(newname, "%s_spin%d%s", name, pct.spinpe, ".vxc");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(mpi_fhand, vxc, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

    sprintf(newname, "%s_spin%d%s", name, pct.spinpe, ".rho");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(mpi_fhand, rho, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

    sprintf(newname, "%s%s", name, ".vh_corr");
    MPI_File_open(pct.grid_comm, newname, amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(mpi_fhand, vh_corr, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

    MPI_Barrier(pct.img_comm);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        int state_permuted = perm_state_index[state];
        sprintf(newname, "%s_spin%d%s%d", name, pct.spinpe, ".orbit_", state_permuted);
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
        {
            dprintf("\n  unable to open file %s", newname);
            exit(0);
        }

        nbytes = read(fhand, states[state].psiR, states[state].size * sizeof(double));
        idx = states[state].size * sizeof(double);
        if (nbytes != (size_t)idx)
        {
            printf("\n read %d is different from %d for state %d", nbytes, idx, state);
            error_handler("Unexpected end of file orbit");
        }


//      read(fhand, &states[state].ixmin, sizeof(int));
//      read(fhand, &states[state].ixmax, sizeof(int));
//      read(fhand, &states[state].iymin, sizeof(int));
//      read(fhand, &states[state].iymax, sizeof(int));
//      read(fhand, &states[state].izmin, sizeof(int));
//      read(fhand, &states[state].izmax, sizeof(int));
//      double hxgrid, hygrid, hzgrid;
//      read(fhand, &hxgrid, sizeof(double));
//      read(fhand, &hygrid, sizeof(double));
//      read(fhand, &hzgrid, sizeof(double));
//      read(fhand, &states[state].crds[0], 3 * sizeof(double));


        close(fhand);
    }


    if(ct.num_ldaU_ions > 0)
    {
        ldaU_on->ReadLdaU(std::string(name), *LocalOrbital);
    }


    MPI_Barrier(pct.img_comm);

    fflush(NULL);


}                               /* end read_data */
