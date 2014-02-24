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



/* Reads the potential, the wavefunctions and various other things
   from a file. 
 */

void read_data(char *name, double *vh, double *vxc, double *vh_old,
        double *vxc_old, double *rho, STATE * states)
{
    int fhand, ion, ii, jj, ione = 1, ithree = 3;
    int state, st, t1;
    int i, ic, flag, temp_int[MAX_LOC_ST];
    off_t newoffset;
    STATE *sp;
    double temp_occ;
    double temp_pos[3];
    unsigned nbytes;
    int ipe;
    int numst;
    char newname[MAX_PATH + 200];
    rmg_double_t tem, tem1, tem2, tem3;
    int idx;
    int idx1, ix, iy, PNX0, PNY0, PNZ0, NX, NY, NZ;
    int pex, pey, pez, position, ix1, iy1;

    numst = get_P0_BASIS();
    /* Wait until everybody gets here */
    my_barrier();

tem1 = my_crtc();

    sprintf(newname, "%s%s", name, ".basis");
    fhand = open(newname, O_RDWR);
    if (fhand < 0)
        error_handler(" Unable to open file ");


    /* Some control information */

    nbytes = read(fhand, &t1, sizeof(int));
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");
    if (pct.gridpe == 0)
        printf(" %d states\n", t1);
    if (t1 != ct.num_states)
        error_handler("Wrong number of states");

    nbytes = read(fhand, &t1, sizeof(int));
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");
    if (pct.gridpe == 0)
        printf(" Read %d ions\n", t1);
    if (t1 != ct.num_ions)
        error_handler("Wrong number of ions");


    nbytes = read(fhand, &t1, sizeof(int));
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");


    nbytes = read(fhand, &t1, sizeof(int));
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");


    nbytes = read(fhand, &t1, sizeof(int));
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");

    nbytes = read(fhand, &t1, sizeof(int));
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");


    nbytes = read(fhand, &t1, sizeof(int));
    if (pct.gridpe == 0)
        printf(" Read functions on %d grid points\n", t1);
    if (nbytes != sizeof(int))
        error_handler("Unexpected end of file");


    /* Read current ionic positions from the file */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        nbytes = read(fhand, &temp_pos[0], 3 * sizeof(double));
        if (nbytes != 3 * sizeof(double))
            error_handler("Unexpected end of file");

        if (ct.override_atoms != 1)
        {
            for (i = 0; i < 3; i++)
            {
                ct.ions[ion].crds[i] = temp_pos[i];
            }                   /* i */
        }                       /* end if "override" */

        to_crystal(ct.ions[ion].xtal, ct.ions[ion].crds);

    }                           /* end for */
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        if (pct.gridpe == 0)
            printf("\n  %d   %10.4f  %10.4f  %10.4f",
                    ct.ions[ion].species + 1,
                    ct.ions[ion].crds[0], ct.ions[ion].crds[1], ct.ions[ion].crds[2]);

    }                           /* end for */


    /* If constrained dynamics is on and you want to overwrite
       the initial positions */



    /* Read initial ionic positions from the file */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        nbytes = read(fhand, &ct.ions[ion].icrds[0], 3 * sizeof(double));
        if (nbytes != 3 * sizeof(double))
            error_handler("Unexpected end of file");

    }
    if (pct.gridpe == 0)
        printf("\n Initial coordinates \n");
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        if (pct.gridpe == 0)
            printf("  %d   %10.4f  %10.4f  %10.4f \n",
                    ct.ions[ion].species + 1,
                    ct.ions[ion].icrds[0], ct.ions[ion].icrds[1], ct.ions[ion].icrds[2]);

    }                           /* end for */



    /* Read ionic velocities from the file */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        nbytes = read(fhand, &ct.ions[ion].velocity[0], 3 * sizeof(double));
        if (nbytes != 3 * sizeof(double))
            error_handler("Unexpected end of file");

    }                           /* end for */

    /* Read forces pointer from the file */
    read(fhand, &ct.fpt[0], 4 * sizeof(int));

    /* Read ionic forces from the file */
    flag = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        nbytes = read(fhand, &ct.ions[ion].force[0][0], 3 * 4 * sizeof(double));

    }


    /* Read current ionic localizations from the file */


    for (ion = 0; ion < ct.num_ions; ion++)
    {

        read(fhand, &ct.ions[ion].n_loc_states, sizeof(int));

    }

    read(fhand, &states[0].eig, sizeof(double));
    read(fhand, &ct.TOTAL, sizeof(double));


    /* Occupations */
    sp = states;
    for (state = 0; state < ct.num_states * ct.num_kpts; state++)
    {

        nbytes = read(fhand, &temp_occ, sizeof(double));
        if (nbytes != sizeof(double))
            error_handler("Unexpected end of file");

        sp->occupation[0] = temp_occ;

        sp++;

    }                           /* end for */


    close(fhand);

    tem2 = my_crtc();
    if(pct.gridpe == 0) printf("\n aaaa read basis %f", tem2 - tem1);
    fflush(NULL);

    sprintf(newname, "%s%s", name, ".pot_rho");

    pe2xyz(pct.gridpe, &pex, &pey, &pez);

    int sizes[3], subsizes[3], starts[3];
    MPI_Info fileinfo;
    MPI_Datatype  filetype; 
    MPI_Status status;
    MPI_Offset disp, offset;


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
    MPI_File_open(MPI_COMM_WORLD, newname, amode, fileinfo, &mpi_fhand);


    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

    MPI_File_read(mpi_fhand, vh, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_read(mpi_fhand, vxc, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_read(mpi_fhand, rho, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_read(mpi_fhand, vh_old, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_read(mpi_fhand, vxc_old, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);

    my_barrier();
    tem1 = my_crtc();
    if(pct.gridpe == 0) printf("\n aaaa read pot_rho %f", tem2 - tem1);
    fflush(NULL);

    for (state = ct.state_begin; state < ct.state_end; state++)
    {
        sprintf(newname, "%s%s%d", name, ".orbit_", state);
        fhand = open(newname, O_RDWR);
        if (fhand < 0)
            error_handler(" Unable to write file ");

        nbytes = read(fhand, states[state].psiR, states[state].size * sizeof(double));
        idx = states[state].size * sizeof(double);
        if (nbytes != idx)
        {
            printf("\n read %d is different from %d for state %d", nbytes, idx, state);
            error_handler("Unexpected end of file orbit");
        }


        close(fhand);
    }

    my_barrier();

    tem1 = my_crtc();
    if(pct.gridpe == 0) printf("\n aaaa read orbit %f", tem2 - tem1);
    fflush(NULL);
    fflush(NULL);


}                               /* end read_data */
