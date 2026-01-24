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
void write_data(char *name, double *vh, double *vxc, double *vh_old,
        double *vxc_old, double *rho, double *vh_corr, STATE * states)
{
    int amode;
    int state;
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
    std::string newname;
    if(pct.spinpe == 0)
    {
        newname = std::format("{}{}", name, ".vh");
        MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
        MPI_File_write_all(mpi_fhand, vh, get_FP0_BASIS(),MPI_DOUBLE, &status);
        MPI_File_close(&mpi_fhand);

        newname = std::format("{}{}", name, ".vh_corr");
        MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
        disp=0;
        MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
        MPI_File_write_all(mpi_fhand, vh_corr, get_FP0_BASIS(),MPI_DOUBLE, &status);
        MPI_File_close(&mpi_fhand);
    }
    MPI_Barrier(pct.grid_comm);

    newname = std::format("{}_spin{}{}", name, pct.spinpe, ".vxc");
    MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(mpi_fhand, vxc, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);

    newname = std::format("{}_spin{}{}", name, pct.spinpe, ".rho");
    MPI_File_open(pct.grid_comm, newname.c_str(), amode, fileinfo, &mpi_fhand);
    disp=0;
    MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(mpi_fhand, rho, get_FP0_BASIS(),MPI_DOUBLE, &status);
    MPI_File_close(&mpi_fhand);
    MPI_Barrier(pct.grid_comm);


    if(pct.gridpe == 0) 
    {
        FILE *fhand_EF;
        newname = std::format("{}{}", name, ".EF");
        fhand_EF = fopen (newname.c_str(), "w");
        fprintf(fhand_EF, "%15.8f\n", ct.efermi * Ha_eV);
        fprintf(fhand_EF, "%15.8f\n", ct.E_lowbound * Ha_eV);
        fclose(fhand_EF);
    }


    if(pct.gridpe == 0) fflush(NULL);


    if(ct.LocalizedOrbitalLayout == LO_projection)
    {
        LocalOrbital->WriteOrbitals(std::string(ct.outfile), *Rmg_G);
    }
    else
    {
        hxgrid = get_hxgrid() * get_xside();
        hygrid = get_hygrid() * get_yside();
        hzgrid = get_hzgrid() * get_zside();

        for (state = ct.state_begin; state < ct.state_end; state++)
        {
            int state_permuted = perm_state_index[state];
            newname = std::format("{}_spin{}{}{}", name, pct.spinpe, ".orbit_", state_permuted);
            amode = S_IREAD | S_IWRITE;
            int fhand = open(newname.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
            if (fhand < 0)
                rmg::error(" Unable to write file ");

            rmg::writefile(fhand, states[state].psiR, states[state].size * sizeof(double));

            rmg::writefile(fhand, &states[state].ixmin, sizeof(int));
            rmg::writefile(fhand, &states[state].ixmax, sizeof(int));
            rmg::writefile(fhand, &states[state].iymin, sizeof(int));
            rmg::writefile(fhand, &states[state].iymax, sizeof(int));
            rmg::writefile(fhand, &states[state].izmin, sizeof(int));
            rmg::writefile(fhand, &states[state].izmax, sizeof(int));
            rmg::writefile(fhand, &hxgrid, sizeof(double));
            rmg::writefile(fhand, &hygrid, sizeof(double));
            rmg::writefile(fhand, &hzgrid, sizeof(double));
            rmg::writefile(fhand, &states[state].crds[0], 3 * sizeof(double));



            close(fhand);
        }
    }
    /* write out the charge density around first atom 
     * mainly for building the initial charge density for FIREBALL start
     */

    if(ct.num_ions > 1) return;

    hxgrid = get_hxxgrid() * get_xside();
    hygrid = get_hyygrid() * get_yside();
    hzgrid = get_hzzgrid() * get_zside();


    ixmin = Rmg_G->default_FG_RATIO * states[0].ixmin;
    ixmax = Rmg_G->default_FG_RATIO * states[0].ixmax;
    iymin = Rmg_G->default_FG_RATIO * states[0].iymin;
    iymax = Rmg_G->default_FG_RATIO * states[0].iymax;
    izmin = Rmg_G->default_FG_RATIO * states[0].izmin;
    izmax = Rmg_G->default_FG_RATIO * states[0].izmax;
    ixdim = ixmax - ixmin;
    iydim = iymax - iymin;
    izdim = izmax - izmin;

    my_malloc_init( rho_tem, ixdim * iydim * izdim, double );
    for(idx = 0; idx < ixdim * iydim * izdim; idx++)
    {
        rho_tem[idx] = 0.0;
    }

    pe2xyz (pct.gridpe, &pex, &pey, &pez);
    PNX0 = pex * get_FPX0_GRID();
    PNY0 = pey * get_FPY0_GRID();
    PNZ0 = pez * get_FPZ0_GRID();
    for (ix = ixmin; ix < ixmax; ix++)
        for (iy = iymin; iy < iymax; iy++)
            for (iz = izmin; iz < izmax; iz++)
            {
                ixx = ix;
                iyy = iy;
                izz = iz;
                if (ixx < 0)
                    ixx += get_FNX_GRID();
                if (iyy < 0)
                    iyy += get_FNY_GRID();
                if (izz < 0)
                    izz += get_FNZ_GRID();

                if (ixx >= get_FNX_GRID())
                    ixx -= get_FNX_GRID();
                if (iyy >= get_FNY_GRID())
                    iyy -= get_FNY_GRID();
                if (izz >= get_FNZ_GRID())
                    izz -= get_FNZ_GRID();

                idx = ixx * get_FNY_GRID() * get_FNZ_GRID() + iyy * get_FNZ_GRID() + izz;
                idx1 = (ix - ixmin) * iydim * izdim + (iy - iymin) * izdim + (iz - izmin);

                ixx -= PNX0;
                iyy -= PNY0;
                izz -= PNZ0;
                idx = ixx *get_FPY0_GRID() * get_FPZ0_GRID() + iyy * get_FPZ0_GRID() + izz;
                if(     ixx >= 0 && ixx < get_FPX0_GRID() &&
                        iyy >= 0 && iyy < get_FPY0_GRID() &&
                        izz >= 0 && izz < get_FPZ0_GRID()  )
                {

                    rho_tem[idx1] = rho[idx];

                }
            }



    idx = ixdim * iydim * izdim ;
    rmg::reduce(rho_tem, idx, pct.grid_comm);
    if (pct.gridpe == 0)
    {
        newname = std::format("{}{}", name, ".rho_firstatom");
        amode = S_IREAD | S_IWRITE;

        int fhand = open(newname.c_str(), O_CREAT | O_TRUNC | O_RDWR, amode);
        if (fhand < 0)
            rmg::error(" Unable to write file ");

        rmg::writefile(fhand, &ixmin, sizeof(int));
        rmg::writefile(fhand, &ixmax, sizeof(int));
        rmg::writefile(fhand, &iymin, sizeof(int));
        rmg::writefile(fhand, &iymax, sizeof(int));
        rmg::writefile(fhand, &izmin, sizeof(int));
        rmg::writefile(fhand, &izmax, sizeof(int));
        rmg::writefile(fhand, &hxgrid, sizeof(double));
        rmg::writefile(fhand, &hygrid, sizeof(double));
        rmg::writefile(fhand, &hzgrid, sizeof(double));
        rmg::writefile(fhand, &states[0].crds[0], 3*sizeof(double));

        rmg::writefile(fhand, rho_tem, idx * sizeof(double));
        close(fhand);
    }



    my_free(rho_tem);


    MPI_Barrier(pct.img_comm);
}


