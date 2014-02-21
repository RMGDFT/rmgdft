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
#include "main_on.h"
#include <math.h>

#include <mpi.h>


/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void write_data(char *name, double *vh, double *vxc, double *vh_old,
        double *vxc_old, double *rho, STATE * states)
{
    int amode;
    int ion;
    int state, i1;
    STATE *sp;
    rmg_double_t *work;
    char newname[MAX_PATH + 20];
    rmg_double_t tem, tem1, tem2, tem3;
    int idx, idx1, st;
    int fhand;
    int ione = 1, st1, st2;
    rmg_double_t dis;
    double hxgrid, hygrid, hzgrid;
    double *rho_tem;
    int ix, iy, iz, ixdim, iydim, izdim;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    int ixx, iyy, izz;
    int ipe, NX, NY, NZ, PNX0, PNY0, PNZ0;
    int ix1, iy1, pex, pey, pez, position;

    time_t tt;

    char *timeptr;



    /* Wait until everyone gets here */
    my_barrier();

    /* Make the new output file name */
    if (pct.gridpe == 0)
    {
		sprintf(newname, "%s%s", name, ".basis");

		fhand = open_wave_file (newname);

		if (fhand < 0)
			error_handler(" Unable to write file ");


		/* Some control information */
		write(fhand, &ct.num_states, sizeof(int));
		write(fhand, &ct.num_ions, sizeof(int));

		i1 = pct.pe_x;
		write(fhand, &i1, sizeof(int));
		i1 = pct.pe_y;
		write(fhand, &i1, sizeof(int));
		i1 = pct.pe_z;
		write(fhand, &i1, sizeof(int));
		i1 = NPES;
		write(fhand, &i1, sizeof(int));
		i1 = get_P0_BASIS();
		write(fhand, &i1, sizeof(int));


		/* Write current ionic positions to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write(fhand, &ct.ions[ion].crds[0], 3 * sizeof(double));

		}                       /* end for */


		/* Write original ionic positions to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write(fhand, &ct.ions[ion].icrds[0], 3 * sizeof(double));

		}                       /* end for */


		/* Write current ionic velocities to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write(fhand, &ct.ions[ion].velocity[0], 3 * sizeof(double));

		}                       /* end for */



		/* Write current ionic forces pointer array to the file */
		write(fhand, &ct.fpt[0], 4 * sizeof(int));

		/* Write current ionic forces to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write(fhand, &ct.ions[ion].force[0][0], 3 * 4 * sizeof(double));

		}                       /* end for */


		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write(fhand, &ct.ions[ion].n_loc_states, sizeof(int));

		}

		write(fhand, &states[0].eig, sizeof(double));
		write(fhand, &ct.TOTAL, sizeof(double));




		/* Occupations */
		sp = states;
		for (state = 0; state < ct.num_states * ct.num_kpts; state++)
		{

			write(fhand, &sp->occupation, sizeof(double));
			sp++;

		}                       /* end for */

		close(fhand);

	}                           /* end if(pct.gridpe == 0) */

	my_barrier();

	time(&tt);
	timeptr = ctime(&tt);

	if(pct.gridpe == 0) printf("\n    Write start at %s\n", timeptr);
	if(pct.gridpe == 0) fflush(NULL);

	sprintf(newname, "%s%s", name, ".pot_rho");

	pe2xyz (pct.gridpe, &pex, &pey, &pez);

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


	amode = MPI_MODE_RDWR|MPI_MODE_CREATE;
	MPI_File mpi_fhand ;
	MPI_File_open(MPI_COMM_WORLD, newname, amode, fileinfo, &mpi_fhand);


	disp=0;
	MPI_File_set_view(mpi_fhand, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

	MPI_File_write_all(mpi_fhand, vh, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_write_all(mpi_fhand, vxc, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_write_all(mpi_fhand, rho, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_write_all(mpi_fhand, vh_old, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_write_all(mpi_fhand, vxc_old, get_FP0_BASIS(),MPI_DOUBLE, &status);
	MPI_File_close(&mpi_fhand);

	my_barrier();

	time(&tt);
	timeptr = ctime(&tt);

	if(pct.gridpe == 0) printf("\n    Write middle at %s\n", timeptr);
	if(pct.gridpe == 0) fflush(NULL);



	/* Force change mode of output file */
	amode = S_IREAD | S_IWRITE;
	if (pct.gridpe == 0)
		chmod(newname, amode);
	if (pct.gridpe == 0)
		close(fhand);

	hxgrid = ct.hxgrid * ct.xside;
	hygrid = ct.hygrid * ct.yside;
	hzgrid = ct.hzgrid * ct.zside;

	for (state = ct.state_begin; state < ct.state_end; state++)
	{
		sprintf(newname, "%s%s%d", name, ".orbit_", state);
		amode = S_IREAD | S_IWRITE;
		fhand = open(newname, O_CREAT | O_TRUNC | O_RDWR, amode);
		if (fhand < 0)
			error_handler(" Unable to write file ");

		write(fhand, states[state].psiR, states[state].size * sizeof(double));

		write(fhand, &states[state].ixmin, sizeof(int));
		write(fhand, &states[state].ixmax, sizeof(int));
		write(fhand, &states[state].iymin, sizeof(int));
		write(fhand, &states[state].iymax, sizeof(int));
		write(fhand, &states[state].izmin, sizeof(int));
		write(fhand, &states[state].izmax, sizeof(int));
		write(fhand, &hxgrid, sizeof(double));
		write(fhand, &hygrid, sizeof(double));
		write(fhand, &hzgrid, sizeof(double));
		write(fhand, &states[state].crds[0], 3 * sizeof(double));



		close(fhand);
	}

	/* write out the charge density around first atom 
	 * mainly for building the initial charge density for FIREBALL start
	 */


	hxgrid = ct.hxxgrid * ct.xside;
	hygrid = ct.hyygrid * ct.yside;
	hzgrid = ct.hzzgrid * ct.zside;


	ixmin = 2 * states[0].ixmin;
	ixmax = 2 * states[0].ixmax;
	iymin = 2 * states[0].iymin;
	iymax = 2 * states[0].iymax;
	izmin = 2 * states[0].izmin;
	izmax = 2 * states[0].izmax;
	ixdim = ixmax - ixmin;
	iydim = iymax - iymin;
	izdim = izmax - izmin;

	my_malloc_init( rho_tem, ixdim * iydim * izdim, rmg_double_t );
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
	global_sums(rho_tem, &idx, pct.grid_comm);
	idx = ixdim * iydim * izdim * sizeof(double);
	if (pct.gridpe == 0)
	{
		sprintf(newname, "%s%s", name, ".rho_firstatom");
		amode = S_IREAD | S_IWRITE;

		fhand = open(newname, O_CREAT | O_TRUNC | O_RDWR, amode);
		if (fhand < 0)
			error_handler(" Unable to write file ");

		write(fhand, &ixmin, sizeof(int));
		write(fhand, &ixmax, sizeof(int));
		write(fhand, &iymin, sizeof(int));
		write(fhand, &iymax, sizeof(int));
		write(fhand, &izmin, sizeof(int));
		write(fhand, &izmax, sizeof(int));
		write(fhand, &hxgrid, sizeof(double));
		write(fhand, &hygrid, sizeof(double));
		write(fhand, &hzgrid, sizeof(double));
		write(fhand, &states[0].crds[0], 3 * sizeof(double));

		write(fhand, rho_tem, idx);
		close(fhand);
	}



	my_free(rho_tem);


	my_barrier();


}                               /* end write_data */
