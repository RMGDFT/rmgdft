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
#include "md.h"





/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
void write_data (char *name, double *vh, double *vxc, double *vh_old, double *vxc_old,
                 double *rho, double *vbias, STATE * states)
{
    int ion;
    int state, i1;
    STATE *sp;
    int idx, st, ntot, i;
    int fhand;
    char newname[MAX_PATH + 20];
    int st1;


    /* Wait until everyone gets here */
    my_barrier ();


    /* Make the new output file name */
    if (pct.gridpe == 0)
    {
        sprintf (newname, "%s%s", name, ".basis");

		fhand = open_wave_file (newname);

		/* Some control information */
		write (fhand, &ct.num_states, sizeof (int));
		write (fhand, &lcr[1].num_states, sizeof (int));
		write (fhand, &lcr[2].num_states, sizeof (int));
		write (fhand, &ct.num_ions, sizeof (int));

		i1 = pct.pe_x;
		write (fhand, &i1, sizeof (int));
		i1 = pct.pe_y;
		write (fhand, &i1, sizeof (int));
		i1 = pct.pe_z;
		write (fhand, &i1, sizeof (int));
		i1 = NPES;
		write (fhand, &i1, sizeof (int));
		i1 = P0_BASIS;
		write (fhand, &i1, sizeof (int));

		/* Write current ionic positions to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write (fhand, &ct.ions[ion].crds[0], 3 * sizeof (double));

		}                       /* end for */

		/* Write original ionic positions to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write (fhand, &ct.ions[ion].icrds[0], 3 * sizeof (double));

		}                       /* end for */

		/* Write current ionic velocities to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write (fhand, &ct.ions[ion].velocity[0], 3 * sizeof (double));

		}


		/* Write current ionic forces pointer array to the file */
		write (fhand, &ct.fpt[0], 4 * sizeof (int));

		/* Write current ionic forces to the file */
		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write (fhand, &ct.ions[ion].force[0][0], 3 * 4 * sizeof (double));

		}                       /* end for */

		for (ion = 0; ion < ct.num_ions; ion++)
		{

			write (fhand, &ct.ions[ion].n_loc_states, sizeof (int));

		}

		write (fhand, &states[0].eig, sizeof (double));
		write (fhand, &ct.TOTAL, sizeof (double));

		/* Occupations */
		sp = states;
		for (state = 0; state < ct.num_states * ct.num_kpts; state++)
		{

			write (fhand, &sp->occupation, sizeof (double));
			sp++;

		}                       /* end for */

		close (fhand);

	}                           /* end if(pct.gridpe == 0) */

	my_barrier ();


	if (pct.gridpe == 0)
	{
		sprintf (newname, "%s%s", name, ".pot_rho");
		my_open (fhand, newname, O_CREAT | O_TRUNC | O_RDWR, S_IREAD | S_IWRITE);
	}

	write_global_data (fhand, vh, FNX_GRID, FNY_GRID, FNZ_GRID);

	write_global_data (fhand, vxc, FNX_GRID, FNY_GRID, FNZ_GRID);

	write_global_data (fhand, rho, FNX_GRID, FNY_GRID, FNZ_GRID);


	idx = ct.num_states * ct.num_states;
	get_distributed_mat (work_matrix, matB);
	global_sums (work_matrix, &idx);
	if (pct.gridpe == 0)
		write (fhand, work_matrix, idx * sizeof (double));

	get_distributed_mat (work_matrix, Hij);
	global_sums (work_matrix, &idx);
	if (pct.gridpe == 0)
		write (fhand, work_matrix, idx * sizeof (double));

	/*        ntot = 0;
	 *        for (i =0; i <ct.num_blocks; i++) ntot += ct.block_dim[i] * ct.block_dim[i];
	 *        for (i =1; i <ct.num_blocks; i++) ntot += ct.block_dim[i-1] * ct.block_dim[i];
	 *  	if(pct.gridpe ==0) write(fhand, lcr[0].Htri, ntot * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[0].Stri, ntot * sizeof(double));
	 *
	 *	idx = lcr[1].num_states * lcr[1].num_states;
	 *  	if(pct.gridpe ==0) write(fhand, lcr[1].H00, idx * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[1].S00, idx * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[1].H01, idx * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[1].S01, idx * sizeof(double));
	 *
	 *	idx = lcr[2].num_states * lcr[2].num_states;
	 *  	if(pct.gridpe ==0) write(fhand, lcr[2].H00, idx * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[2].S00, idx * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[2].H01, idx * sizeof(double));
	 *  	if(pct.gridpe ==0) write(fhand, lcr[2].S01, idx * sizeof(double));
	 */
	my_barrier ();
	/* Force change mode of output file */
	if (pct.gridpe == 0)
		chmod (newname, S_IREAD | S_IWRITE);
	if (pct.gridpe == 0)
		close (fhand);

/* ===================== writing pot ======================= */


    int ix, iy, iz, idx2, FPYZ0, FNXY;
    double *vtot_global, *rho_global;
    FILE *file, *file2;
    int ii, jj, kk;

    for (idx = 0; idx < FP0_BASIS; idx++)
        /*vtot[idx] = vh[idx] + vxc[idx];*/
        vtot[idx] = vh[idx];


    /* apply_potential_drop ( vtot ); */


    FPYZ0 = FPY0_GRID * FPZ0_GRID;

    for (ix = 0; ix < FPX0_GRID; ix++)
    {
        for (iy = 0; iy < FPY0_GRID; iy++)
        {
            idx = iy + ix * FPY0_GRID;
                                                                                                                     
            for (iz = 0; iz < FPZ0_GRID; iz++)
            {
                idx2 = iz + iy * FPZ0_GRID + ix * FPYZ0;
                vtot[idx2] += vbias[idx];
            }
        }
    }


    FNXY = FNX_GRID * FNY_GRID;
    my_malloc_init(vtot_global, FNXY, double);
    my_malloc_init(rho_global, FNXY, double);


    pe2xyz (pct.gridpe, &ii, &jj, &kk);
    ii *= FPX0_GRID;
    jj *= FPY0_GRID;
    kk *= FPZ0_GRID;


    for (iz = 0; iz < FPZ0_GRID; iz++)
    {
        if ((iz + kk) == (FNZ_GRID/6)-1)           /* for a given iz */
        {
            for (ix = 0; ix < FPX0_GRID; ix++)
            {
                for (iy = 0; iy < FPY0_GRID; iy++)
                {
                    idx = ix * FPYZ0 + iy * FPZ0_GRID + iz;
                    idx2 = (ix + ii) * FNY_GRID + (iy + jj);
                    vtot_global[idx2] = vtot[idx];
                    rho_global[idx2] = rho[idx];
                }
            }
        }
    }

    global_sums (vtot_global, &FNXY);
    global_sums (rho_global, &FNXY);


    if (pct.gridpe == 0)
    {          
        file = fopen ("pot.dat", "w");
        file2 = fopen ("rho.dat", "w");

        for (ix = 0; ix < FNX_GRID; ix++)
        {
            for (iy = 0; iy < FNY_GRID; iy++)
            {
                idx = iy + ix * FNY_GRID;

                fprintf (file, " %d %d %f \n", ix, iy, vtot_global[idx]);
                fprintf (file2, " %d %d %f \n", ix, iy, rho_global[idx]);
            }
        }

        fclose (file);
        fclose (file2);
    }

    my_free (vtot_global);
    my_free (rho_global);

/* ================================================= */

	if (ct.runflag < 111)
	{
		for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
		{
			sprintf (newname, "%s%s%d", name, ".orbit_", st1);

			my_open (fhand, newname, O_CREAT | O_TRUNC | O_RDWR, S_IREAD | S_IWRITE);

			idx = states[st1].size * (int) sizeof (double);
			write (fhand, states[st1].psiR, idx);
			idx = sizeof (int);
			write (fhand, &states[st1].ixmin, idx);
			write (fhand, &states[st1].ixmax, idx);
			write (fhand, &states[st1].iymin, idx);
			write (fhand, &states[st1].iymax, idx);
			write (fhand, &states[st1].izmin, idx);
			write (fhand, &states[st1].izmax, idx);

			close (fhand);
		}
	}


}                               /* end write_data */
