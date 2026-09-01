/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#if !(defined(_WIN32) || defined(_WIN64))
    #include <unistd.h>
#else
    #include <io.h>
    #include <BaseTsd.h>
    #define ssize_t SSIZE_T
#endif
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "ZfpCompress.h"

static void read_double (int fhand, double * rp, int count);
static void read_int (int fhand, int *ip, int count);

void read_compressed_buffer(int fh, double *array, int nx, int ny, int nz);

/* Reads the hartree potential, charge density, vxc from RMG results for ON part */
void ReadDataFromRMG (char *name, double * vh, double * rho, double * vxc)
{
    char newname[MAX_PATH + 200];
    int grid[3];
    int fine[3];
    int fpgrid[3];
    int pe[3];
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk;
    int ns;

    fpgrid[0] = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    fpgrid[1] = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    fpgrid[2] = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);

    /* wait until everybody gets here */
    MPI_Barrier(pct.img_comm);	

    /* Make the new output file name */
    rmg_printf("\nspin flag =%d\n", ct.spin_flag);
    
    int kstart = pct.kstart;
    sprintf (newname, "%s_spin%d_kpt%d_gridpe%d", name, pct.spinpe, kstart, pct.gridpe);


    int fhand = open(newname, O_RDWR, S_IREAD | S_IWRITE);
    if (fhand < 0) {
        rmg_printf("Can't open data file %s", newname);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }



    /* read grid info */
    read_int (fhand, grid, 3);
    if (grid[0] != Rmg_G->get_NX_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong NX_GRID");
    if (grid[1] != Rmg_G->get_NY_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong NY_GRID");
    if (grid[2] != Rmg_G->get_NZ_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong NZ_GRID");
    rmg_printf("\n grid %d %d %d\n", grid[0], grid[1], grid[2]);

    /* read grid processor topology */
    read_int (fhand, pe, 3);
    if (pe[0] != Rmg_G->get_PE_X())
        rmg_error_handler (__FILE__, __LINE__,"Wrong PE_X");
    if (pe[1] != Rmg_G->get_PE_Y())
        rmg_error_handler (__FILE__, __LINE__,"Wrong PE_Y");
    if (pe[2] != Rmg_G->get_PE_Z())
        rmg_error_handler (__FILE__, __LINE__,"Wrong PE_Z");

    grid_size = Rmg_G->get_P0_BASIS(0);

    /* read fine grid info */
    read_int (fhand, fine, 3);
    if (fine[0] != Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO) / Rmg_G->get_PX0_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    if (fine[1] != Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO) / Rmg_G->get_PY0_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    if (fine[2] != Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO) / Rmg_G->get_PZ0_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    fgrid_size = grid_size * fine[0] * fine[1] * fine[2];

    /* print out  */
    rmg_printf ("read_data: psi grid = %d %d %d\n", grid[0], grid[1], grid[2]);
    rmg_printf ("read_data: pe grid = %d %d %d\n", pe[0], pe[1], pe[2]);
    rmg_printf ("read_data: grid_size  = %d\n", grid_size);
    rmg_printf ("read_data: fine = %d %d %d\n", fine[0], fine[1], fine[2]);
    rmg_printf ("read_data: fgrid_size = %d\n", fgrid_size);


    /* read wavefunction info */
    read_int (fhand, &gamma, 1);
    //if (gamma != ct.is_gamma)
    //    rmg_error_handler (__FILE__, __LINE__,"Wrong gamma data");


    read_int (fhand, &nk, 1);
    if (nk != ct.num_kpts_pe && ct.forceflag != BAND_STRUCTURE)    /* bandstructure calculation */
        rmg_error_handler (__FILE__, __LINE__,"Wrong number of k points");

    rmg_printf ("read_data: gamma = %d\n", gamma);
    rmg_printf ("read_data: nk = %d\n", ct.num_kpts_pe);

    /* read number of states */  
    read_int (fhand, &ns, 1);


    /* read the hartree potential, electronic density and xc potential */
    if(ct.compressed_infile)
    {
        read_compressed_buffer(fhand, vh, fpgrid[0], fpgrid[1], fpgrid[2]);
        rmg_printf ("read_data: read 'vh'\n");
        read_compressed_buffer(fhand, rho, fpgrid[0], fpgrid[1], fpgrid[2]);
        rmg_printf ("read_data: read 'rho'\n");
        read_compressed_buffer(fhand, vxc, fpgrid[0], fpgrid[1], fpgrid[2]);
        rmg_printf ("read_data: read 'vxc'\n");
    }
    else
    {
        read_double (fhand, vh, fgrid_size);
        rmg_printf ("read_data: read 'vh'\n");
        read_double (fhand, rho, fgrid_size);
        rmg_printf ("read_data: read 'rho'\n");
        read_double (fhand, vxc, fgrid_size);
        rmg_printf ("read_data: read 'vxc'\n");
    }


    close (fhand);


}                               /* end read_data */


static void read_double (int fhand, double * rp, int count)
{
    ssize_t wanted = sizeof (double) * (ssize_t)count;
    ssize_t size = read (fhand, rp, wanted);
    if(size != wanted)
        rmg_error_handler (__FILE__, __LINE__,"error reading");


}

static void read_int (int fhand, int *ip, int count)
{
    int size = count * sizeof (int);
    if (size != read (fhand, ip, size))
        rmg_error_handler (__FILE__, __LINE__,"error reading");
}



void read_compressed_buffer(int fh, double *array, int nx, int ny, int nz)
{

    ZfpCompress C;
    size_t csize;

    double *in = new double[2*nx*ny*nz];

    size_t wsize = read (fh, &csize, sizeof(csize));
    if(wsize != sizeof(csize))
        rmg_error_handler (__FILE__,__LINE__,"error reading");

    if(csize > sizeof(double)*nx*ny*nz)
        rmg_error_handler (__FILE__,__LINE__,"error reading input buffer too small");

    wsize = read (fh, in, csize);
    if(wsize != csize)
        rmg_error_handler (__FILE__,__LINE__,"error reading");

    csize = C.decompress_buffer(array, in, nx, ny, nz, RESTART_TOLERANCE, 2*nx*ny*nz*sizeof(double));
    delete [] in;

}
/******/
/******/
