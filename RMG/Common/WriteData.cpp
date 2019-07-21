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
#endif
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"
#include "ZfpCompress.h"

static size_t totalsize;


static void write_double (int fh, double * rp, int count);
static void write_int (int fh, int *ip, int count);


template void WriteData (int, double *, double *, double *, Kpoint<double> **);
template void WriteData (int, double *, double *, double *, Kpoint<std::complex<double> > **);

void write_compressed_buffer(int fh, double *array, int nx, int ny, int nz);

/* Writes the hartree potential, the wavefunctions, the */
/* compensating charges and various other things to a file. */
template <typename KpointType>
void WriteData (int fhand, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr)
{
    int fine[3];
    int grid[3];
    int pgrid[3];
    int fpgrid[3];
    int pe[3];
    int npe;
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;
    double time0, write_time;

    time0 = my_crtc ();
    totalsize = 0;

    pgrid[0] = Kptr[0]->G->get_PX0_GRID(1);
    pgrid[1] = Kptr[0]->G->get_PY0_GRID(1);
    pgrid[2] = Kptr[0]->G->get_PZ0_GRID(1);
    fpgrid[0] = Kptr[0]->G->get_PX0_GRID(Kptr[0]->G->default_FG_RATIO);
    fpgrid[1] = Kptr[0]->G->get_PY0_GRID(Kptr[0]->G->default_FG_RATIO);
    fpgrid[2] = Kptr[0]->G->get_PZ0_GRID(Kptr[0]->G->default_FG_RATIO);

    /* write grid info */
    grid[0] = Kptr[0]->G->get_NX_GRID(1);
    grid[1] = Kptr[0]->G->get_NY_GRID(1);
    grid[2] = Kptr[0]->G->get_NZ_GRID(1);
    write_int (fhand, grid, 3);

    /* write grid processor topology */
    pe[0] = Kptr[0]->G->get_PE_X();
    pe[1] = Kptr[0]->G->get_PE_Y();
    pe[2] = Kptr[0]->G->get_PE_Z();
    write_int (fhand, pe, 3);

    npe = (pe[0] * pe[1] * pe[2]);
    grid_size = Kptr[0]->pbasis;

    /* write fine grid info */
    fine[0] = Kptr[0]->G->get_PX0_GRID(Kptr[0]->G->default_FG_RATIO) / Kptr[0]->G->get_PX0_GRID(1);
    fine[1] = Kptr[0]->G->get_PY0_GRID(Kptr[0]->G->default_FG_RATIO) / Kptr[0]->G->get_PY0_GRID(1);
    fine[2] = Kptr[0]->G->get_PZ0_GRID(Kptr[0]->G->default_FG_RATIO) / Kptr[0]->G->get_PZ0_GRID(1);
    write_int (fhand, fine, 3);
    fgrid_size = grid_size * fine[0] * fine[1] * fine[2];


    /* write wavefunction info */
    gamma = ct.is_gamma;
    nk = ct.num_kpts_pe;
    write_int (fhand, &gamma, 1);
    write_int (fhand, &nk, 1); 

    ns = ct.num_states;
    write_int (fhand, &ns, 1); 




    /* write the hartree potential, electronic density and xc potentials */
    if(ct.compressed_outfile)
    {
        write_compressed_buffer(fhand, vh, fpgrid[0], fpgrid[1], fpgrid[2]);
        write_compressed_buffer(fhand, rho, fpgrid[0], fpgrid[1], fpgrid[2]);
        write_compressed_buffer(fhand, vxc, fpgrid[0], fpgrid[1], fpgrid[2]);
    }
    else
    {
        write_double (fhand, vh, fgrid_size);
        write_double (fhand, rho, fgrid_size);
        write_double (fhand, vxc, fgrid_size);
    }

    /* write wavefunctions */
    {
        double *psi_R = new double[grid_size];
        double *psi_I = new double[grid_size];

        int wvfn_size = (gamma) ? grid_size : 2 * grid_size;

        for (ik = 0; ik < ct.num_kpts_pe; ik++)
        {
            for (is = 0; is < ns; is++)
            {
                if(ct.compressed_outfile)
                {
                    if(gamma)
                    {
                        write_compressed_buffer(fhand, (double *)Kptr[ik]->Kstates[is].psi, pgrid[0], pgrid[1], pgrid[2]);
                    }
                    else
                    {
                        for(int idx=0;idx < grid_size;idx++) psi_R[idx] = std::real(Kptr[ik]->Kstates[is].psi[idx]);
                        for(int idx=0;idx < grid_size;idx++) psi_I[idx] = std::imag(Kptr[ik]->Kstates[is].psi[idx]);
                        write_compressed_buffer(fhand, psi_R, pgrid[0], pgrid[1], pgrid[2]);
                        write_compressed_buffer(fhand, psi_I, pgrid[0], pgrid[1], pgrid[2]);
                    }
                }
                else
                {
                    write_double (fhand, (double *)Kptr[ik]->Kstates[is].psi, wvfn_size);
                }
            }
        }

        delete [] psi_I;
        delete [] psi_R;
    }


    /* write the state occupations, in spin-polarized calculation, 
     * it's occupation for the processor's own spin */ 
    {
	for (ik = 0; ik < ct.num_kpts_pe; ik++)
	    for (is = 0; is < ns; is++)
	    {
		write_double (fhand, &Kptr[ik]->Kstates[is].occupation[0], 1); 
	    }
    }
    

    /* write the state eigenvalues, while in spin-polarized case, 
     * it's eigenvalues of processor's own spin */
    {
	for (ik = 0; ik < ct.num_kpts_pe; ik++)
	    for (is = 0; is < ns; is++)
	    {
		write_double (fhand, &Kptr[ik]->Kstates[is].eig[0], 1);
	    }

    }


    write_time = my_crtc () - time0;
    
    rmg_printf ("WriteData: total size of each of the %d files = %.1f Mb\n", npe,
	    ((double) totalsize) / (1024 * 1024));
    rmg_printf ("WriteData: writing took %.1f seconds, writing speed %.3f Mbps \n", write_time,
	    ((double) totalsize) / (1024 * 1024) / write_time);




}                               /* end write_data */



static void write_double (int fh, double * rp, int count)
{
    int size;

    size = count * sizeof (double);
    if (size != write (fh, rp, size))
        rmg_error_handler (__FILE__,__LINE__,"error writing");

    totalsize += size;
}


static void write_int (int fh, int *ip, int count)
{
    int size;

    size = count * sizeof (int);
    if (size != write (fh, ip, size))
        rmg_error_handler (__FILE__, __LINE__, "error writing");

    totalsize += size;
}

void write_compressed_buffer(int fh, double *array, int nx, int ny, int nz)
{

    ZfpCompress C;
    double *out = new double[2*nx*ny*nz];

    size_t csize = C.compress_buffer(array, out, nx, ny, nz, RESTART_TOLERANCE, 2*nx*ny*nz*sizeof(double));
    size_t wsize = write (fh, &csize, sizeof(csize));
    if(wsize != sizeof(csize))
        rmg_error_handler (__FILE__,__LINE__,"error writing");

    wsize = write (fh, out, csize);
    if(wsize != csize)
        rmg_error_handler (__FILE__,__LINE__,"error writing");

    totalsize += csize;
    delete [] out;

}
