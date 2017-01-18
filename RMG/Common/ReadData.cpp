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
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


static void read_double (int fhand, double * rp, int count);
static void read_int (int fhand, int *ip, int count);

template void ReadData(char *, double *, double *, double *, Kpoint<double> **);
template void ReadData(char *, double *, double *, double *, Kpoint<std::complex<double> > **);

/* Reads the hartree potential, the wavefunctions, the */
/* compensating charges and various other things from a file. */
template <typename KpointType>
void ReadData (char *name, double * vh, double * rho, double * vxc, Kpoint<KpointType> ** Kptr)
{
    char newname[MAX_PATH + 200];
    int grid[3];
    int fine[3];
    int pe[3];
    int grid_size;
    int fgrid_size;
    int gamma;
    int nk, ik;
    int ns, is;

    /* wait until everybody gets here */
    /* my_barrier (); */
    MPI_Barrier(pct.img_comm);	

    /* Make the new output file name */
    rmg_printf("\nspin flag =%d\n", ct.spin_flag);
    
    int kstart = pct.kstart;
    if (ct.forceflag == BAND_STRUCTURE) kstart = 0;
    sprintf (newname, "%s_spin%d_kpt%d_gridpe%d", name, pct.spinpe, kstart, pct.gridpe);


    int fhand = open(newname, O_RDWR, S_IREAD | S_IWRITE);
    if (fhand < 0) {
        rmg_printf("Can't open data file %s", newname);
        rmg_error_handler(__FILE__, __LINE__, "Terminating.");
    }



    /* read grid info */
    read_int (fhand, grid, 3);
    if (grid[0] != Kptr[0]->G->get_NX_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong NX_GRID");
    if (grid[1] != Kptr[0]->G->get_NY_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong NY_GRID");
    if (grid[2] != Kptr[0]->G->get_NZ_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong NZ_GRID");

    /* read grid processor topology */
    read_int (fhand, pe, 3);
    if (pe[0] != Kptr[0]->G->get_PE_X())
        rmg_error_handler (__FILE__, __LINE__,"Wrong PE_X");
    if (pe[1] != Kptr[0]->G->get_PE_Y())
        rmg_error_handler (__FILE__, __LINE__,"Wrong PE_Y");
    if (pe[2] != Kptr[0]->G->get_PE_Z())
        rmg_error_handler (__FILE__, __LINE__,"Wrong PE_Z");

    grid_size = Kptr[0]->pbasis;

    /* read fine grid info */
    read_int (fhand, fine, 3);
    if (fine[0] != Kptr[0]->G->get_PX0_GRID(Kptr[0]->G->default_FG_RATIO) / Kptr[0]->G->get_PX0_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    if (fine[1] != Kptr[0]->G->get_PY0_GRID(Kptr[0]->G->default_FG_RATIO) / Kptr[0]->G->get_PY0_GRID(1))
        rmg_error_handler (__FILE__, __LINE__,"Wrong fine grid info");
    if (fine[2] != Kptr[0]->G->get_PZ0_GRID(Kptr[0]->G->default_FG_RATIO) / Kptr[0]->G->get_PZ0_GRID(1))
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
    if (ns > ct.num_states) {
        rmg_printf ("Wrong number of states: read %d from wave file, but ct.num_states is %d",ns, ct.num_states);
        rmg_error_handler (__FILE__, __LINE__,"Terminating.");
    }

    rmg_printf ("read_data: ns = %d\n", ns);


    /* read the hartree potential */
    read_double (fhand, vh, fgrid_size);
    rmg_printf ("read_data: read 'vh'\n");

    /* read density */
    read_double (fhand, rho, fgrid_size);
    rmg_printf ("read_data: read 'rho'\n");

    /* read Vxc */
    read_double (fhand, vxc, fgrid_size);
    rmg_printf ("read_data: read 'vxc'\n");

    /* read wavefunctions */
    {
        int wvfn_size = (gamma) ? grid_size : 2 * grid_size;
        double *tbuf = new double[wvfn_size];
        std::complex<double> *tptr;

        for (ik = 0; ik < ct.num_kpts_pe; ik++)
        {
            for (is = 0; is < ns; is++)
            {

                if(gamma == ct.is_gamma) {
                    read_double (fhand, (double *)Kptr[ik]->Kstates[is].psi, wvfn_size);
                }
                else {
                    // If wavefunctions on disk are complex but current calc is real then throw error
                    if(ct.is_gamma)
                         rmg_error_handler (__FILE__, __LINE__,"Can't convert complex wavefunctions to real.");

                    // Wavefunctions on disk are real but current calc is complex so convert them
                    ssize_t wanted = sizeof (double) * (ssize_t)wvfn_size;
                    ssize_t size = read (fhand, tbuf, wanted);

                    if(size != wanted)
                        rmg_error_handler (__FILE__, __LINE__,"error reading");

                    tptr = (std::complex<double> *)&Kptr[ik]->Kstates[is].psi;
                    for(int ix=0;ix < wvfn_size;ix++) tptr[ix] = std::complex<double>(tbuf[ix], 0.0);

                }

            }

            // for band structure calculation, just read wave functions for first kpoints
            if(ct.forceflag == BAND_STRUCTURE) break;
        }

        rmg_printf ("read_data: read 'wfns'\n");

        delete [] tbuf;
    }

    // If we have added unoccupied orbitals initialize them to a random state
    if(ct.num_states > ns) {

        for (ik = 0; ik < ct.num_kpts_pe; ik++){

            int PX0_GRID = Kptr[0]->G->get_PX0_GRID(1);
            int PY0_GRID = Kptr[0]->G->get_PY0_GRID(1);
            int PZ0_GRID = Kptr[0]->G->get_PZ0_GRID(1);

            int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
            double *tmp_psiR = new double[pbasis];
            double *tmp_psiI = new double[pbasis];

            double *xrand = new double[2 * Kptr[0]->G->get_NX_GRID(1)];
            double *yrand = new double[2 * Kptr[0]->G->get_NY_GRID(1)];
            double *zrand = new double[2 * Kptr[0]->G->get_NZ_GRID(1)];

            int factor = 2;
            if(ct.is_gamma) factor = 1;

            long int idum = 7493;
            int xoff = Kptr[0]->G->get_PX_OFFSET(1);
            int yoff = Kptr[0]->G->get_PY_OFFSET(1);
            int zoff = Kptr[0]->G->get_PZ_OFFSET(1);

            /* Initialize the random number generator */
            rand0 (&idum);

            for (int state = ns; state < ct.num_states; state++)
            {


                /* Generate x, y, z random number sequences */
                for (int idx = 0; idx < factor*Kptr[0]->G->get_NX_GRID(1); idx++)
                    xrand[idx] = rand0 (&idum) - 0.5;
                for (int idx = 0; idx < factor*Kptr[0]->G->get_NY_GRID(1); idx++)
                    yrand[idx] = rand0 (&idum) - 0.5;
                for (int idx = 0; idx < factor*Kptr[0]->G->get_NZ_GRID(1); idx++)
                    zrand[idx] = rand0 (&idum) - 0.5;


                int idx = 0;
                for (int ix = 0; ix < PX0_GRID; ix++)
                {

                    for (int iy = 0; iy < PY0_GRID; iy++)
                    {

                        for (int iz = 0; iz < PZ0_GRID; iz++)
                        {


                            tmp_psiR[idx] = xrand[xoff + ix] * 
                                yrand[yoff + iy] * 
                                zrand[zoff + iz];
                            tmp_psiR[idx] = tmp_psiR[idx] * tmp_psiR[idx];


                            if(!ct.is_gamma) {

                                tmp_psiI[idx] = xrand[Kptr[0]->G->get_NX_GRID(1) + xoff + ix] * 
                                    yrand[Kptr[0]->G->get_NY_GRID(1) + yoff + iy] * 
                                    zrand[Kptr[0]->G->get_NZ_GRID(1) + zoff + iz];
                                tmp_psiI[idx] = tmp_psiI[idx] * tmp_psiI[idx];

                            }

                            idx++;

                        }               /* end for */
                    }                   /* end for */
                }                       /* end for */

                // Copy data from tmp_psi into orbital storage
                for(idx = 0;idx < pbasis;idx++) {
                    Kptr[ik]->Kstates[state].psi[idx] = tmp_psiR[idx];
                }
                if(typeid(KpointType) == typeid(std::complex<double>)) {
                    for(idx = 0;idx < pbasis;idx++) {
                        double *a = (double *)&Kptr[ik]->Kstates[state].psi[idx];
                        if(!ct.is_gamma)
                            a[1] = tmp_psiI[idx];

                    }

                }

                // Hit the orbital with the right hand mehrstellen operator which should smooth it a bit
                //        CPP_app_cir_driver (this->L, this->T, this->Kstates[state].psi, this->Kstates[state].psi, PX0_GRID, PY0_GRID, PZ0_GRID, APP_CI_FOURTH);

            }                           /* end for */

            delete [] zrand;
            delete [] yrand;
            delete [] xrand;
            delete [] tmp_psiI;
            delete [] tmp_psiR;
        }                           /* end for */
        if(ct.forceflag == BAND_STRUCTURE) return;
    }


    /* read state occupations */
    {
        double *occ = new double[nk * ct.num_states]();

        read_double (fhand, occ, (nk * ns));

        rmg_printf ("read_data: read 'occupations'\n"); 


        if (ct.forceflag != BAND_STRUCTURE)
        {
            double occ_total = 0.0; 

            for (ik = 0; ik < nk; ik++)
                for (is = 0; is < ct.num_states; is++)
                {
                    occ_total += ( Kptr[ik]->Kstates[is].occupation[0] = occ[ik * ns + is] );
                }



            /* 
               since we are using floats on the data file the precision
               of the occupations is worse than 1e-10 required by the fill() routine
               therefore we need to 'renormalize' the occupations so
               that they add up to an integer
               it's a hack I know, but whatever... not a biggie
               */

            {
                double iocc_total = (double) (int) (occ_total + 0.5);
                double fac = iocc_total / occ_total;

                for (ik = 0; ik < nk; ik++)
                    for (is = 0; is < ct.num_states; is++)
                    {
                        Kptr[ik]->Kstates[is].occupation[0] *= fac;
                    }
                /* end of normalization*/
            }

        }             /* end if */

        delete [] occ;

    }           /* end of read occupations */




    /* read state eigenvalues, not needed really */
    {

        /* Read eigenvalue in pairwised case, while in polarized case, 
         * it's the eigenvalue for proceesor's own spin  */ 
        for (ik = 0; ik < nk; ik++)
            for (is = 0; is < ns; is++)
            {
                read_double (fhand, &Kptr[ik]->Kstates[is].eig[0], 1);
            }

        rmg_printf ("read_data: read 'eigenvalues'\n");

    }      /* end of read eigenvalues */





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
    int size;

    size = count * sizeof (int);

    if (size != read (fhand, ip, size))
        rmg_error_handler (__FILE__, __LINE__,"error reading");

}


/******/
