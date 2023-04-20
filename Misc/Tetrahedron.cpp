/*
 *
 * Copyright 2023 The RMG Project Developers. See the COPYRIGHT file 
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

#include <math.h>
#include "transition.h"
#include "rmg_error.h"
#include "rmg_mangling.h"
#include "Symmetry.h"
#include "Tetrahedron.h"

Tetrahedron *Tetra;

template double Tetrahedron::FillTetra (Kpoint<double> **Kptr);
template double Tetrahedron::FillTetra (Kpoint<std::complex<double>> **Kptr);

extern "C" void tetra_set_method(int *method); 
extern "C" void tetra_init( int *nsym, int *s, bool *time_reversal, int *t_rev, 
                 double *at, double *bg, int *npk,
                 int *k1, int *k2, int *k3, 
                 int *nk1, int *nk2, int *nk3,
                 int *nks, double *xk);

extern "C" void opt_tetra_init( int *nsym, int *s, bool *time_reversal, int *t_rev, 
                double *at, double *bg, int *npk,
                int *k1, int *k2, int *k3,
                int *nk1, int *nk2, int *nk3,
                int *nks, double *xk, int *kstep );

extern "C" void tetra_weights( int *nks, int *nspin, int *nbnd, double *nelec, double *et,
                          double *ef, double *wg, int *is, int *isk );
extern "C" void tetra_weights_only( int *nks, int *nspin, int *is, int *isk, int *nbnd,
                          double *nelec, double *et, double *ef, double *wg); 
extern "C" void opt_tetra_weights( int *nks, int *nspin, int *nbnd, double *nelec, double *et, double *ef,
                              double *wg, int *is, int *isk );
extern "C" void opt_tetra_weights_only( int *nks, int *nspin, int *nbnd, double *et, double *ef,
                                   double *wg, int *is, int *isk ); 
extern "C" void deallocate_tetra(void); 

Tetrahedron::Tetrahedron(void)
{
    tetra_set_method(&ct.tetra_method);
    int *so = Rmg_Symm->sym_rotate.data();
    int s[48][3][3];
    double at[9], bg[9];
    double *xk = new double[3*ct.num_kpts];
    int *t_rev = new int[48]();
    double alat = Rmg_L.celldm[0];

    for(int is=0;is < ct.nsym;is++) t_rev[is] = (int)Rmg_Symm->time_rev[is];

    for(int is=0;is < ct.nsym;is++)
    {
        for(int i=0;i < 3;i++)
        {
            for(int j=0;j < 3;j++)
            {
                s[is][i][j] = so[is*9 + 3*i + j];
            }
        }
    }

    // Setup lattice vectors and reciprocal lattice vectors for f90 routines
    for(int i=0;i < 3;i++)
    {
        at[i] = Rmg_L.a0[i] / alat;
        at[i+3] = Rmg_L.a1[i] / alat;
        at[i+6] = Rmg_L.a2[i] / alat;
        bg[i] = alat * Rmg_L.b0[i];
        bg[i+3] = alat * Rmg_L.b1[i];
        bg[i+6] = alat * Rmg_L.b2[i];
    }

    double kunit = twoPI /Rmg_L.celldm[0];
    for (int kpt = 0; kpt < ct.num_kpts; kpt++)
    {
        xk[kpt*3] = ct.kp[kpt].kvec[0]/kunit;
        xk[kpt*3+1] = ct.kp[kpt].kvec[1]/kunit;
        xk[kpt*3+2] = ct.kp[kpt].kvec[2]/kunit;
    }

    if(ct.tetra_method == 0)
    {
        tetra_init(&ct.nsym, &s[0][0][0], &Rmg_Symm->time_reversal, t_rev,
               at, bg, &ct.num_kpts,
               &ct.kpoint_is_shift[0], &ct.kpoint_is_shift[1], &ct.kpoint_is_shift[2],
               &ct.kpoint_mesh[0], &ct.kpoint_mesh[1], &ct.kpoint_mesh[2],
               &ct.num_kpts, xk);
    }
    else
    {
        int ione = 1; 
        opt_tetra_init(&ct.nsym, &s[0][0][0], &Rmg_Symm->time_reversal, t_rev,
               at, bg, &ct.num_kpts,
               &ct.kpoint_is_shift[0], &ct.kpoint_is_shift[1], &ct.kpoint_is_shift[2],
               &ct.kpoint_mesh[0], &ct.kpoint_mesh[1], &ct.kpoint_mesh[2],
               &ct.num_kpts, xk, &ione);
    }

    delete [] t_rev;
    delete [] xk;

}

template <typename KpointType>
double Tetrahedron::FillTetra (Kpoint<KpointType> **Kptr)
{
    size_t knum_states = ct.num_kpts*ct.num_states;
    int nspin = (ct.spin_flag + 1);

    std::vector<double> eigs(nspin * ct.num_kpts * ct.num_states);
    std::vector<double> wg(nspin * ct.num_kpts * ct.num_states);
    eigs.assign(nspin * ct.num_kpts * ct.num_states, 0.0);

    // Need to collect eigenvalues into a contiguous array. If spin polarized first
    // ct.num_states*ct.num_kpts is spin up and next ct.num_states*ct.num_kpts is down
    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
        int kpt1 = pct.kstart + kpt;  // absolute kpoint index
        for(int st = 0;st < ct.num_states;st++)
        {   
            eigs[kpt1*ct.num_states + st] = Kptr[kpt]->Kstates[st].eig[0];
            if(nspin == 2)
            {
                eigs[kpt1*ct.num_states + st + ct.num_kpts * ct.num_states] = Kptr[kpt]->Kstates[st].eig[1];
            }
        }
    }
    MPI_Allreduce (MPI_IN_PLACE, eigs.data(), nspin * ct.num_kpts * ct.num_states, MPI_DOUBLE, MPI_SUM, pct.kpsub_comm);
    MPI_Bcast(eigs.data(), nspin * ct.num_kpts * ct.num_states, MPI_DOUBLE, 0, pct.kpsub_comm);

    int is = 0;
    std::vector<int> isk(nspin*ct.num_kpts);
    isk.assign(ct.num_kpts, 1);
    if(nspin == 2) for(int ix=0;ix < ct.num_kpts;ix++) isk[ix+ct.num_kpts] = 2;
    int nks;
    nks = nspin*ct.num_kpts;
    // The next bit is used to match the QE routine expectation for the noncollinear case
    int nspin1 = nspin;
    if(ct.noncoll_factor) nspin1 = 4;
    if(ct.tetra_method == 0)  // Bloechl
    {
        tetra_weights( &nks, &nspin1, &ct.num_states, &ct.nel,
               eigs.data(), &this->ef, wg.data(), &is, isk.data() );
    }
    else  // Linear=1, Optimized=2
    {
        opt_tetra_weights( &nks, &nspin1, &ct.num_states, &ct.nel,
               eigs.data(), &this->ef, wg.data(), &is, isk.data() );
    }

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
        int kpt1 = kpt + pct.kstart;
        for(int st = 0;st < ct.num_states;st++)
        {
            Kptr[kpt]->Kstates[st].occupation[0] = wg[kpt1*ct.num_states + st]/ct.kp[kpt1].kweight;
            if(nspin == 2)
                Kptr[kpt]->Kstates[st].occupation[1] = wg[kpt1*ct.num_states + st + ct.num_kpts * ct.num_states]/ct.kp[kpt1].kweight;
            //printf("WEIGHT for KPT %d  STATE %d = %f\n",kpt, st, Kptr[kpt]->Kstates[st].occupation[0]);
        }
    }

    //printf("EFERMI  %f\n", this->ef);fflush(NULL);
    return this->ef;
}

Tetrahedron::~Tetrahedron(void)
{
    deallocate_tetra();
}
