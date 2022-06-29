/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iterator>
#include <omp.h>


#include "const.h"
#include "Exxbase.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "transition.h"
#include "rmgtypedefs.h"
#include "pe_control.h"
#include "GlobalSums.h"
#include "vdw_Grimme.h"


using namespace vdw_Grimme;
double vdw_d2_energy(Lattice &L,  std::vector<ION>& Atoms_in)
{
    double Rij, energy, c6i, c6j, c6ij, dist, f6d;
    double xtal[3], xtal_tot[3];
    energy = 0.0;
    int maxnx = 1 + rcut/L.celldm[0];
    int maxny = 1 + rcut/(L.celldm[0] * L.celldm[1]);
    int maxnz = 1 + rcut/(L.celldm[0] * L.celldm[2]);
    for(size_t i = pct.gridpe; i < Atoms_in.size(); i+=pct.grid_npes)
        for(size_t j = 0; j < Atoms_in.size(); j++)
        {

            c6i = C6[Atoms_in[i].symbol];
            c6j = C6[Atoms_in[j].symbol];
            // C6 is in unit of Ry/au^6, change to Ha/au^6
            c6ij = std::sqrt(c6i * c6j) * 0.5; 
            Rij = Ri[Atoms_in[i].symbol] + Ri[Atoms_in[j].symbol];
            xtal[0] = Atoms_in[i].xtal[0] - Atoms_in[j].xtal[0];
            xtal[1] = Atoms_in[i].xtal[1] - Atoms_in[j].xtal[1];
            xtal[2] = Atoms_in[i].xtal[2] - Atoms_in[j].xtal[2];
            for(int ix = -maxnx; ix <= maxnx; ix++)
                for(int iy = -maxny; iy <= maxny; iy++)
                    for(int iz = -maxnz; iz <= maxnz; iz++)
                    {
                        if(i == j && ix == 0 && iy == 0 && iz == 0) continue;
                        xtal_tot[0] = xtal[0] + ix;
                        xtal_tot[1] = xtal[1] + iy;
                        xtal_tot[2] = xtal[2] + iz;
                        dist = L.metric(xtal_tot);
                        if(dist < rcut)
                        {
                            f6d = scale6/(1.0 + std::exp(-damp * (dist/Rij -1.0)));
                            energy -=  0.5 * c6ij/std::pow(dist, 6) * f6d;
                        }
                    }

        }

    GlobalSums(&energy, 1, pct.grid_comm);

    return energy;
}

void vdw_d2_forces(Lattice &L,  std::vector<ION>& Atoms_in, double *forces)
{
    double Rij, c6i, c6j, c6ij, dist, f6d, exptmp;
    double xtal[3], xtal_tot[3], xcrt[3];
    for(size_t i = 0; i < 3*Atoms_in.size(); i++)
        forces[i] = 0.0;

    int maxnx = 1 + rcut/L.celldm[0];
    int maxny = 1 + rcut/(L.celldm[0] * L.celldm[1]);
    int maxnz = 1 + rcut/(L.celldm[0] * L.celldm[2]);
    for(size_t i = pct.gridpe; i < Atoms_in.size(); i+=pct.grid_npes)
        for(size_t j = 0; j < Atoms_in.size(); j++)
        {

            c6i = C6[Atoms_in[i].symbol];
            c6j = C6[Atoms_in[j].symbol];
            // C6 is in unit of Ry/au^6, change to Ha/au^6
            c6ij = std::sqrt(c6i * c6j) * 0.5; 
            Rij = Ri[Atoms_in[i].symbol] + Ri[Atoms_in[j].symbol];
            xtal[0] = Atoms_in[i].xtal[0] - Atoms_in[j].xtal[0];
            xtal[1] = Atoms_in[i].xtal[1] - Atoms_in[j].xtal[1];
            xtal[2] = Atoms_in[i].xtal[2] - Atoms_in[j].xtal[2];
            for(int ix = -maxnx; ix <= maxnx; ix++)
                for(int iy = -maxny; iy <= maxny; iy++)
                    for(int iz = -maxnz; iz <= maxnz; iz++)
                    {
                        if(i == j && ix == 0 && iy == 0 && iz == 0) continue;
                        xtal_tot[0] = xtal[0] + ix;
                        xtal_tot[1] = xtal[1] + iy;
                        xtal_tot[2] = xtal[2] + iz;
                        dist = L.metric(xtal_tot);
                        L.to_cartesian(xtal_tot, xcrt);
                        if(dist < rcut)
                        {
                            exptmp = std::exp(-damp * (dist/Rij -1.0));
                            f6d = scale6/(1.0 + exptmp );
                            double tem = c6ij/std::pow(dist, 6) * f6d * 
                                (6.0/dist - damp/Rij * exptmp/(1.0 + exptmp)) /dist;
                            forces[i * 3+ 0] += tem * xcrt[0];
                            forces[i * 3+ 1] += tem * xcrt[1];
                            forces[i * 3+ 2] += tem * xcrt[2];
                        }
                    }

        }

}

void vdw_d2_stress(Lattice &L,  std::vector<ION>& Atoms_in, double *stress_d2)
{
    double Rij, c6i, c6j, c6ij, dist, f6d, exptmp;
    double xtal[3], xtal_tot[3], xcrt[3];
    for(size_t i = 0; i < 9; i++)
        stress_d2[i] = 0.0;

    int maxnx = 1 + rcut/L.celldm[0];
    int maxny = 1 + rcut/(L.celldm[0] * L.celldm[1]);
    int maxnz = 1 + rcut/(L.celldm[0] * L.celldm[2]);
    for(size_t i = pct.gridpe; i < Atoms_in.size(); i+=pct.grid_npes)
        for(size_t j = 0; j < Atoms_in.size(); j++)
        {

            c6i = C6[Atoms_in[i].symbol];
            c6j = C6[Atoms_in[j].symbol];
            // C6 is in unit of Ry/au^6, change to Ha/au^6
            c6ij = std::sqrt(c6i * c6j) * 0.5; 
            Rij = Ri[Atoms_in[i].symbol] + Ri[Atoms_in[j].symbol];
            xtal[0] = Atoms_in[i].xtal[0] - Atoms_in[j].xtal[0];
            xtal[1] = Atoms_in[i].xtal[1] - Atoms_in[j].xtal[1];
            xtal[2] = Atoms_in[i].xtal[2] - Atoms_in[j].xtal[2];
            for(int ix = -maxnx; ix <= maxnx; ix++)
                for(int iy = -maxny; iy <= maxny; iy++)
                    for(int iz = -maxnz; iz <= maxnz; iz++)
                    {
                        if(i == j && ix == 0 && iy == 0 && iz == 0) continue;
                        xtal_tot[0] = xtal[0] + ix;
                        xtal_tot[1] = xtal[1] + iy;
                        xtal_tot[2] = xtal[2] + iz;
                        dist = L.metric(xtal_tot);
                        L.to_cartesian(xtal_tot, xcrt);
                        if(dist < rcut)
                        {
                            exptmp = std::exp(-damp * (dist/Rij -1.0));
                            f6d = scale6/(1.0 + exptmp );
                            double tem = c6ij/std::pow(dist, 6) * f6d * 
                                (6.0/dist - damp/Rij * exptmp/(1.0 + exptmp)) /dist;

                            for(int id1 = 0; id1 < 3; id1++)
                            for(int id2 = 0; id2 < 3; id2++)
                                stress_d2[id1 * 3 + id2] += tem * xcrt[id1] * xcrt[id2];

                        }
                    }

        }

    for(int idx = 0; idx < 9; idx++) stress_d2[idx] = -stress_d2[idx]/(2.0 * L.omega);
    GlobalSums(stress_d2, 9, pct.grid_comm);

}



