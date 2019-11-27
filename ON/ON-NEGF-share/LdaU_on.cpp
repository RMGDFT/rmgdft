/*
 *
 * Copyright (c) 2019, Wenchang Lu
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

#include "LdaU_on.h"



#include <float.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include "const.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>

#include "LocalObject.h"
#include "BaseGrid.h"
#include "prototypes_on.h"
#include "blas.h"
#include "transition.h"
#include "GpuAlloc.h"


LdaU_on::~LdaU_on(void)
{
    delete []  this->ns_occ;
    delete []  this->Upsi_mat;
    delete this->AtomicOrbital;

}

LdaU_on::LdaU_on(LocalObject<double> &LO, BaseGrid &BG)
{
    InitDelocalizedOrbital ();

    int tot_orbitals_ldaU = 0;

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        SPECIES *sp = &Species[Atoms[ion].species];
        double Ueff = sp->Hubbard_U / 2.0;       // FIXME: Have to deal with more complicated cases later
        

        for (int ip = 0; ip < sp->num_orbitals; ip++)
        {
            if(sp->awave_is_ldaU[ip])
            {
                this->Hubbard_U.push_back(Ueff);
                tot_orbitals_ldaU++;
            }
        }
    }

    this->Ehub = 0.0;
    this->Ecorrect = 0.0;

    if(tot_orbitals_ldaU == 0) return;
    int *ixmin=NULL, *iymin=NULL, *izmin=NULL, *dimx=NULL, *dimy=NULL, *dimz=NULL;
    int density = 1;
    this->tot_orbitals_ldaU = tot_orbitals_ldaU;
    this->AtomicOrbital = new LocalObject<double>(tot_orbitals_ldaU, ixmin, iymin, izmin,
            dimx, dimy, dimz, 1, BG, density, pct.grid_comm);
    this->AtomicOrbital->GetAtomicOrbitals(ct.num_ions, BG);
    size_t size = tot_orbitals_ldaU * LO.num_tot * sizeof(double);
    this->Upsi_mat = (double *)GpuMallocManaged(size);
    size = this->AtomicOrbital->num_thispe * LO.num_thispe * sizeof(double);
    this->Upsi_mat_local = (double *)GpuMallocManaged(size);

    // now only works for atomic orbital expanded in whole space.
    assert(this->AtomicOrbital->num_thispe == tot_orbitals_ldaU);
    this->ns_occ = new double[tot_orbitals_ldaU * tot_orbitals_ldaU];
    for(int i = 0; i < tot_orbitals_ldaU * tot_orbitals_ldaU; i++) this->ns_occ[i] = 0.0;
    for(int i = 0; i < tot_orbitals_ldaU; i++) this->ns_occ[i * tot_orbitals_ldaU+i] = 1.0;
    if(ct.spin_flag) 
        for(int i = 0; i < tot_orbitals_ldaU; i++) this->ns_occ[i * tot_orbitals_ldaU+i] = 0.5;

    for(int i = 0; i <  this->AtomicOrbital->num_thispe * LO.num_thispe; i++) this->Upsi_mat_local[i] = 0.0;

    //Upsi_mat_local = new double[LocalAtomicOrbital->num_thispe * LocalOrbital->num_thispe];

}


void LdaU_on::calc_ns_occ(LocalObject<double> &LocalOrbital, double *mat_X, BaseGrid &BG)
{
    LO_x_LO(*this->AtomicOrbital, LocalOrbital, this->Upsi_mat_local, BG);

    mat_local_to_glob(this->Upsi_mat_local, this->Upsi_mat, *this->AtomicOrbital, LocalOrbital, 
            0, this->AtomicOrbital->num_tot, 0, LocalOrbital.num_tot);

    int nldaU = this->tot_orbitals_ldaU;
    int norb = LocalOrbital.num_tot;
    double *mat_X_glob = new double[norb*norb];
    int nmax = std::max(nldaU, norb);
    double *temA = new double[nldaU * nmax];

    mat_dist_to_global(mat_X, pct.desca, mat_X_glob);
    double one = 1.0, zero = 0.0;
    dgemm("N", "N", &nldaU, &norb, &norb, &one, this->Upsi_mat, &nldaU, mat_X_glob, &norb,
            &zero, temA, &nldaU);
    dgemm("N", "T", &nldaU, &nldaU, &norb, &one, temA, &nldaU, this->Upsi_mat, &nldaU,
            &zero, this->ns_occ, &nldaU);

    if(!ct.spin_flag)
        for(int i = 0; i < nldaU*nldaU; i++) this->ns_occ[i] *=0.5;
    
    mat_global_to_local(*this->AtomicOrbital, LocalOrbital, this->Upsi_mat, this->Upsi_mat_local);
    for(int j = 0; j < LocalOrbital.num_thispe; j++)
        for(int i = 0; i < nldaU; i++)
            this->Upsi_mat_local[j * nldaU + i] *= this->Hubbard_U[i] * (1.0-2.0*this->ns_occ[i * nldaU + i]);

    this->Ehub = 0.0;
    this->Ecorrect = 0.0;
    for(int i = 0; i < nldaU; i++)
    {
        this->Ehub += this->Hubbard_U[i] * this->ns_occ[i * nldaU + i] * (1.0-this->ns_occ[i * nldaU + i]);
        this->Ecorrect += this->Hubbard_U[i] * this->ns_occ[i * nldaU + i] * this->ns_occ[i * nldaU + i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &this->Ehub, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);
    MPI_Allreduce(MPI_IN_PLACE, &this->Ecorrect, 1, MPI_DOUBLE, MPI_SUM, pct.spin_comm);


    //  for(int i = 0; i < nldaU; i++)
    //  {
    //      printf("\n");
    //for(int j = 0; j < nldaU; j++)
    //{
    //    printf("\n %f ", this->ns_occ[j*nldaU + j]);
    //} 

    //  }

    //      printf("\n");
    //  for(int i = 0; i < norb; i++)
    //  {
    //      printf("\n");
    //      for(int j = 0; j < nldaU; j++)
    //      {
    //          printf(" %f ", this->Upsi_mat[i*nldaU + j]);
    //      } 

    //  }

    delete [] mat_X_glob;
    delete [] temA;
}

void LdaU_on::app_vhubbard(LocalObject<double> &HLO, BaseGrid &BG)
{

    int pbasis = BG.get_P0_BASIS(1);
    int nldaU = this->AtomicOrbital->num_thispe;
    int norb = HLO.num_thispe;
    double one = 1.0;
    dgemm("N", "N", &pbasis, &norb, &nldaU, &one, this->AtomicOrbital->storage_proj, &pbasis, 
            this->Upsi_mat_local, &nldaU, &one, HLO.storage_proj, &pbasis);

}
void LdaU_on::WriteLdaU(std::string file_prefix, LocalObject<double> &LO)
{
    std::string filename = file_prefix + "_spin"+std::to_string(pct.spinpe)+".LdaU";
    std::ofstream fhand(filename.c_str(), std::ofstream::binary);
    size_t size= this->tot_orbitals_ldaU * this->tot_orbitals_ldaU * sizeof(double);
    fhand.write((char *)this->ns_occ, size);

    size = LO.num_tot * this->tot_orbitals_ldaU * sizeof(double);
    fhand.write((char *)this->Upsi_mat, size);
    fhand.close();

}
void LdaU_on::ReadLdaU(std::string file_prefix, LocalObject<double> &LO)
{
    std::string filename = file_prefix + "_spin"+std::to_string(pct.spinpe)+".LdaU";
    std::ifstream fhand(filename.c_str(), std::ifstream::binary);
    int nldaU = this->tot_orbitals_ldaU;

    size_t size= nldaU * nldaU * sizeof(double);
    fhand.read((char *)this->ns_occ, size);

    size = LO.num_tot * nldaU * sizeof(double);
    fhand.read((char *)this->Upsi_mat, size);
    fhand.close();

    mat_global_to_local(*this->AtomicOrbital, LO, this->Upsi_mat, this->Upsi_mat_local);

    for(int j = 0; j < LO.num_thispe; j++)
        for(int i = 0; i < nldaU; i++)
            this->Upsi_mat_local[j * nldaU + i] *= this->Hubbard_U[i] * (1.0-2.0*this->ns_occ[i * nldaU + i]);
}


