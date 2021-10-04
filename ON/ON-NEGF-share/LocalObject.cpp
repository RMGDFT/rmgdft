/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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

#include "GpuAlloc.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex>
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "RmgParallelFft.h"

#include "LocalObject.h"
#include "GpuAlloc.h"
#include "blas.h"


template LocalObject<double>::~LocalObject(void);
template LocalObject<std::complex<double>>::~LocalObject(void);

template LocalObject<double>::LocalObject(int, int*, int*, int*, int*, int*, int*, bool, BaseGrid&, int, MPI_Comm);
template LocalObject<std::complex<double>>::LocalObject(int, int*, int*, int*, int*, int*, int*, bool, BaseGrid&, int, MPI_Comm);
template <class KpointType> LocalObject<KpointType>::LocalObject(int num_objects, 
        int *ixmin, int *iymin, int *izmin, int *dimx, int *dimy, int *dimz, bool delocalized,
        BaseGrid &BG, int density, MPI_Comm comm)
{
    this->num_tot = num_objects;
    this->num_thispe = 0;
    this->comm = comm;
    this->density = density;
    this->delocalized = delocalized;

    this->index_global_to_proj = new int[num_objects];
    this->index_proj_to_global = new int[num_objects];
    this->ixmin = new int[num_objects];
    this->iymin = new int[num_objects];
    this->izmin = new int[num_objects];
    this->dimx = new int[num_objects];
    this->dimy = new int[num_objects];
    this->dimz = new int[num_objects];

    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    int PX_OFFSET = BG.get_PX_OFFSET(density);
    int PY_OFFSET = BG.get_PY_OFFSET(density);
    int PZ_OFFSET = BG.get_PZ_OFFSET(density);
    int NX_GRID = BG.get_NX_GRID(density);
    int NY_GRID = BG.get_NY_GRID(density);
    int NZ_GRID = BG.get_NZ_GRID(density);
    int ilow = PX_OFFSET;
    int ihigh = ilow + PX0_GRID;
    int jlow = PY_OFFSET;
    int jhigh = jlow + PY0_GRID;
    int klow = PZ_OFFSET;
    int khigh = klow + PZ0_GRID;

    this->pbasis = P0_BASIS;
    if(this->delocalized)
    {
        this->num_thispe = this->num_tot;
        for (int i = 0; i < this->num_tot; i++)
        {
            this->index_global_to_proj[i] = i;
            this->index_proj_to_global[i] = i;
            this->ixmin[i] = 0;
            this->iymin[i] = 0; 
            this->izmin[i] = 0;
            this->dimx[i] = NX_GRID;
            this->dimy[i] = NY_GRID; 
            this->dimz[i] = NZ_GRID; 
        }
        this->max_dimx = NX_GRID;
        this->max_dimy = NY_GRID;
        this->max_dimz = NZ_GRID;

    }
    else
    {
        this->max_dimx = *std::max_element(dimx, dimx+this->num_tot);
        this->max_dimy = *std::max_element(dimy, dimy+this->num_tot);
        this->max_dimz = *std::max_element(dimz, dimz+this->num_tot);
        for (int i = 0; i < this->num_tot; i++)
        {
            this->index_global_to_proj[i] = -1;
            this->index_proj_to_global[i] = -1;
            this->ixmin[i] = ixmin[i];
            this->iymin[i] = iymin[i];
            this->izmin[i] = izmin[i];
            this->dimx[i] = dimx[i];
            this->dimy[i] = dimy[i];
            this->dimz[i] = dimz[i];

        }

        bool xo, yo, zo;
        for (int i = 0; i < this->num_tot; i++)
        {
            xo = false;
            if(ixmin[i] < 0) 
            {
                if (ilow < ixmin[i] + dimx[i]) xo = true;
                if (ihigh > ixmin[i] + NX_GRID) xo = true;
            }
            else if(ixmin[i] +dimx[i] > NX_GRID)
            {
                if (ilow < ixmin[i] + dimx[i] - NX_GRID) xo = true;
                if (ihigh > ixmin[i]) xo = true;
            }         
            {
                if( (ilow < ixmin[i] + dimx[i]) && (ihigh >ixmin[i]) ) xo = true;
            }

            yo = false;
            if(iymin[i] < 0) 
            {
                if (jlow < iymin[i] + dimy[i]) yo = true;
                if (jhigh > iymin[i] + NY_GRID) yo = true;
            }
            else if(iymin[i] +dimy[i] > NY_GRID)
            {
                if (jlow < iymin[i] + dimy[i] - NY_GRID) yo = true;
                if (jhigh > iymin[i]) yo = true;
            }         
            {
                if( (jlow < iymin[i] + dimy[i]) && (jhigh >iymin[i]) ) yo = true;
            }

            zo = false;
            if(izmin[i] < 0) 
            {
                if (klow < izmin[i] + dimz[i]) zo = true;
                if (khigh > izmin[i] + NZ_GRID) zo = true;
            }
            else if(izmin[i] +dimz[i] > NZ_GRID)
            {
                if (klow < izmin[i] + dimz[i] - NZ_GRID) zo = true;
                if (khigh > izmin[i]) zo = true;
            }         
            {
                if( (klow < izmin[i] + dimz[i]) && (khigh >izmin[i]) ) zo = true;
            }

            if(xo && yo && zo)
            {
                this->index_global_to_proj[i] = this->num_thispe;
                this->index_proj_to_global[this->num_thispe] = i;
                this->num_thispe++;

            }

        }
    }

    size_t size = this->num_thispe * P0_BASIS *sizeof(KpointType) +8;
    this->storage_size = size;
    this->storage_cpu = (KpointType *) RmgMallocHost(size);
    gpuMalloc((void **)&this->storage_gpu, size);
    this->storage_ptr = MemoryPtrHostDevice(this->storage_cpu, this->storage_gpu);

}

template <class KpointType> LocalObject<KpointType>::~LocalObject(void)
{
    delete [] this->index_global_to_proj;
    delete [] this->index_proj_to_global;
    delete [] this->ixmin;
    delete [] this->iymin;
    delete [] this->izmin;
    delete [] this->dimx;
    delete [] this->dimy;
    delete [] this->dimz;
    RmgFreeHost(this->storage_cpu);
    gpuFree(this->storage_gpu);

}

template void LocalObject<double>::ReadOrbitalsFromSingleFiles(std::string filename, BaseGrid &BG);
template void LocalObject<std::complex<double>>::ReadOrbitalsFromSingleFiles(std::string filename, BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::ReadOrbitalsFromSingleFiles(std::string filename, BaseGrid &BG)
{

    int density = this->density;
    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    int PX_OFFSET = BG.get_PX_OFFSET(density);
    int PY_OFFSET = BG.get_PY_OFFSET(density);
    int PZ_OFFSET = BG.get_PZ_OFFSET(density);
    int NX_GRID = BG.get_NX_GRID(density);
    int NY_GRID = BG.get_NY_GRID(density);
    int NZ_GRID = BG.get_NZ_GRID(density);
    int ilow = PX_OFFSET;
    int ihigh = ilow + PX0_GRID;
    int jlow = PY_OFFSET;
    int jhigh = jlow + PY0_GRID;
    int klow = PZ_OFFSET;
    int khigh = klow + PZ0_GRID;

    int fhand;

    for(int idx = 0; idx < this->num_thispe * P0_BASIS; idx++) this->storage_cpu[idx] = 0.0;
    for(int st = 0; st < this->num_thispe; st++)
    {

        int st_glob = this->index_proj_to_global[st];
        if(st_glob < 0) continue;
        size_t size = this->dimx[st_glob] * this->dimy[st_glob] * this->dimz[st_glob] *sizeof(KpointType);

        std::string newname= filename + "_spin"+std::to_string(pct.spinpe)+".orbit_"+std::to_string(st_glob);
        fhand = open(newname.c_str(), O_RDWR, S_IREAD | S_IWRITE);
        if(fhand < 0)
        {
            printf ("\n %s \n", newname.c_str());
            fflush(NULL);
            exit(0);
        }

        KpointType *psi = new KpointType[size/sizeof(KpointType)];
        size_t nbytes = read(fhand, psi, size);
        if (nbytes != size)
        {
            printf ("\n read %d is different from %d for state %d", (int) nbytes, (int)size, st_glob);
            fflush(NULL);
            exit(0);
        }
        close(fhand);

        for(int ix = 0; ix < this->dimx[st_glob]; ix++)
            for(int iy = 0; iy < this->dimy[st_glob]; iy++)
                for(int iz = 0; iz < this->dimz[st_glob]; iz++)
                {
                    int ixx = ix + this->ixmin[st_glob]; 
                    int iyy = iy + this->iymin[st_glob]; 
                    int izz = iz + this->izmin[st_glob]; 

                    if (ixx < 0) ixx += NX_GRID;
                    if (iyy < 0) iyy += NY_GRID;
                    if (izz < 0) izz += NZ_GRID;
                    if (ixx >= NX_GRID) ixx -=NX_GRID;
                    if (iyy >= NY_GRID) iyy -=NY_GRID;
                    if (izz >= NZ_GRID) izz -=NZ_GRID;

                    if(ixx >=ilow && ixx < ihigh && iyy >=jlow && iyy <jhigh && izz >=klow && izz < khigh)
                    {
                        int idx = (ixx - ilow) * PY0_GRID *PZ0_GRID + (iyy-jlow)*PZ0_GRID + izz-klow;
                        int idx0 = ix * this->dimy[st_glob] * this->dimz[st_glob] + iy * this->dimz[st_glob] + iz;
                        this->storage_cpu[st * P0_BASIS + idx] = psi[idx0];
                    }

                }

        delete []psi;
    }
}

template void LocalObject<double>::ReadProjectors(int, int, int *,BaseGrid &BG);
template void LocalObject<std::complex<double>>::ReadProjectors(int, int, int*, BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::ReadProjectors(int num_ions, int max_nlpoints,
        int *num_prj_perion, BaseGrid &BG)
{

    int density = this->density;
    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    int PX_OFFSET = BG.get_PX_OFFSET(density);
    int PY_OFFSET = BG.get_PY_OFFSET(density);
    int PZ_OFFSET = BG.get_PZ_OFFSET(density);
    int NX_GRID = BG.get_NX_GRID(density);
    int NY_GRID = BG.get_NY_GRID(density);
    int NZ_GRID = BG.get_NZ_GRID(density);
    int ilow = PX_OFFSET;
    int ihigh = ilow + PX0_GRID;
    int jlow = PY_OFFSET;
    int jhigh = jlow + PY0_GRID;
    int klow = PZ_OFFSET;
    int khigh = klow + PZ0_GRID;

    int fhand;
    int proj_count = 0;
    for (int i = 0; i < this->num_thispe; i++)
    {
        for (int j = 0; j < P0_BASIS; j++)
            this->storage_cpu[i * P0_BASIS + j] = 0.0;
    }


    for(int ion = 0; ion < num_ions; ion++)
    {
        if ( this->index_global_to_proj[proj_count] == -1)
        {
            proj_count += num_prj_perion[ion];
            continue;
        }


        int proj_local_index;

        std::string newname;
        newname = std::string("PROJECTORS/ion_") + std::to_string(ion);

        fhand = open(newname.c_str(), O_RDWR, S_IREAD | S_IWRITE);
        ssize_t size = (ssize_t)num_prj_perion[ion] * (ssize_t)max_nlpoints;
        double *beta = new double[size];
        ssize_t read_size = read(fhand, beta, size * sizeof(double));
        if(read_size != (ssize_t)(size * sizeof(double)))
            rmg_error_handler (__FILE__, __LINE__,"error reading");


        close(fhand);

        for (int ip = 0; ip < num_prj_perion[ion]; ip++)
        {
            proj_local_index = this->index_global_to_proj[proj_count];

            double *beta_ip = &beta[ip * max_nlpoints];
            for(int ix = 0; ix < this->dimx[proj_count]; ix++)
                for(int iy = 0; iy < this->dimy[proj_count]; iy++)
                    for(int iz = 0; iz < this->dimz[proj_count]; iz++)
                    {
                        int ixx = ix + this->ixmin[proj_count]; 
                        int iyy = iy + this->iymin[proj_count]; 
                        int izz = iz + this->izmin[proj_count]; 

                        if (ixx < 0) ixx += NX_GRID;
                        if (iyy < 0) iyy += NY_GRID;
                        if (izz < 0) izz += NZ_GRID;
                        if (ixx >= NX_GRID) ixx -=NX_GRID;
                        if (iyy >= NY_GRID) iyy -=NY_GRID;
                        if (izz >= NZ_GRID) izz -=NZ_GRID;

                        if(ixx >=ilow && ixx < ihigh && iyy >=jlow && iyy <jhigh && izz >=klow && izz < khigh)
                        {
                            int idx = (ixx - ilow) * PY0_GRID *PZ0_GRID + (iyy-jlow)*PZ0_GRID + izz-klow;
                            int idx0 = ix * this->dimy[proj_count] * this->dimz[proj_count] + iy * this->dimz[proj_count] + iz;
                            this->storage_cpu[proj_local_index * P0_BASIS + idx] = (KpointType) beta_ip[idx0];
                        }

                    }
            proj_count++;
        }

        delete [] beta;
    }

}

template void LocalObject<double>::GetAtomicOrbitals(int, BaseGrid &BG);
template void LocalObject<std::complex<double>>::GetAtomicOrbitals(int, BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::GetAtomicOrbitals(int num_ions, BaseGrid &BG)
{

    // Used to generate LDA+U orbital projectors that span the full space
    if (!this->delocalized)
        rmg_error_handler (__FILE__, __LINE__,"error GetAtomicOrbitals");

    KpointType *weight;
    std::complex<double> I_t(0.0, 1.0);
    int pbasis = BG.get_P0_BASIS(density);

    /*Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;


    /*Get memory to store the phase factor applied to the forward Fourier transform 
     * and to store the backwards transform*/
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);

    if ((beptr == NULL) || (gbptr == NULL))
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];


    double kvec[3];

    kvec[0] = 0.0;
    kvec[1] = 0.0;
    kvec[2] = 0.0;

    /* Loop over ions */
    weight = (KpointType *)this->storage_cpu;
    for (int ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        ION *iptr = &Atoms[ion];

        /* Get species type */
        SPECIES *sp = &Species[iptr->species];

        int nlxdim = get_NX_GRID();
        int nlydim = get_NY_GRID();
        int nlzdim = get_NZ_GRID();

        double vect[3], nlcrds[3];

        /* Find nlcdrs, vector that gives shift of ion from center of its ionic box */
        /* for delocalized case it's just half the cell dimensions */
        vect[0] = iptr->xtal[0] - 0.5;
        vect[1] = iptr->xtal[1] - 0.5;
        vect[2] = iptr->xtal[2] - 0.5;

        /*The vector we are looking for should be */
        to_cartesian (vect, nlcrds);

        /*Calculate the phase factor */
        FindPhaseKpoint (kvec, nlxdim, nlydim, nlzdim, nlcrds, fftw_phase, false);


        /*Temporary pointer to the already calculated forward transform */
        int kpt = 0;
        fptr = (std::complex<double> *)&sp->forward_orbital[kpt * sp->num_orbitals * pbasis];


        /* Loop over radial projectors */
        for (int ip = 0; ip < sp->num_orbitals; ip++)
        {

            // This ranges over all orbitals including the m-dependence
            if(sp->awave_is_ldaU[ip])
            {


                // Apply the phase factor.
                for (int idx = 0; idx < pbasis; idx++) gbptr[idx] =  fptr[idx] * std::conj(fftw_phase[idx]);

                /*Do the backwards transform */
                coarse_pwaves->FftInverse(gbptr, beptr);

                if(ct.is_gamma)
                {
                    double *weight_R = (double *)weight;
                    for (int idx = 0; idx < pbasis; idx++) weight_R[idx] = std::real(beptr[idx]);
                }
                else
                {
                    std::complex<double> *weight_C = (std::complex<double> *)weight;
                    for (int idx = 0; idx < pbasis; idx++) weight_C[idx] = beptr[idx];
                }

                weight += pbasis;
            }

            /*Advance the temp pointers */
            fptr += pbasis;

        }


    }                           /* end for */


    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);

}

//localized orbitals are projected on 3D domain decompostion, writing for each orbital here is not very 
//efficient  so we write the orbitals for each processor 
template void LocalObject<double>::WriteOrbitals(std::string filename, BaseGrid &BG);
template void LocalObject<std::complex<double>>::WriteOrbitals(std::string filename, BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::WriteOrbitals(std::string filename, BaseGrid &BG)
{

    int density = this->density;
    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;

    int fhand;

    std::string newname= filename + "_spin" + std::to_string(pct.spinpe); 
    newname =  newname + "_gridpe" + std::to_string(pct.gridpe);
    fhand = open(newname.c_str(), O_CREAT |O_TRUNC| O_RDWR, S_IREAD | S_IWRITE);
    if(fhand < 0)
    {
        printf ("\n cannot open file: %s \n", newname.c_str());
        fflush(NULL);
        exit(0);
    }

    size_t size = this->num_thispe * P0_BASIS * sizeof(KpointType);
    write(fhand, this->storage_cpu, size);

    close(fhand);
}

template void LocalObject<double>::SetZeroBoundary(BaseGrid&, int kh_level, int fd_order);
template void LocalObject<std::complex<double>>::SetZeroBoundary(BaseGrid&, int kh_level, int fd_order);
template <class KpointType> void LocalObject<KpointType>::SetZeroBoundary(BaseGrid &BG, int kh_level, int fd_order)
{

    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    int PX_OFFSET = BG.get_PX_OFFSET(density);
    int PY_OFFSET = BG.get_PY_OFFSET(density);
    int PZ_OFFSET = BG.get_PZ_OFFSET(density);
    int NX_GRID = BG.get_NX_GRID(density);
    int NY_GRID = BG.get_NY_GRID(density);
    int NZ_GRID = BG.get_NZ_GRID(density);
    int ilow = PX_OFFSET;
    int ihigh = ilow + PX0_GRID;
    int jlow = PY_OFFSET;
    int jhigh = jlow + PY0_GRID;
    int klow = PZ_OFFSET;
    int khigh = klow + PZ0_GRID;

    //int dim = 1 << kh_level;
    //if( (PX0_GRID % dim != 0) || (PY0_GRID % dim != 0) || (PZ0_GRID % dim != 0) )
    // {
    //    std::cout << " multigrid level of " << kh_level << " require grid on each processor divisible by " << dim << std::endl;
    //}

    this->mask = new char[this->num_thispe * P0_BASIS];
    for(int i = 0; i < this->num_thispe * P0_BASIS; i++) mask[i] = 0;
    char *mask_x = new char[PX0_GRID];
    char *mask_y = new char[PY0_GRID];
    char *mask_z = new char[PZ0_GRID];

    for (int i_local = 0; i_local < this->num_thispe; i_local++)
    {
        for(int ix = 0; ix < PX0_GRID; ix++) mask_x[ix] = 0;
        for(int iy = 0; iy < PY0_GRID; iy++) mask_y[iy] = 0;
        for(int iz = 0; iz < PZ0_GRID; iz++) mask_z[iz] = 0;

        int i = this->index_proj_to_global[i_local];
        if(i < 0) continue;

        for(int ix = ixmin[i] + fd_order/2; ix < ixmin[i] + dimx[i] - fd_order/2; ix++) 
        {
            int ixx = (ix+NX_GRID) % NX_GRID;
            if( ixx >= ilow && ixx < ihigh)  mask_x[ixx-ilow] = 1;
        }
        for(int iy = iymin[i] + fd_order/2; iy < iymin[i] + dimy[i] - fd_order/2; iy++) 
        {
            int iyy = (iy+NY_GRID) % NY_GRID;
            if( iyy >= jlow && iyy < jhigh)  mask_y[iyy-jlow] = 1;
        }
        for(int iz = izmin[i] + fd_order/2; iz < izmin[i] + dimz[i] - fd_order/2; iz++) 
        {
            int izz = (iz+NZ_GRID) % NZ_GRID;
            if( izz >= klow && izz < khigh)  mask_z[izz-klow] = 1;
        }

        for(int ix = 0; ix < PX0_GRID; ix++)
            for(int iy = 0; iy < PY0_GRID; iy++)
                for(int iz = 0; iz < PZ0_GRID; iz++)
                {
                    int idx = ix * PY0_GRID * PZ0_GRID + iy * PZ0_GRID + iz;
                    mask[i_local * P0_BASIS + idx] = mask_x[ix] * mask_y[iy] * mask_z[iz];
                }
    }

    delete [] mask_x;
    delete [] mask_y;
    delete [] mask_z;


}

template void LocalObject<double>::ReAssign(BaseGrid &BG);
template void LocalObject<std::complex<double>>::ReAssign(BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::ReAssign(BaseGrid &BG)
{

    // reassign orrbitals so that one orbital's local index on different proceessors are the same.
    // we can perform multigrid preconditioning with distributed domain decomposioitn. 
    if(this->delocalized) return;

    int rank,npes;

    MPI_Comm_size(this->comm, &npes);
    MPI_Comm_rank(this->comm, &rank);
    int *orbital_proj = new int[npes * this->num_tot];
    char *onerow = new char[npes];
    int *assigned = new int[this->num_tot];

    for(int idx = 0; idx < npes * this->num_tot; idx++) orbital_proj[idx] = 0;

    for(int i = 0; i < this->num_thispe; i++)
    {
        int i_glob = this->index_proj_to_global[i];
        orbital_proj[i_glob * npes + rank] = 1;
    }

    
    
    int size_mat = npes * this->num_tot;
    MPI_Allreduce(MPI_IN_PLACE, orbital_proj, size_mat, MPI_INT, MPI_SUM, this->comm);


    int num_proj = 0;
    bool can_be_assigned, no_left;
    for(int i = 0; i < this->num_tot; i++) assigned[i] = 0;

    for(int i = 0; i < this->num_tot; i++) 
    {
        this->index_global_to_proj[i] = -1;
        this->index_proj_to_global[i] = -1;
    }
    for(int i = 0; i < this->num_tot; i++)   // maximu number of projected orbitals will be the num_tot
    {
        for (int ip = 0; ip < npes; ip++) onerow[ip] = 0;

        no_left = true;
        for(int j = i; j < this->num_tot; j++)  
        {
            if( assigned[j] ) continue;
            no_left = false;
            // check if this orbital can assign to this processor's nmu_proj orbital
            can_be_assigned = true;
            for (int ip = 0; ip < npes; ip++) 
            {
                if( onerow[ip] + orbital_proj[j* npes + ip] >1) 
                {
                    can_be_assigned = false;
                    break; 
                }
            }

            if( can_be_assigned ) 
            {

                for (int ip = 0; ip < npes; ip++) 
                {
                    if( orbital_proj[j * npes + ip] ) 
                    {
                        onerow[ip] = 1;
                    }
                }

                if(orbital_proj[j * npes + rank])
                {
                    this->index_proj_to_global[i] = j;
                    this->index_global_to_proj[j] = i;
                }

                assigned[j] = 1;
            }
        }

        if(no_left) break;

        num_proj++;

    }

    this->num_thispe = num_proj;


    for(int i = 0; i < this->num_thispe; i++)
    {
        //        std::cout << rank << "  " <<i<< "  aaaa  " << this->index_proj_to_global[i] << std::endl;
    }

    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;

    RmgFreeHost(this->storage_cpu);
    gpuFree(this->storage_gpu);
    size_t size = this->num_thispe * P0_BASIS *sizeof(KpointType) + 8;
    this->storage_size = size;
    this->storage_cpu = (KpointType *) RmgMallocHost(size);
    gpuMalloc((void **)&this->storage_gpu, size);
    this->storage_ptr = MemoryPtrHostDevice(this->storage_cpu, this->storage_gpu);

    delete [] orbital_proj;
    delete [] onerow;
    delete [] assigned;

}



template void LocalObject<double>::Normalize();
template void LocalObject<std::complex<double>>::Normalize();
template <class KpointType> void LocalObject<KpointType>::Normalize()
{

    double t1 = Rmg_G->get_NX_GRID(this->density);
    t1 *= Rmg_G->get_NY_GRID(this->density);
    t1 *= Rmg_G->get_NZ_GRID(this->density);

    double vol = Rmg_L.get_omega() /t1;

    int P0_BASIS = Rmg_G->get_P0_BASIS(this->density);
    double *norm_coef = new double[this->num_tot]();
    for(int st = 0; st < this->num_thispe; st++)
    {
        int st_glob = this->index_proj_to_global[st];
        for(int idx = 0; idx < P0_BASIS; idx++)
            norm_coef[st_glob] += std::norm(this->storage_cpu[st*P0_BASIS + idx]);
    }

    MPI_Allreduce(MPI_IN_PLACE, norm_coef, this->num_tot, MPI_DOUBLE, MPI_SUM, this->comm);

    for(int st = 0; st < this->num_thispe; st++)
    {
        int st_glob = this->index_proj_to_global[st];
        double alpha = std::sqrt(norm_coef[st_glob] * vol);
        for(int idx = 0; idx < P0_BASIS; idx++)
            this->storage_cpu[st*P0_BASIS + idx] /= alpha;
    }

    delete [] norm_coef;
}
template void LocalObject<double>::AssignOrbital(int st, double *psi);
template void LocalObject<std::complex<double>>::AssignOrbital(int st, std::complex<double> *psi);
template <class KpointType> void LocalObject<KpointType>::AssignOrbital(int st, KpointType *psi)
{

    int density = this->density;
    int PX0_GRID = Rmg_G->get_PX0_GRID(density);
    int PY0_GRID = Rmg_G->get_PY0_GRID(density);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    int PX_OFFSET = Rmg_G->get_PX_OFFSET(density);
    int PY_OFFSET = Rmg_G->get_PY_OFFSET(density);
    int PZ_OFFSET = Rmg_G->get_PZ_OFFSET(density);
    int NX_GRID = Rmg_G->get_NX_GRID(density);
    int NY_GRID = Rmg_G->get_NY_GRID(density);
    int NZ_GRID = Rmg_G->get_NZ_GRID(density);
    int ilow = PX_OFFSET;
    int ihigh = ilow + PX0_GRID;
    int jlow = PY_OFFSET;
    int jhigh = jlow + PY0_GRID;
    int klow = PZ_OFFSET;
    int khigh = klow + PZ0_GRID;


    for(int idx = 0; idx < P0_BASIS; idx++) this->storage_cpu[st * P0_BASIS + idx] = 0.0;

    int st_glob = this->index_proj_to_global[st];

    for(int ix = 0; ix < this->dimx[st_glob]; ix++)
        for(int iy = 0; iy < this->dimy[st_glob]; iy++)
            for(int iz = 0; iz < this->dimz[st_glob]; iz++)
            {
                int ixx = ix + this->ixmin[st_glob]; 
                int iyy = iy + this->iymin[st_glob]; 
                int izz = iz + this->izmin[st_glob]; 

                if (ixx < 0) ixx += NX_GRID;
                if (iyy < 0) iyy += NY_GRID;
                if (izz < 0) izz += NZ_GRID;
                if (ixx >= NX_GRID) ixx -=NX_GRID;
                if (iyy >= NY_GRID) iyy -=NY_GRID;
                if (izz >= NZ_GRID) izz -=NZ_GRID;

                if(ixx >=ilow && ixx < ihigh && iyy >=jlow && iyy <jhigh && izz >=klow && izz < khigh)
                {
                    int idx = (ixx - ilow) * PY0_GRID *PZ0_GRID + (iyy-jlow)*PZ0_GRID + izz-klow;
                    int idx0 = ix * this->dimy[st_glob] * this->dimz[st_glob] + iy * this->dimz[st_glob] + iz;
                    this->storage_cpu[st * P0_BASIS + idx] = psi[idx0];
                }

            }

}

template void LocalObject<double>::ReadProjectedOrbitals(std::string filename, BaseGrid &BG);
template void LocalObject<std::complex<double>>::ReadProjectedOrbitals(std::string filename, BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::ReadProjectedOrbitals(std::string filename, BaseGrid &BG)
{

    int density = this->density;
    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;

    int fhand;


    std::string newname= filename + "_spin" + std::to_string(pct.spinpe);
    newname =  newname + "_gridpe" + std::to_string(pct.gridpe);
    fhand = open(newname.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
    if(fhand < 0)
    {
        printf ("\n cannot open file: %s \n", newname.c_str());
        fflush(NULL);
        exit(0);
    }

    size_t size = this->num_thispe * P0_BASIS * sizeof(KpointType);
    read(fhand, this->storage_cpu, size);

    close(fhand);
}

template void LocalObject<double>::SetOrbitalProjToSingle(BaseGrid &BG);
template void LocalObject<std::complex<double>>::SetOrbitalProjToSingle(BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::SetOrbitalProjToSingle(BaseGrid &BG)
{
    int rank,npes;

    MPI_Comm_size(this->comm, &npes);
    MPI_Comm_rank(this->comm, &rank);
    int *orbital_proj = new int[npes * this->num_tot];
    for(int idx = 0; idx < npes * this->num_tot; idx++) orbital_proj[idx] = 0;

    for(int i = 0; i < this->num_thispe; i++)
    {
        int i_glob = this->index_proj_to_global[i];
        orbital_proj[i_glob * npes + rank] = 1;
    }



    int size_mat = npes * this->num_tot;
    MPI_Allreduce(MPI_IN_PLACE, orbital_proj, size_mat, MPI_INT, MPI_SUM, this->comm);

    int NX_GRID = BG.get_NX_GRID(this->density);
    int NY_GRID = BG.get_NY_GRID(this->density);
    int NZ_GRID = BG.get_NZ_GRID(this->density);
    int PX_OFFSET, PY_OFFSET, PZ_OFFSET;

    this->OrbitalOwnedProcs_pelist = new std::vector<int>[this->num_tot];
    this->OrbitalOwnedProcs_xoffset = new std::vector<int>[this->num_tot];
    this->OrbitalOwnedProcs_yoffset = new std::vector<int>[this->num_tot];
    this->OrbitalOwnedProcs_zoffset = new std::vector<int>[this->num_tot];
    for (int st = 0; st < this->num_tot; st++)
    {
        for(int pe = 0; pe < npes; pe++)
        {
            if(orbital_proj[st * npes + pe] == 1)
            {
                BG.find_node_offsets(pe, NX_GRID, NY_GRID, NZ_GRID,
                        &PX_OFFSET, &PY_OFFSET, &PZ_OFFSET);
                int xoff = this->ixmin[st] - PX_OFFSET;
                int yoff = this->iymin[st] - PY_OFFSET;
                int zoff = this->izmin[st] - PZ_OFFSET;

                if( xoff >= NX_GRID/2) xoff -= NX_GRID;
                if( yoff >= NY_GRID/2) yoff -= NY_GRID;
                if( zoff >= NZ_GRID/2) zoff -= NZ_GRID;

                if( xoff <= -NX_GRID/2) xoff += NX_GRID;
                if( yoff <= -NY_GRID/2) yoff += NY_GRID;
                if( zoff <= -NZ_GRID/2) zoff += NZ_GRID;

                this->OrbitalOwnedProcs_pelist[st].push_back(pe);
                this->OrbitalOwnedProcs_xoffset[st].push_back(xoff);
                this->OrbitalOwnedProcs_yoffset[st].push_back(yoff);
                this->OrbitalOwnedProcs_zoffset[st].push_back(zoff);
            }
        }
        if(pct.gridpe == 0 && 0)
        {
            printf("\n orbital %d: proc  xoffset  yoffset  zoffset", st);
            for(size_t i = 0; i < this->OrbitalOwnedProcs_pelist[st].size(); i++)
                printf("\n             %d    %d   %d   %d", 
                        this->OrbitalOwnedProcs_pelist[st][i], 
                        this->OrbitalOwnedProcs_xoffset[st][i], 
                        this->OrbitalOwnedProcs_yoffset[st][i], 
                        this->OrbitalOwnedProcs_zoffset[st][i]); 
        }
    }
    delete [] orbital_proj;
}


template void LocalObject<double>::WriteOrbitalsToSingleFiles(std::string filename, BaseGrid &BG);
template void LocalObject<std::complex<double>>::WriteOrbitalsToSingleFiles(std::string filename, BaseGrid &BG);
template <class KpointType> void LocalObject<KpointType>::WriteOrbitalsToSingleFiles(std::string filename, BaseGrid &BG)
{

    int density = this->density;
    int PX0_GRID = BG.get_PX0_GRID(density);
    int PY0_GRID = BG.get_PY0_GRID(density);
    int PZ0_GRID = BG.get_PZ0_GRID(density);
    int P0_BASIS = PX0_GRID * PY0_GRID * PZ0_GRID;
    int NX_GRID = BG.get_NX_GRID(this->density);
    int NY_GRID = BG.get_NY_GRID(this->density);
    int NZ_GRID = BG.get_NZ_GRID(this->density);

    // max number of procs in one direction
    int xb = (this->max_dimx + PX0_GRID -1)/PX0_GRID +1;
    int yb = (this->max_dimy + PY0_GRID -1)/PY0_GRID +1;
    int zb = (this->max_dimz + PZ0_GRID -1)/PZ0_GRID +1;
    xb = std::min(xb, BG.get_PE_X());
    yb = std::min(yb, BG.get_PE_Y());
    zb = std::min(zb, BG.get_PE_Z());


    MPI_Request *request = new MPI_Request[xb*yb*zb];
    KpointType *phi = new KpointType[this->max_dimx * this->max_dimy * this->max_dimz];

    KpointType *phi_proj = new KpointType[xb*yb*zb*P0_BASIS];

    int fhand;
    int factor = sizeof(KpointType)/sizeof(double);

    for(int st_glob = 0; st_glob < this->num_tot; st_glob++)
    {

        int st = this->index_global_to_proj[st_glob];

        int num_pelist = this->OrbitalOwnedProcs_pelist[st_glob].size();
        if(num_pelist > xb * yb * zb)
        {
            std::cout << "something wrong with pelist "<< std::endl;
            exit(0);
        }
        int recv_pe = st_glob%pct.grid_npes;
        if(recv_pe== pct.gridpe)
        {

            for(int irecv = 0; irecv < num_pelist; irecv++)
            {

                int pe =  this->OrbitalOwnedProcs_pelist[st_glob][irecv];
                int tag = recv_pe * pct.grid_npes + pe;

                MPI_Irecv(&phi_proj[irecv * P0_BASIS], factor *P0_BASIS, MPI_DOUBLE, pe, tag, this->comm, &request[irecv]); 
            }
        }

        for(int irecv = 0; irecv < num_pelist; irecv++)
        {

            int pe =  this->OrbitalOwnedProcs_pelist[st_glob][irecv];
            if(pe == pct.gridpe )
            {

                KpointType *sendbuf = &this->storage_cpu[st * P0_BASIS];
                int tag = recv_pe * pct.grid_npes + pe;
                MPI_Send(sendbuf, factor * P0_BASIS, MPI_DOUBLE, recv_pe, tag, this->comm);
            }

        }

        if(recv_pe== pct.gridpe)
        {
            for(int idx = 0; idx < this->dimx[st_glob] * this->dimy[st_glob] * this->dimz[st_glob]; idx++)
                phi[idx] = 0.0;
            for(int i = 0; i < num_pelist; i++)
            {
                int irecv;
                MPI_Waitany(num_pelist, request, &irecv, MPI_STATUS_IGNORE);
                int xoff = this->OrbitalOwnedProcs_xoffset[st_glob][irecv];
                int yoff = this->OrbitalOwnedProcs_yoffset[st_glob][irecv];
                int zoff = this->OrbitalOwnedProcs_zoffset[st_glob][irecv];
                for(int ix = 0; ix < PX0_GRID; ix++)
                {
                    int ixx = ix - xoff;
                    if( ixx < 0 || ixx >= this->dimx[st_glob])
                    {
                        ixx = (ixx + NX_GRID) % NX_GRID;
                        if( ixx < 0 || ixx >= this->dimx[st_glob]) continue;
                    }

                    for(int iy = 0; iy < PY0_GRID; iy++)
                    {
                        int iyy = iy - yoff;
                        if(iyy < 0 || iyy >= this->dimy[st_glob]) 
                        {
                            iyy = (iyy + NY_GRID) % NY_GRID;
                            if( iyy < 0 || iyy >= this->dimy[st_glob]) continue;
                        }

                        for(int iz = 0; iz < PZ0_GRID; iz++)
                        {
                            int izz = iz - zoff;
                            if(izz < 0 || izz >= this->dimz[st_glob])
                            {
                                izz = (izz + NZ_GRID) % NZ_GRID;
                                if( izz < 0 || izz >= this->dimz[st_glob])continue;
                            }

                            int idx1 = ix * PY0_GRID * PZ0_GRID;
                            int idx2 = ixx * this->dimy[st_glob] * this->dimz[st_glob];
                            idx1 += iy * PZ0_GRID;
                            idx2 += iyy * this->dimz[st_glob];

                            phi[idx2 + izz] = phi_proj[irecv * P0_BASIS + idx1 + iz];
                        } 
                    }
                }

            }

            std::string newname= filename + "_spin"+std::to_string(pct.spinpe)+".orbit_"+std::to_string(st_glob);
            fhand = open(newname.c_str(),  O_CREAT | O_TRUNC |O_RDWR, S_IREAD | S_IWRITE);
            if(fhand < 0)
            {
                printf ("\n cannot open file: %s \n", newname.c_str());
                fflush(NULL);
                exit(0);
            }

            size_t size = this->dimx[st_glob] * this->dimy[st_glob] * this->dimz[st_glob] *sizeof(KpointType);
            write(fhand, phi, size);

            int ix0 = this->ixmin[st_glob];
            int ix1 = this->ixmin[st_glob] + this->dimx[st_glob];
            int iy0 = this->iymin[st_glob];
            int iy1 = this->iymin[st_glob] + this->dimy[st_glob];
            int iz0 = this->izmin[st_glob];
            int iz1 = this->izmin[st_glob] + this->dimz[st_glob];
            write(fhand, &ix0, sizeof(int));
            write(fhand, &ix1, sizeof(int));
            write(fhand, &iy0, sizeof(int));
            write(fhand, &iy1, sizeof(int));
            write(fhand, &iz0, sizeof(int));
            write(fhand, &iz1, sizeof(int));

            close(fhand);
        }
    }
    delete [] request;
    delete [] phi_proj;
    delete [] phi;
}


template LocalObject<double>::LocalObject(const LocalObject<double> &);
template LocalObject<std::complex<double>>::LocalObject(const LocalObject<std::complex<double>> &);
template <class KpointType> LocalObject<KpointType>::LocalObject(const LocalObject<KpointType> &Old_LO)
{
    int num_tot = Old_LO.num_tot;
    this->index_global_to_proj = new int[num_tot];
    this->index_proj_to_global = new int[num_tot];
    for(int i = 0; i < num_tot; i++)
    {
        this->index_global_to_proj[i]  = Old_LO.index_global_to_proj[i];
        this->index_proj_to_global[i]  = Old_LO.index_proj_to_global[i];

    }
    size_t size = Old_LO.num_thispe * Old_LO.pbasis *sizeof(KpointType) +8;
    this->storage_size = size;
    this->storage_cpu = (KpointType *) RmgMallocHost(size);
    gpuMalloc( (void **)&this->storage_gpu, size);
    this->storage_ptr = MemoryPtrHostDevice(this->storage_cpu, this->storage_gpu);
}
