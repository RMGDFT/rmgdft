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

#include "LocalObject.h"


template LocalObject<double>::~LocalObject(void);
template LocalObject<std::complex<double>>::~LocalObject(void);

template LocalObject<double>::LocalObject(int, int*, int*, int*, int*, int*, int*, BaseGrid*, MPI_Comm);
template LocalObject<std::complex<double>>::LocalObject(int, int*, int*, int*, int*, int*, int*, BaseGrid*, MPI_Comm);
template <class KpointType> LocalObject<KpointType>::LocalObject(int num_objects, 
        int *ixmin, int *iymin, int *izmin, int *dimx, int *dimy, int *dimz,
        BaseGrid *Rmg_G, MPI_Comm comm)
{
    this->num_tot = num_objects;
    this->num_thispe = 0;

    this->index_global_to_proj = new int[num_objects];
    this->index_proj_to_global = new int[num_objects];
    this->ixmin = new int[num_objects];
    this->iymin = new int[num_objects];
    this->izmin = new int[num_objects];
    this->dimx = new int[num_objects];
    this->dimy = new int[num_objects];
    this->dimz = new int[num_objects];

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

    int density = 1;
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

    this->storage_proj = new KpointType[this->num_thispe * P0_BASIS];
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
    delete [] this->storage_proj;

}

template void LocalObject<double>::ReadOrbitals(std::string filename, BaseGrid *Rmg_G);
template void LocalObject<std::complex<double>>::ReadOrbitals(std::string filename, BaseGrid *Rmg_G);
template <class KpointType> void LocalObject<KpointType>::ReadOrbitals(std::string filename, BaseGrid *Rmg_G)
{
    
    int density = 1;
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

    int fhand;
    for(int st = 0; st > this->num_thispe; st++)
    {

        int st_glob = this->index_proj_to_global[st];
        size_t size = this->dimx[st] * this->dimy[st] * this->dimz[st] *sizeof(double);

        std::string newname= filename + ".orbit_"+std::to_string(st_glob);
        fhand = open(newname.c_str(), O_RDWR, S_IREAD | S_IWRITE);
        if(fhand < 0)
        {
            printf ("\n %s \n", newname.c_str());
            fflush(NULL);
            exit(0);
        }

        double *psi = new double[size/sizeof(double)];
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
                this->storage_proj[st * P0_BASIS + idx] = (KpointType) psi[idx0];
            }
    
        }
        
        delete []psi;
    }
}

