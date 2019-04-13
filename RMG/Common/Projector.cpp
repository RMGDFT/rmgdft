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
#include <stdlib.h>
#include <string.h>
#include <complex>
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
//#include "rmg_error.h"
#include "Kpoint.h"
#include "RmgSumAll.h"
#include "Projector.h"



template Projector<double>::Projector(Kpoint<double> *, int);
template Projector<std::complex<double>>::Projector(Kpoint<std::complex<double>> *, int);
template Projector<double>::~Projector(void);
template Projector<std::complex<double>>::~Projector(void);
template void Projector<double>::set_storage(double *);
template void Projector<std::complex<double>>::set_storage(std::complex<double> *);


template <class KpointType> Projector<KpointType>::Projector(Kpoint<KpointType> *K, int projector_type)
{
    this->type = projector_type;
    this->kptr = K;
    this->np = 0;
    this->num_nonloc_ions = 0;
    this->num_owned_ions = 0;
    this->num_owned_pe = 0;
    this->num_owners = 0;
    this->owned_pe_list = new int[pct.grid_npes]();
    this->num_owned_ions_per_pe = new int[pct.grid_npes]();
    this->owners_list = new int[pct.grid_npes]();
    this->num_nonowned_ions_per_pe = new int[ct.num_ions]();
    this->nonloc_ions_list = new int[ct.num_ions]();
    this->nonloc_ion_ownflag = new bool[ct.num_ions]();

    this->idxptrlen = new int[ct.num_ions](); 

    size_t NX_GRID = get_NX_GRID();
    size_t NY_GRID = get_NY_GRID();
    size_t NZ_GRID = get_NZ_GRID();

    size_t FNX_GRID = get_FNX_GRID();
    size_t FNY_GRID = get_FNY_GRID();
    size_t FNZ_GRID = get_FNZ_GRID();

    size_t PX0_GRID = get_PX0_GRID();
    size_t PY0_GRID = get_PY0_GRID();
    size_t PZ0_GRID = get_PZ0_GRID();

    size_t PX_OFFSET = get_PX_OFFSET();
    size_t PY_OFFSET = get_PY_OFFSET();
    size_t PZ_OFFSET = get_PZ_OFFSET();

    // The serial version of this is a little faster for small systems
    // but is horrendously slow when there are thousands of ions. We
    // parallelize over ions but keep the number of threads used limited
    // to a maximum of 4 due to the large amounts of memory required
    // per thread.
    int nthreads = std::min(ct.OMP_THREADS_PER_NODE, 4);
#pragma omp parallel
    {

        int *Aix = new int[NX_GRID];
        int *Aiy = new int[NY_GRID];
        int *Aiz = new int[NZ_GRID];

        #pragma omp parallel for schedule(static, 1) num_threads(nthreads)
        for (int ion = 0; ion < ct.num_ions; ion++)
        /* Loop over ions */
        {
            int ilow, jlow, klow, ihi, jhi, khi, map;
            double vect[3];

            /* Generate ion pointer */
            ION *iptr = &ct.ions[ion];

            /* Get species type */
            SPECIES *sp = &ct.sp[iptr->species];

            int icenter = sp->nldim / 2;
            int icut = (icenter + 1) * (icenter + 1);

            /* Determine mapping indices or even if a mapping exists */
            int nlxdim = sp->nldim;
            int nlydim = sp->nldim;
            int nlzdim = sp->nldim;
            if(this->type == DELOCALIZED) {
                map = true;
                nlxdim = NX_GRID;
                nlydim = NY_GRID;
                nlzdim = NZ_GRID;
                xcstart = ycstart = zcstart = 0.0;
                ilow = PX_OFFSET;
                ihi = ilow + get_PX0_GRID() -1;
                jlow = PY_OFFSET;
                jhi = jlow + get_PY0_GRID() -1;
                klow = PZ_OFFSET;
                khi = klow + get_PZ0_GRID() -1;
            }
            else
            {
                map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                             sp->nldim, PX0_GRID, PY0_GRID, PZ0_GRID,
                             NX_GRID, NY_GRID, NZ_GRID,
                             &this->xcstart, &this->ycstart, &this->zcstart);
            }

            /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
            /*xtal vector between ion and left bottom corner of the box */
            vect[0] = iptr->xtal[0] - this->xcstart;
            vect[1] = iptr->xtal[1] - this->ycstart;
            vect[2] = iptr->xtal[2] - this->zcstart;

            /*Substract vector between left bottom corner of the box and center of the box */
            vect[0] -= (nlxdim / 2) / (double) NX_GRID;
            vect[1] -= (nlydim / 2) / (double) NY_GRID;
            vect[2] -= (nlzdim / 2) / (double) NZ_GRID;

            /*The vector we are looking for should be */
            to_cartesian (vect, iptr->nlcrds);

            /* If there is a mapping for this ion then we have to generate */
            /* the projector.                                              */
            int icount = 0;
            if (map)
            {

                /* Generate index arrays */
                int idx = 0;
                for (int ix = 0; ix < nlxdim; ix++)
                {
                    for (int iy = 0; iy < nlydim; iy++)
                    {
                        for (int iz = 0; iz < nlzdim; iz++)
                        {
                            if(this->type == LOCALIZED)
                            {
                                if ((((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                                     ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                                     ((Aiz[iz] >= klow) && (Aiz[iz] <= khi))))
                                {
                                    /* Cut it off if required */
                                    int itmp =
                                        (ix - icenter) * (ix - icenter) +
                                        (iy - icenter) * (iy - icenter) + (iz - icenter) * (iz - icenter);

                                    if (icut >= itmp) icount++;
                                }
                            }
                            else 
                            {
                                if ((((ix >= ilow) && (ix <= ihi)) &&
                                     ((iy >= jlow) && (iy <= jhi)) &&
                                     ((iz >= klow) && (iz <= khi)))) {

                                        icount++;
                                } 
                            }
                            idx++;
                        }
                    }
                }

                /* Save number of points */
                this->idxptrlen[ion] = icount;

            }                       /* end if (map) */

        }                           /* end for (ion = 0; ion < ct.num_ions; ion++) */

        
        delete [] Aiz;
        delete [] Aiy;
        delete [] Aix;

    } // end omp private area


    for (int ion = 0; ion < ct.num_ions; ion++)
    {

        /* Generate ion pointer */
        ION *iptr = &ct.ions[ion];

        /*Add ion into list of nonlocal ions if it has overlap with given processor */
        if (this->idxptrlen[ion])
        {
            this->nonloc_ions_list[this->num_nonloc_ions] = ion;

            /*Ownership flag for current ion */
            this->nonloc_ion_ownflag[this->num_nonloc_ions] = false;

            /*See if this processor owns current ion */
            if (pct.gridpe == claim_ion (iptr->xtal, PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID))
            {
                this->owned_ions_list[this->num_owned_ions] = ion;
                this->num_owned_ions++;
                this->nonloc_ion_ownflag[this->num_nonloc_ions] = true;
            }
            this->num_nonloc_ions++;
        }

    } // end for (ion = 0; ion < ct.num_ions; ion++)


    int *Aix = new int[NX_GRID]();
    int *Aiy = new int[NY_GRID]();
    int *Aiz = new int[NZ_GRID]();

    int *Aix2 = new int[FNX_GRID]();
    int *Aiy2 = new int[FNY_GRID]();
    int *Aiz2 = new int[FNZ_GRID]();

    this->np = this->num_nonloc_ions * ct.max_nl;

    int known_nonowner, nonwner_index, known_owner, owner_index, owned_ions_per_pe, nonowned_ions_per_pe; 
    int owned, owner;

    /*Make sure that ownership of ions is properly established
     * This conditional can be removed if it is found that claim_ions works reliably*/
    if (int_sum_all (this->num_owned_ions, pct.grid_comm) != ct.num_ions)
        rmg_error_handler (__FILE__, __LINE__, "Problem with claimimg ions.");
  
    
    /* Loop over all ions to obtain the lists necessary for communication */
    for (int nlion = 0; nlion < this->num_nonloc_ions; nlion++)
    {
        int ion = this->nonloc_ions_list[nlion];
        int map, map2;

        /* Generate ion pointer */
        ION *iptr = &ct.ions[ion];

        /* Get species type */
        SPECIES *sp = &ct.sp[iptr->species];

        /*See if this processor owns current ion */
        owner = claim_ion (iptr->xtal, PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID);
	owned = false;
        if (pct.gridpe == owner) owned = true;

	/*Case when current PE owns an ion, in such a case we need to assemble
	 * list of other PEs that have overlap - those non-owning PEs will be sending data*/
	if (owned)
	{
	    /*Loop over all processors in caculation to determine if it overlaps with current ion */
	    /* This should only loop over processors in a grid */
	    for (int pe = 0; pe < pct.grid_npes; pe++)
	    {
		/*Skip the case when pe is equal to the rank of current processor, in that case no send or
		 * receive is necessary*/
		if (pe == pct.gridpe)
		    continue;

		/* Determine if ion has overlap with a given PE becasue of beta functions */
                if(!ct.localize_projectors) map = true;
                else
                    map = test_overlap (pe, iptr, Aix, Aiy, Aiz, sp->nldim, PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID);

                /* Determine if ion has overlap with a given PE becasue of Q function */
                if(ct.norm_conserving_pp) {
                    map2 = FALSE;
                }
                else {
                    map2 = test_overlap (pe, iptr, Aix2, Aiy2, Aiz2, sp->qdim, 
                            get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), FNX_GRID, FNY_GRID, FNZ_GRID);
                }


                if (map || map2)
                {
                    /* See if this non-owning pe has already been recorded*/
                    known_nonowner = 0;
                    for (int i = 0; i < this->num_owned_pe; i++)
                    {
                        if (this->owned_pe_list[i] == pe)
                        {
                            known_nonowner++;
                            nonwner_index = i;
                            break;
                        }
                    }

                    /* Record pe index into communication list, if it has not been recorded yet*/
                    if (!known_nonowner)
                    {
                        this->owned_pe_list[this->num_owned_pe] = pe;
                        this->num_owned_ions_per_pe[this->num_owned_pe] = 1;
                        this->list_owned_ions_per_pe[this->num_owned_pe][0] = nlion;
                        this->num_owned_pe++;
                    }
                    else
                    {
                        /*Number of ions to communicate with current non-owner PE*/
                        owned_ions_per_pe = this->num_owned_ions_per_pe[nonwner_index];
                        this->list_owned_ions_per_pe[nonwner_index][owned_ions_per_pe] = nlion;
                        this->num_owned_ions_per_pe[nonwner_index]++;
                    }
                }
            }  // End of for (pe=0; pe < pct.grid_npes; pe++)
        }            
        /*For non-owned ion we need to record rank of owner */
        else	
        {
            /*Loop over list of processors that we already have to see if owning "pe" has alredy been recorded */
            known_owner = 0;
            for (int i = 0; i < this->num_owners; i++)
            {
                if (this->owners_list[i] == owner)
                {
                    known_owner++;
                    owner_index = i;
                    break;
                }
            }

            if (!known_owner)
            {
                this->owners_list[this->num_owners] = owner;
                this->num_nonowned_ions_per_pe[this->num_owners] = 1;
                this->list_ions_per_owner[this->num_owners][0] = nlion;
                this->num_owners++;
            }
            else
            {
                nonowned_ions_per_pe = this->num_nonowned_ions_per_pe[owner_index];
                this->list_ions_per_owner[owner_index][nonowned_ions_per_pe] = nlion;
                this->num_nonowned_ions_per_pe[owner_index]++;
            }
        }
    } // for (nlion = 0; nlion < pct.num_nonloc_ions; nlion++)


    /* Release temporary memory */
    delete [] Aiz2;
    delete [] Aiy2;
    delete [] Aix2;
    delete [] Aiz;
    delete [] Aiy;
    delete [] Aix;



}

template <class KpointType> void Projector<KpointType>::set_storage(KpointType *storage)
{
    this->weight = storage;
}

// Applies projectors to orbitals associated with kpoint kptr
template <class KpointType> void Projector<KpointType>::project(KpointType *p, int offset, int n)
{

    // Do delocalized case first
    if((this->type == DELOCALIZED) || (pct.grid_npes == 1))
    {
        int factor = 1;
        KpointType rzero(0.0);
        KpointType alpha(get_vel());

        if(typeid(KpointType) == typeid(std::complex<double>)) factor = 2;
        char *transt = "t", *transn = "n", *transc = "c";
        char *transa;
        transa = transc;
        if(typeid(KpointType) == typeid(double)) transa = transt;

//        int length = factor * ct.num_ions * this->kptr->nstates * ct.max_nl;
// check this
        int length = factor * this->np * this->kptr->nstates;
        RmgGemm (transa, transn, this->np, this->kptr->nstates, this->kptr->pbasis, alpha,
            this->weight, this->kptr->pbasis, &kptr->orbital_storage[offset*this->kptr->pbasis], this->kptr->pbasis,
            rzero, p, this->np);

        if(pct.grid_npes != 1)
            MPI_Allreduce(MPI_IN_PLACE, (double *)p, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        return;
    }

    // And now localized

}



// Destructor
template <class KpointType> Projector<KpointType>::~Projector(void)
{
    delete [] this->idxptrlen;
    delete [] this->nonloc_ion_ownflag;
    delete [] this->nonloc_ions_list;
    delete [] this->num_nonowned_ions_per_pe;
    delete [] this->owners_list;
    delete [] this->num_owned_ions_per_pe;
    delete [] this->owned_pe_list;
}
