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
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "Kpoint.h"
#include "RmgSumAll.h"
#include "RmgGemm.h"
#include "GpuAlloc.h"
#include "Projector.h"



template Projector<double>::Projector(int, int, int, int, int);
template Projector<std::complex<double>>::Projector(int, int, int, int, int);
template Projector<double>::~Projector(void);
template Projector<std::complex<double>>::~Projector(void);

template void Projector<double>::project(Kpoint<double> *, double *, int, int, double *);
template void Projector<std::complex<double>>::project(Kpoint<std::complex<double>> *, std::complex<double> *, int, int, std::complex<double> *);

template int Projector<double>::get_num_nonloc_ions(void);
template int Projector<double>::get_num_owned_ions(void);
template int * Projector<double>::get_owned_ions_list(void);
template int Projector<double>::get_num_tot_proj(void);
template int * Projector<double>::get_nonloc_ions_list(void);
template int Projector<double>::get_nldim(int);
template int Projector<double>::get_pstride(void);
template int Projector<std::complex<double>>::get_num_nonloc_ions(void);
template int Projector<std::complex<double>>::get_num_owned_ions(void);
template int * Projector<std::complex<double>>::get_owned_ions_list(void);
template int Projector<std::complex<double>>::get_num_tot_proj(void);
template int * Projector<std::complex<double>>::get_nonloc_ions_list(void);
template int Projector<std::complex<double>>::get_nldim(int);
template int Projector<std::complex<double>>::get_pstride(void);



template <class KpointType> Projector<KpointType>::Projector(int projector_type, int num_pes, int num_ions, int stride, int projector_kind) : 
                        list_owned_ions_per_pe(boost::extents[num_pes][num_ions]),
                        list_ions_per_owner(boost::extents[num_pes][num_ions])

{

    this->type = projector_type;    // LOCALIZED or DELOCALIZED
    this->kind = projector_kind;    // ORBITAL_PROJECTOR OR BETA_PROJECTOR
    this->num_tot_proj = 0;
    this->num_nonloc_ions = 0;
    this->num_owned_ions = 0;
    this->num_owned_pe = 0;
    this->num_owners = 0;
    this->owned_ions_list = new int[ct.num_ions]();
    this->owned_pe_list = new int[pct.grid_npes]();
    this->num_owned_ions_per_pe = new int[pct.grid_npes]();
    this->owners_list = new int[pct.grid_npes]();
    this->num_nonowned_ions_per_pe = new int[ct.num_ions]();
    this->nonloc_ions_list = new int[ct.num_ions]();
    this->idxptrlen = new int[ct.num_ions](); 
    this->nldims = new int[ct.num_species](); 
    this->pstride = stride;
    this->nlcrds.reserve(ct.num_ions);

    // Go through the list of atomic species and set up the nldims array
    for(int isp = 0;isp < ct.num_species;isp++)
    {

        SPECIES *sp = &ct.sp[isp];

        // All the same for beta projectors
        if(this->kind == BETA_PROJECTOR)
        {
            this->nldims[isp] = sp->nldim;
        }
        else if(this->kind == ORBITAL_PROJECTOR)
        {
            this->nldims[isp] = sp->adim_wave;    
        }
    }

    for(int i=0;i < num_pes;i++)
    {
        for(int j=0;j < num_ions;j++)
        {
            list_owned_ions_per_pe[i][j] = 0;
            list_ions_per_owner[i][j] = 0;
        }
    }
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
            SPECIES *sp = &ct.sp[iptr->species];
       
            int icenter = this->nldims[iptr->species] / 2;
            int icut = (icenter + 1) * (icenter + 1);

            /* Determine mapping indices or even if a mapping exists */
            int nlxdim = this->nldims[iptr->species];
            int nlydim = this->nldims[iptr->species];
            int nlzdim = this->nldims[iptr->species];
            double nlxcstart = 0.0;
            double nlycstart = 0.0;
            double nlzcstart = 0.0;

            if(this->type == DELOCALIZED) {
                map = true;
                nlxdim = NX_GRID;
                nlydim = NY_GRID;
                nlzdim = NZ_GRID;
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
                             this->nldims[iptr->species], PX0_GRID, PY0_GRID, PZ0_GRID,
                             NX_GRID, NY_GRID, NZ_GRID,
                             &nlxcstart, &nlycstart, &nlzcstart);
            }

            // Make map false if ionic type is not included in lda+u
            if(this->kind == ORBITAL_PROJECTOR)
            {
                if(sp->num_ldaU_orbitals == 0) map = false;
            }

            /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
            /*xtal vector between ion and left bottom corner of the box */
            vect[0] = iptr->xtal[0] - nlxcstart;
            vect[1] = iptr->xtal[1] - nlycstart;
            vect[2] = iptr->xtal[2] - nlzcstart;

            /*Substract vector between left bottom corner of the box and center of the box */
            vect[0] -= (nlxdim / 2) / (double) NX_GRID;
            vect[1] -= (nlydim / 2) / (double) NY_GRID;
            vect[2] -= (nlzdim / 2) / (double) NZ_GRID;

            /*The vector we are looking for should be */
            to_cartesian (vect, this->nlcrds[ion].data());

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
        SPECIES *sp = &ct.sp[iptr->species];

        /*Add ion into list of nonlocal ions if it has overlap with given processor */
        if (((this->kind == BETA_PROJECTOR) && (this->idxptrlen[ion] || pct.Qidxptrlen[ion])) ||
            ((this->kind == ORBITAL_PROJECTOR) && this->idxptrlen[ion] && (sp->num_ldaU_orbitals > 0)))
        {
            this->nonloc_ions_list[this->num_nonloc_ions] = ion;

            /*See if this processor owns current ion */
            if (pct.gridpe == claim_ion (iptr->xtal, PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID))
            {
                this->owned_ions_list[this->num_owned_ions] = ion;
                this->num_owned_ions++;
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

    this->num_tot_proj = this->num_nonloc_ions * this->pstride;

    int known_nonowner, nonwner_index, known_owner, owner_index, owned_ions_per_pe, nonowned_ions_per_pe; 
    int owned, owner;

    /*Make sure that ownership of ions is properly established
     * This conditional can be removed if it is found that claim_ions works reliably*/
    if (int_sum_all (this->num_owned_ions, pct.grid_comm) != num_ions)
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
		if (pe == pct.gridpe) continue;

		/* Determine if ion has overlap with a given PE becasue of beta functions */
                if(this->type == DELOCALIZED) 
                    map = true;
                else
                    map = test_overlap (pe, iptr, Aix, Aiy, Aiz, this->nldims[iptr->species], PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID);

                /* Determine if ion has overlap with a given PE becasue of Q function */
                if(ct.norm_conserving_pp) {
                    map2 = FALSE;
                }
                else {
                    map2 = test_overlap (pe, iptr, Aix2, Aiy2, Aiz2, sp->qdim, 
                            get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), FNX_GRID, FNY_GRID, FNZ_GRID);
                }

                if((this->kind == ORBITAL_PROJECTOR) && (sp->num_ldaU_orbitals == 0))
                {
                    map = false;
                    map2 = false;
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
    } // for (nlion = 0; nlion < this->num_nonloc_ions; nlion++)


    /* Release temporary memory */
    delete [] Aiz2;
    delete [] Aiy2;
    delete [] Aix2;
    delete [] Aiz;
    delete [] Aiy;
    delete [] Aix;

}

template <class KpointType> int Projector<KpointType>::get_num_tot_proj(void)
{
    return this->num_tot_proj;
}

template <class KpointType> int Projector<KpointType>::get_pstride(void)
{
    return this->pstride;
}

template <class KpointType> int Projector<KpointType>::get_num_nonloc_ions(void)
{
    return this->num_nonloc_ions;
}

template <class KpointType> int Projector<KpointType>::get_num_owned_ions(void)
{
    return this->num_owned_ions;
}

template <class KpointType> int * Projector<KpointType>::get_owned_ions_list(void)
{
    return this->owned_ions_list;
}

template <class KpointType> int * Projector<KpointType>::get_nonloc_ions_list(void)
{
    return this->nonloc_ions_list;
}


template <class KpointType> int Projector<KpointType>::get_nldim(int species)
{
    return this->nldims[species];
}

// Applies projectors to orbitals associated with kpoint kptr
template <class KpointType> void Projector<KpointType>::project(Kpoint<KpointType> *kptr, KpointType *p, int offset, int nstates, KpointType *weight)
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

        int length = factor * ct.num_ions * nstates * this->pstride;
        RmgGemm (transa, transn, this->num_tot_proj, nstates, kptr->pbasis, alpha,
            weight, kptr->pbasis, &kptr->orbital_storage[offset*kptr->pbasis], kptr->pbasis,
            rzero, p, this->num_tot_proj);

        if(pct.grid_npes != 1)
            MPI_Allreduce(MPI_IN_PLACE, (double *)p, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

        return;
    }

    // And now localized
    KpointType *own_buff = NULL, *nown_buff = NULL;
    KpointType *send_buff, *recv_buff;
    MPI_Request *req_send, *req_recv;
    MPI_Request *req_own = NULL, *req_nown = NULL;
    int size_own, size_nown, pe;

    /*Allocate memory for communication */
    size_own = 0;
    for (pe = 0; pe < this->num_owned_pe; pe++)
        size_own += this->num_owned_ions_per_pe[pe];
    size_own *= nstates * this->pstride;

    size_nown = 0;
    for (pe = 0; pe < this->num_owners; pe++)
        size_nown += this->num_nonowned_ions_per_pe[pe];
    size_nown *= nstates * this->pstride;


    if (size_own) {
        own_buff = new KpointType[size_own]();
    }
    if (size_nown) {
        nown_buff = new KpointType[size_nown]();
    }

    if (this->num_owned_pe) {
        req_own = new MPI_Request[this->num_owned_pe];
    }
    if (this->num_owners) {
        req_nown = new MPI_Request[this->num_owners];
    }

    /*First owning cores will receive data from non-owners */
    send_buff = nown_buff;
    recv_buff = own_buff;

    req_send = req_nown;
    req_recv = req_own;

    /*First post non-blocking receives for data about owned ions from cores who do not own the ions */
    betaxpsi_receive (recv_buff, this->num_owned_pe, this->owned_pe_list,
                       this->num_owned_ions_per_pe, req_recv, nstates);


    // Set sint array
    //sint = &kptr->newsint_local[offset*this->num_nonloc_ions*this->pstride];
    KpointType *sint = new KpointType[this->num_nonloc_ions * nstates * this->pstride]();


    /*Loop over ions and calculate local projection between beta functions and wave functions */
    betaxpsi_calculate (kptr, sint, &kptr->orbital_storage[offset*kptr->pbasis], nstates, weight);


    /*Pack data for sending */
    betaxpsi_pack (sint, send_buff, this->num_owners,
                    this->num_nonowned_ions_per_pe, this->list_ions_per_owner, nstates);

    /*Send <beta|psi> contributions  to the owning PE */
    betaxpsi_send (send_buff, this->num_owners, this->owners_list,
                    this->num_nonowned_ions_per_pe, req_send, nstates);

    /*Wait until all data is received */
    if(this->num_owned_pe)
        MPI_Waitall (this->num_owned_pe, req_recv, MPI_STATUSES_IGNORE);
    if(this->num_owners)
        MPI_Waitall (this->num_owners, req_send, MPI_STATUSES_IGNORE);


    /*Unpack received data and sum contributions from all pes for owned ions */
    betaxpsi_sum_owned (recv_buff, sint, this->num_owned_pe,
                         this->num_owned_ions_per_pe, this->list_owned_ions_per_pe, nstates);

    /*In the second stage, owning cores will send summed data to non-owners */
    send_buff = own_buff;
    recv_buff = nown_buff;

    req_send = req_own;
    req_recv = req_nown;
    
    
    /*Receive summed data for non-owned ions from owners */
    betaxpsi_receive (recv_buff, this->num_owners, this->owners_list,
                       this->num_nonowned_ions_per_pe, req_recv, nstates);

    /*Pack summed data for owned ions to send to non-owners */
    betaxpsi_pack (sint, send_buff, this->num_owned_pe,
                    this->num_owned_ions_per_pe, this->list_owned_ions_per_pe, nstates);

    /*Send packed data for owned ions to non-owners */
    betaxpsi_send (send_buff, this->num_owned_pe, this->owned_pe_list,
                    this->num_owned_ions_per_pe, req_send, nstates);

    /*Wait until all data is received */
    if(this->num_owned_pe)
        MPI_Waitall (this->num_owned_pe, req_send, MPI_STATUSES_IGNORE);
    if(this->num_owners)
        MPI_Waitall (this->num_owners, req_recv, MPI_STATUSES_IGNORE);

    /*Finaly, write received data about non-owned ions into sint array */
    betaxpsi_write_non_owned (sint, recv_buff, this->num_owners,
                               this->num_nonowned_ions_per_pe, this->list_ions_per_owner, nstates);


    if (this->num_owned_pe)
        delete [] req_own;
    if (this->num_owners)
        delete [] req_nown;
    if (size_nown)
        delete [] nown_buff;
    if (size_own)
        delete [] own_buff;


    int idx= 0 ;
    for(int st = 0;st < nstates;st++) {
        for(int ion = 0;ion < this->num_nonloc_ions;ion++) {
            for(int ip = 0;ip < this->pstride;ip++) {
                p[idx] = sint[ion*nstates*this->pstride + st*this->pstride + ip];
                idx++;
            }
        }
    }
    delete [] sint;
}


template <class KpointType>
void Projector<KpointType>::betaxpsi_calculate (Kpoint<KpointType> *kptr, KpointType * sint_ptr, KpointType *psi, int num_states, KpointType *weight)
{
    KpointType rzero(0.0);
    char *transt = "t", *transn = "n", *transc = "c";
    char *transa;
   
    transa = transc;
    if(ct.is_gamma) transa = transt;

    KpointType alpha(get_vel());

    if(this->num_tot_proj == 0) return;
    int pbasis = kptr->pbasis;

#if GPU_ENABLED
    KpointType *nlarray = (KpointType *)GpuMallocManaged(sizeof(KpointType) * this->num_tot_proj * num_states);
#else
    KpointType *nlarray = new KpointType[this->num_tot_proj * num_states]();
#endif
    RmgGemm (transa, transn, this->num_tot_proj, num_states, pbasis, alpha, 
            weight, pbasis, psi, pbasis, rzero, nlarray, this->num_tot_proj);

    for (int nion = 0; nion < this->num_nonloc_ions; nion++)
    {
        for(int istate = 0; istate < num_states; istate++)
        {
            for(int ip = 0; ip < this->pstride; ip++)
            {
                int proj_index = nion * this->pstride + ip;
                int idx1 = nion * num_states * this->pstride + istate * this->pstride + ip;
                //idx2 = proj_index * num_states + istate;
                int idx2 = istate * this->num_tot_proj + proj_index;
                sint_ptr[idx1] = nlarray[idx2];
            }
        }
    }

#if GPU_ENABLED
    GpuFreeManaged(nlarray);
#else
    delete [] nlarray;
#endif

}



/*This receives data from other PEs for ions owned by current PE*/
template <class KpointType>
void Projector<KpointType>::betaxpsi_receive (KpointType * recv_buff, int num_pes,
        int *pe_list, int *num_ions_per_pe,
        MPI_Request * req_recv, int num_states)
{
    KpointType *tpr;
    int tag, pe, source, size;

    tpr = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        source = pe_list[pe];
        /*Tag is sender_rank * pct.grid_npes + receiver_rank */
        tag = 111;
        size = num_ions_per_pe[pe] * num_states * this->pstride;
        int transfer_size = size;
        if(!ct.is_gamma) transfer_size *= 2;
        MPI_Irecv (tpr, transfer_size, MPI_DOUBLE, source, tag, pct.grid_comm, &req_recv[pe]);
        tpr += size;
    }
}

template <class KpointType>
void Projector<KpointType>::betaxpsi_send (KpointType * send_buff, int num_pes,
        int *pe_list, int *num_ions_per_pe,
        MPI_Request * req_send, int num_states)
{
    KpointType *tpr;
    int target, num_ions, size, tag, pe;

    tpr = send_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        target = pe_list[pe];
        /*Tag is sender_rank * pct.grid_npes + receiver_rank */
        tag = 111;
        num_ions = num_ions_per_pe[pe];
        size = num_ions * num_states * this->pstride;
        int transfer_size = size;
        if(!ct.is_gamma) transfer_size *= 2;
        MPI_Isend (tpr, transfer_size, MPI_DOUBLE, target, tag, pct.grid_comm, &req_send[pe]);
        tpr += size;
    }
}


template <class KpointType>
void Projector<KpointType>::betaxpsi_pack (KpointType * sint, KpointType * fill_buff, 
        int num_pes, int *num_ions_per_pe,
        int_2d_array &list_ions_per_pe, int num_states)
{
    KpointType *tpr_buff, *sint_tpr;
    int size, num_ions, ion, nlion, pe;

    tpr_buff = fill_buff;

    size = num_states * this->pstride;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        /*Loop over ions that need to be sent to pe */
        for (ion = 0; ion < num_ions; ion++)
        {
            nlion = list_ions_per_pe[pe][ion];
            sint_tpr = &sint[nlion * num_states * this->pstride];

            /*Pack data into send array */
            for(int idx=0;idx < size;idx++) tpr_buff[idx] = sint_tpr[idx];
            tpr_buff += size;
        }
    }

}


template <class KpointType>
void Projector<KpointType>::betaxpsi_sum_owned (KpointType * recv_buff, KpointType * sint,
        int num_pes, int *num_ions_per_pe,
        int_2d_array &list_ions_per_pe, int num_states)
{
    KpointType *tpr1, *tpr2;
    int size, num_ions, ion_index, pe, ion;

    size = num_states * this->pstride;
    tpr1 = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];
        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];
            tpr2 = &sint[ion_index * num_states * this->pstride];
            for(int idx=0;idx < size;idx++) tpr2[idx] += tpr1[idx];
            tpr1 += size;
        }
    }


}



template <class KpointType>
void Projector<KpointType>::betaxpsi_write_non_owned (KpointType * sint, KpointType * recv_buff, int num_pes,
        int *num_ions_per_pe,
        int_2d_array &list_ions_per_pe, int num_states)
{
    KpointType *tpr1, *tpr2;
    int size, num_ions, ion_index, pe, ion;

    size = num_states * this->pstride;
    tpr1 = recv_buff;

    for (pe = 0; pe < num_pes; pe++)
    {
        num_ions = num_ions_per_pe[pe];

        for (ion = 0; ion < num_ions; ion++)
        {
            ion_index = list_ions_per_pe[pe][ion];

            tpr2 = &sint[ion_index * num_states * this->pstride];
            for(int idx=0;idx < size;idx++) tpr2[idx] = tpr1[idx];

            tpr1 += size;
        }
    }
}



// Destructor
template <class KpointType> Projector<KpointType>::~Projector(void)
{
    this->nlcrds.empty();
    delete [] this->nldims;
    delete [] this->idxptrlen;
    delete [] this->nonloc_ions_list;
    delete [] this->num_nonowned_ions_per_pe;
    delete [] this->owners_list;
    delete [] this->num_owned_ions_per_pe;
    delete [] this->owned_pe_list;
    delete [] this->owned_ions_list;
}
