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


#include "make_conf.h"
#include "const.h"
#include "grid.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include <complex> 
#include "Kpoint.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"
#include "GpuAlloc.h"
#include "rmg_error.h"
#include "ErrorFuncs.h"

/*Call to this function needs to be preceeded by get_QI, since we use pct.Qidxptrlen,
 * which is setup in that function*/

static void reset_pct_arrays (int ion);

template void GetNlop<double> (Kpoint<double> **);
template void GetNlop<std::complex<double> > (Kpoint<std::complex<double>> **);

template <typename KpointType>
void GetNlop (Kpoint<KpointType> **Kptr)
{

    int *Aix, *Aiy, *Aiz;
    int *Aix2, *Aiy2, *Aiz2;
    int P0_BASIS = get_P0_BASIS();

    /*Reset number of nonlocal ions */
    pct.num_nonloc_ions = 0;
    pct.num_owned_ions = 0;
    pct.num_owned_pe = 0;
    pct.num_owners = 0;

    // Allocate memory for some lists
    if(!pct.owned_pe_list) pct.owned_pe_list = new int[NPES]();
    if(!pct.num_owned_ions_per_pe) pct.num_owned_ions_per_pe = new int[NPES]();
    if(!pct.owners_list) pct.owners_list = new int[NPES]();
    if(!pct.num_nonowned_ions_per_pe) pct.num_nonowned_ions_per_pe = new int[ct.num_ions]();
    if(!pct.nonloc_ions_list) pct.nonloc_ions_list = new int[ct.num_ions]();
    if(!pct.nonloc_ion_ownflag) pct.nonloc_ion_ownflag = new int[ct.num_ions]();


    /* Grab some memory for temporary storage */
    int alloc = ct.max_nlpoints;

    if (alloc <get_NX_GRID() * get_NY_GRID() * get_NZ_GRID())
        alloc =get_NX_GRID() * get_NY_GRID() * get_NZ_GRID();

    reset_pct_arrays (ct.num_ions);
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
    int nthreads = ct.THREADS_PER_NODE;
    nthreads = std::min(ct.THREADS_PER_NODE, 4);

#pragma omp parallel private(Aix, Aiy, Aiz, Aix2, Aiy2, Aiz2)
{

    Aix = new int[NX_GRID];
    Aiy = new int[NY_GRID];
    Aiz = new int[NZ_GRID];

#pragma omp parallel for schedule(static, 1) num_threads(nthreads)
    /* Loop over ions */
    for (int ion = 0; ion < ct.num_ions; ion++)
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
        if(!ct.localize_projectors) {
            map = true;
            nlxdim = NX_GRID;
            nlydim = NY_GRID;
            nlzdim = NZ_GRID;
            iptr->nlxcstart = iptr->nlycstart = iptr->nlzcstart = 0.0;
            ilow = PX_OFFSET;
            ihi = ilow + get_PX0_GRID() -1;
            jlow = PY_OFFSET;
            jhi = jlow + get_PY0_GRID() -1;
            klow = PZ_OFFSET;
            khi = klow + get_PZ0_GRID() -1;
        }
        else
            map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         sp->nldim, PX0_GRID, PY0_GRID, PZ0_GRID,
                         NX_GRID, NY_GRID, NZ_GRID,
                         &iptr->nlxcstart, &iptr->nlycstart, &iptr->nlzcstart);

        /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
        /*xtal vector between ion and left bottom corner of the box */
        vect[0] = iptr->xtal[0] - iptr->nlxcstart;
        vect[1] = iptr->xtal[1] - iptr->nlycstart;
        vect[2] = iptr->xtal[2] - iptr->nlzcstart;

        /*Substract vector between left bottom corner of the box and center of the box */
        vect[0] -= (nlxdim / 2) / (double) NX_GRID;
        vect[1] -= (nlydim / 2) / (double) NY_GRID;
        vect[2] -= (nlzdim / 2) / (double) NZ_GRID;

        /*The vector we are looking for should be */
        to_cartesian (vect, iptr->nlcrds);

            pct.nl_flag[ion] = map;

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
                        if(ct.localize_projectors) {

                            if ((((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                                 ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                                 ((Aiz[iz] >= klow) && (Aiz[iz] <= khi))))
                            {
                                /* Cut it off if required */
                                int itmp =
                                    (ix - icenter) * (ix - icenter) +
                                    (iy - icenter) * (iy - icenter) + (iz - icenter) * (iz - icenter);

                                if (icut >= itmp)
                                {

                                    icount++;

                                }

                            }
                        }
                        else {
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
            pct.idxptrlen[ion] = icount;



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
        if (pct.idxptrlen[ion] || pct.Qidxptrlen[ion])
        {

            if (pct.num_nonloc_ions >= MAX_NONLOC_IONS)
                rmg_error_handler(__FILE__, __LINE__, "Too many nonlocal ions. pct.nonloc_ions_list will overflow");

            pct.nonloc_ions_list[pct.num_nonloc_ions] = ion;

            /*Ownership flag for current ion */
            pct.nonloc_ion_ownflag[pct.num_nonloc_ions] = 0;

            /*See if this processor owns current ion */
            if (pct.gridpe == claim_ion (iptr->xtal, PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID))
            {
                if (pct.num_owned_ions >= MAX_NONLOC_IONS)
                    rmg_error_handler (__FILE__, __LINE__, "Too many owned ions, pct.owned_ions_list will overflow");

                pct.owned_ions_list[pct.num_owned_ions] = ion;
                pct.num_owned_ions++;

                pct.nonloc_ion_ownflag[pct.num_nonloc_ions] = 1;
            }

            pct.num_nonloc_ions++;
        }

    }                           /* end for (ion = 0; ion < ct.num_ions; ion++) */

    Aix = new int[NX_GRID]();
    Aiy = new int[NY_GRID]();
    Aiz = new int[NZ_GRID]();

    Aix2 = new int[FNX_GRID]();
    Aiy2 = new int[FNY_GRID]();
    Aiz2 = new int[FNZ_GRID]();

    pct.num_tot_proj = pct.num_nonloc_ions * ct.max_nl;

    int known_nonowner, nonwner_index, known_owner, owner_index, owned_ions_per_pe, nonowned_ions_per_pe; 
    size_t weight_size = pct.num_tot_proj * P0_BASIS + 128;
    int owned, owner;

#if GPU_ENABLED
    if( cudaSuccess != cudaMallocHost((void **)&pct.weight, weight_size * sizeof(double) ))
        rmg_error_handler(__FILE__,__LINE__,"Error: cudaMallocHost failed for: get_nlop \n");
    for(int idx = 0;idx < weight_size;idx++) pct.weight[idx] = 0.0;

    if(ct.need_Bweight) 
    {
        if( cudaSuccess != cudaMallocHost((void **)&pct.Bweight, weight_size * sizeof(double) ))
            rmg_error_handler(__FILE__,__LINE__,"Error: cudaMallocHost failed for: get_nlop \n");
        for(int idx = 0;idx < weight_size;idx++) pct.Bweight[idx] = 0.0;
    }
    else {
        pct.Bweight = pct.weight;
    }
#else
    pct.weight = new double[weight_size]();
    if(ct.need_Bweight) {
        pct.Bweight = new double[weight_size]();
    }
    else {
        pct.Bweight = pct.weight;
    }
#endif



    /*Make sure that ownership of ions is properly established
     * This conditional can be removed if it is found that claim_ions works reliably*/
    if (int_sum_all (pct.num_owned_ions, pct.grid_comm) != ct.num_ions)
        rmg_error_handler (__FILE__, __LINE__, "Problem with claimimg ions.");

  
    if (ct.verbose)
    {
	rmg_printf ("\n PE %d: Number of nonlocal ions is %d", pct.gridpe, pct.num_nonloc_ions);

	for (int i = 0; i < pct.num_nonloc_ions; i++)
	    rmg_printf (" %d", pct.nonloc_ions_list[i]);

	rmg_printf ("\n PE %d: Number of claimed ions is %d", pct.gridpe, pct.num_owned_ions);
	for (int i = 0; i < pct.num_owned_ions; i++)
	    rmg_printf (" %d", pct.owned_ions_list[i]);
    }

    
    
    // Set storage sequentially for real and imaginary components so we can transform storage pattern
#if GPU_ENABLED
    if (pct.newsintR_local)
        cudaFreeHost(pct.newsintR_local);
#else
    if (pct.newsintR_local)
        delete [] pct.newsintR_local;
#endif
   
    int factor = 2;
    if(ct.is_gamma) factor = 1; 
    size_t sint_alloc = (size_t)(factor * ct.num_kpts_pe * pct.num_nonloc_ions * ct.max_nl);
    sint_alloc *= (size_t)ct.max_states;
#if GPU_ENABLED
    cudaError_t custat;
    custat = cudaMallocHost((void **)&pct.newsintR_local , sint_alloc * sizeof(double));
    RmgCudaError(__FILE__, __LINE__, custat, "Error: cudaHostMalloc failed in InitGpuHostMalloc\n");
#else
    pct.newsintR_local = new double[sint_alloc]();
#endif

    KpointType *tsintnew_ptr = (KpointType *)pct.newsintR_local;

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++){
        Kptr[kpt]->sint_size = (size_t)pct.num_nonloc_ions * (size_t)ct.max_states * (size_t)ct.max_nl;
        Kptr[kpt]->newsint_local = tsintnew_ptr;

        tsintnew_ptr += Kptr[kpt]->sint_size;
    }
    

    /* Loop over all ions to obtain the lists necessary for communication */
    for (int nlion = 0; nlion < pct.num_nonloc_ions; nlion++)
    {

        int ion = pct.nonloc_ions_list[nlion];
        int map, map2;

        /* Generate ion pointer */
        ION *iptr = &ct.ions[ion];

        /* Get species type */
        SPECIES *sp = &ct.sp[iptr->species];


        /*See if this processor owns current ion */
        owner = claim_ion (iptr->xtal, PX0_GRID, PY0_GRID, PZ0_GRID, NX_GRID, NY_GRID, NZ_GRID);
	//Dprintf ("Owner of ion %d nlion %d is PE %d", ion, nlion, owner); 
        
	owned = 0;
        if (pct.gridpe == owner)
            owned = 1;
	//Dprintf("owned set to %d for ion %d", owned, ion);



	/*Case when current PE owns an ion, in such a case we need to assemble
	 * list of other PEs that have overlap - those non-owning PEs will be sending data*/
	if (owned)
	{
	    //Dprintf("ion %d is owned", ion);
	    
	    /*Loop over all processors in caculation to determine if it overlaps with current ion */
	    /* This should only loop over processors in a grid */
	    for (int pe = 0; pe < NPES; pe++)
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

                //Dprintf("Overlap condition for ion %d and PE %d is %d, map is %d, map2 is %d ", ion, pe, map || map2, map, map2); 

                if (map || map2)
                {


                    /* See if this non-owning pe has already been recorded*/
                    known_nonowner = 0;
                    for (int i = 0; i < pct.num_owned_pe; i++)
                    {
                        if (pct.owned_pe_list[i] == pe)
                        {
                            known_nonowner++;
                            nonwner_index = i;
                            break;
                        }
                    }

                    //Dprintf("known_nonowner is %d and nonwner_index is %d for PE %d and ion %d", known_nonowner, nonwner_index, pe, ion);

                    /* Record pe index into communication list, if it has not been recorded yet*/
                    if (!known_nonowner)
                    {
                        if (pct.num_owned_pe >= MAX_NONLOC_PROCS)
                            rmg_error_handler (__FILE__, __LINE__, "pct.num_nonloc_pes is too large.");

                        pct.owned_pe_list[pct.num_owned_pe] = pe;

                        pct.num_owned_ions_per_pe[pct.num_owned_pe] = 1;

                        pct.list_owned_ions_per_pe[pct.num_owned_pe][0] = nlion;

                        pct.num_owned_pe++;

                        //Dprintf("Recording previously unknown PE %d, pct.owned_pe_list[%d]=%d,  pct.num_owned_ions_per_pe[%d]=%d, pct.list_owned_ions_per_pe[%d][0]=%d,  pct.num_owned_pe = %d for ion %d", 
                        //		pe, pct.num_owned_pe, pct.owned_pe_list[pct.num_owned_pe],pct.num_owned_pe, pct.num_owned_ions_per_pe[pct.num_owned_pe], pct.num_owned_pe, pct.list_owned_ions_per_pe[pct.num_owned_pe][0], pct.num_owned_pe, ion);

                    }

                    else
                    {
                        /*Number of ions to communicate with current non-owner PE*/
                        owned_ions_per_pe = pct.num_owned_ions_per_pe[nonwner_index];

                        if (owned_ions_per_pe >= MAX_NONLOC_IONS)
                            rmg_error_handler (__FILE__, __LINE__, "pct.num_owned_ions_per_pe too large.");

                        pct.list_owned_ions_per_pe[nonwner_index][owned_ions_per_pe] = nlion;
                        pct.num_owned_ions_per_pe[nonwner_index]++;

                        //Dprintf("ion %d: PE %d is already known, pct.list_owned_ions_per_pe[%d][%d]=%d, pct.num_owned_ions_per_pe[%d]=%d",
                        //	ion, pe, nonwner_index, owned_ions_per_pe, pct.list_owned_ions_per_pe[nonwner_index][owned_ions_per_pe],
                        //	nonwner_index, pct.num_owned_ions_per_pe[nonwner_index]);
                    }
                }

            }  /* End of for (pe=0; pe < NPES; pe++) */

        }            

        /*For non-owned ion we need to record rank of owner */
        else	
        {


            /*Loop over list of processors that we already have to see if owning "pe" has alredy been recorded */
            known_owner = 0;
            for (int i = 0; i < pct.num_owners; i++)
            {
                if (pct.owners_list[i] == owner)
                {
                    known_owner++;
                    owner_index = i;
                    break;
                }
            }

            //Dprintf("known_owner is %d and owner_index is %d for ion %d", known_owner, owner_index, ion);


            if (!known_owner)
            {
                if (pct.num_owners >= MAX_NONLOC_PROCS) {
                    rmg_error_handler (__FILE__, __LINE__, "pct.num_owners (%d) is too large (max: %d)");
                }

                pct.owners_list[pct.num_owners] = owner;

                pct.num_nonowned_ions_per_pe[pct.num_owners] = 1;

                pct.list_ions_per_owner[pct.num_owners][0] = nlion;

                pct.num_owners++;

                //rmg_printf("Owner of ion %d PE %d has not been previously recorded: pct.owners_list[%d]=%d, pct.num_nonowned_ions_per_pe[%d]=%d, pct.list_ions_per_owner[%d][0]=%d, pct.num_owners=%d", 
                //	ion, owner, pct.num_owners, pct.owners_list[pct.num_owners], pct.num_owners, pct.num_nonowned_ions_per_pe[pct.num_owners], pct.num_owners, pct.list_ions_per_owner[pct.num_owners][0], pct.num_owners); 
            }

            else
            {

                nonowned_ions_per_pe = pct.num_nonowned_ions_per_pe[owner_index];

                if (nonowned_ions_per_pe >= MAX_NONLOC_IONS)
                    rmg_error_handler(__FILE__, __LINE__, "pct.num_nonowned_ions_per_pe too large.");

                pct.list_ions_per_owner[owner_index][nonowned_ions_per_pe] = nlion;
                pct.num_nonowned_ions_per_pe[owner_index]++;

                //Dprintf("ion %d owned by %d: pct.list_ions_per_owner[%d][%d]=%d, pct.num_nonowned_ions_per_pe[%d]=%d", 
                //		ion, owner, owner_index, nonowned_ions_per_pe, pct.list_ions_per_owner[owner_index][nonowned_ions_per_pe], owner_index, pct.num_nonowned_ions_per_pe[owner_index]);
            }

        }

    }/*for (nlion = 0; nlion < pct.num_nonloc_ions; nlion++)*/


    /* Release temporary memory */
    delete [] Aiz2;
    delete [] Aiy2;
    delete [] Aix2;
    delete [] Aiz;
    delete [] Aiy;
    delete [] Aix;

} 



/*Resets pct projector arrays, deallocates memory*/
static void reset_pct_arrays (int num_ions)
{

    int ion;
    for(ion = 0; ion < num_ions; ion++)
    {
        pct.idxptrlen[ion] = 0;

    }



    if (pct.weight != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.weight);
#else
        delete [] pct.weight;
#endif
    }
    if ((pct.Bweight != NULL) && ct.need_Bweight) {
#if GPU_ENABLED
        cudaFreeHost(pct.Bweight);
#else
        delete [] pct.Bweight;
#endif
    }

}
