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


/*Call to this function needs to be preceeded by get_QI, since we use pct.Qidxptrlen,
 * which is setup in that function*/

static void reset_pct_arrays (int ion);

template void GetNlop<double> (Kpoint<double> **);
template void GetNlop<std::complex<double> > (Kpoint<std::complex<double>> **);

template <typename KpointType>
void GetNlop (Kpoint<KpointType> **Kptr)
{

    int ion, idx, i, pe, owned, nlion, owner;
    int ix, iy, iz, prj_per_ion;
    int *pvec, *dvec, *ivec;
    int ilow, jlow, klow, ihi, jhi, khi, map, map2, icount;
    int alloc;
    int *Aix, *Aiy, *Aiz;
    int *Aix2, *Aiy2, *Aiz2;
    int icut, itmp, icenter;
    double vect[3];
    SPECIES *sp;
    ION *iptr;
    int known_nonowner, nonwner_index, known_owner, owner_index, owned_ions_per_pe, nonowned_ions_per_pe; 


    /*Reset number of nonlocal ions */
    pct.num_nonloc_ions = 0;
    pct.num_owned_ions = 0;
    pct.num_owned_pe = 0;
    pct.num_owners = 0;


    /* Grab some memory for temporary storage */
    alloc = ct.max_nlpoints;
    if (alloc <get_P0_BASIS())
        alloc =get_P0_BASIS();
    pvec = new int[2*alloc];
    dvec = pvec + alloc;

    Aix = new int[get_NX_GRID()];
    Aiy = new int[get_NY_GRID()];
    Aiz = new int[get_NZ_GRID()];

    Aix2 = new int[get_FNX_GRID()];
    Aiy2 = new int[get_FNY_GRID()];
    Aiz2 = new int[get_FNZ_GRID()];


    reset_pct_arrays (ct.num_ions);

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Release memory and reset other parameters for given ion */

        /* Generate ion pointer */
        iptr = &ct.ions[ion];


        /* Get species type */
        sp = &ct.sp[iptr->species];

        prj_per_ion = sp->nh;


        icenter = sp->nldim / 2;
        icut = (icenter + 1) * (icenter + 1);


        /* Determine mapping indices or even if a mapping exists */
        map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         sp->nldim, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(),
                         get_NX_GRID(), get_NY_GRID(), get_NZ_GRID(),
                         &iptr->nlxcstart, &iptr->nlycstart, &iptr->nlzcstart);

        /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
        /*xtal vector between ion and left bottom corner of the box */
        vect[0] = iptr->xtal[0] - iptr->nlxcstart;
        vect[1] = iptr->xtal[1] - iptr->nlycstart;
        vect[2] = iptr->xtal[2] - iptr->nlzcstart;

        /*Substract vector between left bottom corner of the box and center of the box */
        vect[0] -= (sp->nldim / 2) / (double) get_NX_GRID();
        vect[1] -= (sp->nldim / 2) / (double) get_NY_GRID();
        vect[2] -= (sp->nldim / 2) / (double) get_NZ_GRID();

        /*The vector we are looking for should be */
        to_cartesian (vect, iptr->nlcrds);


        /* If there is a mapping for this ion then we have to generate */
        /* the projector.                                              */

        for (idx = 0; idx < alloc; idx++)
            dvec[idx] = pvec[idx] = 0;

        icount = 0;
        if (map)
        {

            /* Generate index arrays */
            icount = idx = 0;
            for (ix = 0; ix < sp->nldim; ix++)
            {
                for (iy = 0; iy < sp->nldim; iy++)
                {
                    for (iz = 0; iz < sp->nldim; iz++)
                    {
                        dvec[idx] = FALSE;
                        if ((((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                             ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                             ((Aiz[iz] >= klow) && (Aiz[iz] <= khi))))
                        {
                            /* Cut it off if required */
                            itmp =
                                (ix - icenter) * (ix - icenter) +
                                (iy - icenter) * (iy - icenter) + (iz - icenter) * (iz - icenter);

                            if (icut >= itmp)
                            {

                                pvec[icount] =
                                    get_PY0_GRID() * get_PZ0_GRID() * ((Aix[ix]-get_PX_OFFSET()) % get_PX0_GRID()) +
                                    get_PZ0_GRID() * ((Aiy[iy]-get_PY_OFFSET()) % get_PY0_GRID()) + 
                                    ((Aiz[iz]-get_PZ_OFFSET()) % get_PZ0_GRID());

                                dvec[idx] = TRUE;

                                icount++;

                            }

                        }
                        idx++;
                    }
                }
            }


            /* Save number of points */
            pct.idxptrlen[ion] = icount;


            /* Now we have to allocate memory for the index array */
            if (icount)
            {
                pct.idxflag[ion] = new int[idx];

                ivec = pct.idxflag[ion];
                for (ix = 0; ix < idx; ix++)
                    ivec[ix] = (int) dvec[ix];
            }

            else
                pct.idxflag[ion] = NULL;



            if (icount)
            {
                pct.nlindex[ion] = new int[icount + 128]();

                ivec = pct.nlindex[ion];
                for (idx = 0; idx < icount; idx++)
                    ivec[idx] = (int) pvec[idx];
            }
            else
                pct.nlindex[ion] = NULL;



            /* Allocate memory for the phase array */
            if ((icount * prj_per_ion)) {
                itmp = sp->nldim * sp->nldim * sp->nldim;
                pct.phaseptr[ion] = new double[2 * itmp * ct.num_kpts + 128]();
            }
            else {
                pct.phaseptr[ion] = NULL;
            }
            get_phase (iptr, pct.phaseptr[ion], get_P0_BASIS(), dvec);

        }                       /* end if (map) */


        /*Add ion into list of nonlocal ions if it has overlap with given processor */
        if (icount || pct.Qidxptrlen[ion])
        {

            if (pct.num_nonloc_ions >= MAX_NONLOC_IONS)
                rmg_error_handler(__FILE__, __LINE__, "Too many nonlocal ions. pct.nonloc_ions_list will overflow");

            pct.nonloc_ions_list[pct.num_nonloc_ions] = ion;

            /*Ownership flag for current ion */
            pct.nonloc_ion_ownflag[pct.num_nonloc_ions] = 0;

            /*See if this processor owns current ion */
            if (pct.gridpe ==
                claim_ion (iptr->xtal, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_NX_GRID(), get_NY_GRID(),
                           get_NZ_GRID()))
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

    pct.num_tot_proj = pct.num_nonloc_ions * ct.max_nl;

//    pct.num_tot_proj = 0;
//    for(ion = 0; ion <pct.num_nonloc_ions; ion++)
//    {
//        ion1 = pct.nonloc_ions_list[ion];
//
//        iptr = &ct.ions[ion1];
//
        /* Get species type */
//        sp = &ct.sp[iptr->species];

//        prj_per_ion = sp->nh;
//        pct.num_tot_proj += prj_per_ion;
//    }
        
    size_t weight_size = pct.num_tot_proj * get_P0_BASIS() + 128;

#if GPU_ENABLED
    if( cudaSuccess != cudaMallocHost((void **)&pct.weight, weight_size * sizeof(double) ))
        rmg_error_handler(__FILE__,__LINE__,"Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.weight[idx] = 0.0;
    if( cudaSuccess != cudaMallocHost((void **)&pct.Bweight, weight_size * sizeof(double) ))
        rmg_error_handler(__FILE__,__LINE__,"Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.Bweight[idx] = 0.0;

    weight_size = pct.num_tot_proj * pct.num_tot_proj + 128;
    if( cudaSuccess != cudaMallocHost((void **)&pct.M_dnm, weight_size * sizeof(double) ))
        rmg_error_handler(__FILE__,__LINE__,"Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.M_dnm[idx] = 0.0;

    if( cudaSuccess != cudaMallocHost((void **)&pct.M_qqq, weight_size * sizeof(double) ))
        rmg_error_handler(__FILE__,__LINE__,"Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.M_dnm[idx] = 0.0;

#else
    pct.weight = new double[weight_size]();
    pct.Bweight = new double[weight_size]();

    weight_size = pct.num_tot_proj * pct.num_tot_proj + 128;
    pct.M_dnm = new double[weight_size]();
    pct.M_qqq = new double[weight_size]();
#endif



    /*Make sure that ownership of ions is properly established
     * This conditional can be removed if it is found that claim_ions works reliably*/
    if (int_sum_all (pct.num_owned_ions, pct.grid_comm) != ct.num_ions)
        rmg_error_handler (__FILE__, __LINE__, "Problem with claimimg ions.");

    rmg_printf ("\n PE %d: Number of nonlocal ions is %d", pct.gridpe, pct.num_nonloc_ions);

    for (i = 0; i < pct.num_nonloc_ions; i++)
        rmg_printf (" %d", pct.nonloc_ions_list[i]);

    rmg_printf ("\n PE %d: Number of claimed ions is %d", pct.gridpe, pct.num_owned_ions);
    for (i = 0; i < pct.num_owned_ions; i++)
        rmg_printf (" %d", pct.owned_ions_list[i]);

    
    
    // Set storage sequentially for real and imaginary components so we can transform storage pattern
    if (pct.newsintR_local)
        delete [] pct.newsintR_local;
    if (pct.oldsintR_local)
        delete [] pct.oldsintR_local;
    
    int factor = 2;
    if(ct.is_gamma) factor = 1; 
    pct.newsintR_local = new double[factor * ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl]();
    pct.oldsintR_local = new double[factor * ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl]();

    KpointType *tsintnew_ptr = (KpointType *)pct.newsintR_local;
    KpointType *tsintold_ptr = (KpointType *)pct.oldsintR_local;
    for(int kpt = 0;kpt < ct.num_kpts;kpt++) {
        Kptr[kpt]->sint_size = pct.num_nonloc_ions * ct.num_states * ct.max_nl;
        Kptr[kpt]->newsint_local = tsintnew_ptr;
        Kptr[kpt]->oldsint_local = tsintold_ptr;
        tsintnew_ptr += pct.num_nonloc_ions * ct.num_states * ct.max_nl;
        tsintold_ptr += pct.num_nonloc_ions * ct.num_states * ct.max_nl;
    }
    

#if 1

    /* Loop over all ions to obtain the lists necessary for communication */
    for (nlion = 0; nlion < pct.num_nonloc_ions; nlion++)
    {

        ion = pct.nonloc_ions_list[nlion];

        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];


        /*See if this processor owns current ion */
        owner =
            claim_ion (iptr->xtal, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_NX_GRID(), get_NY_GRID(),
                       get_NZ_GRID());
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
	    for (pe = 0; pe < NPES; pe++)
	    {
		
		/*Skip the case when pe is equal to the rank of current processor, in that case no send or
		 * receive is necessary*/
		if (pe == pct.gridpe)
		    continue;

		/* Determine if ion has overlap with a given PE becasue of beta functions */
		map = test_overlap (pe, iptr, Aix, Aiy, Aiz, sp->nldim,
			get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_NX_GRID(), get_NY_GRID(), get_NZ_GRID());

		/* Determine if ion has overlap with a given PE becasue of Q function */
                if(ct.norm_conserving_pp) {
                    map2 = FALSE;
                }
                else {
                    map2 = test_overlap (pe, iptr, Aix2, Aiy2, Aiz2, sp->qdim, 
                        get_FPX0_GRID(), get_FPY0_GRID(), get_FPZ0_GRID(), get_FNX_GRID(), get_FNY_GRID(), get_FNZ_GRID());
                }

		//Dprintf("Overlap condition for ion %d and PE %d is %d, map is %d, map2 is %d ", ion, pe, map || map2, map, map2); 

		if (map || map2)
		{


		    /* See if this non-owning pe has already been recorded*/
		    known_nonowner = 0;
		    for (i = 0; i < pct.num_owned_pe; i++)
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
	    for (i = 0; i < pct.num_owners; i++)
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

#endif

    }/*for (nlion = 0; nlion < pct.num_nonloc_ions; nlion++)*/


    /* Release temporary memory */
    delete [] Aiz2;
    delete [] Aiy2;
    delete [] Aix2;
    delete [] Aiz;
    delete [] Aiy;
    delete [] Aix;

    delete [] pvec;

} 



/*Resets pct projector arrays, deallocates memory*/
static void reset_pct_arrays (int num_ions)
{

    int ion;
    for(ion = 0; ion < num_ions; ion++)
    {
        pct.idxptrlen[ion] = 0;

        if (pct.idxflag[ion]) {
            delete [] pct.idxflag[ion];
            pct.idxflag[ion] = NULL;
        }

        if (pct.nlindex[ion]) {
            delete [] pct.nlindex[ion];
            pct.nlindex[ion] = NULL;
        }

        if (pct.phaseptr[ion]) {
            delete [] pct.phaseptr[ion];
            pct.phaseptr[ion] = NULL;
        }
    }



    if (pct.weight != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.weight);
#else
        delete [] pct.weight;
#endif
    }
    if (pct.Bweight != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.Bweight);
#else
        delete [] pct.Bweight;
#endif
    }

    if (pct.M_dnm != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.M_dnm);
#else
        delete [] pct.M_dnm;
#endif
    }

    if (pct.M_qqq != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.M_qqq);
#else
        delete [] pct.M_qqq;
#endif
    }


}
