/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "grid.h"
#include "common_prototypes.h"
#include "main.h"

/*Call to this function needs to be preceeded by get_QI, since we use pct.Qidxptrlen,
 * which is setup in that function*/

static void reset_pct_arrays (int ion);

void get_nlop (void)
{

    int ion, ion1, idx, i, pe, owned, nlion, owner;
    int ix, iy, iz, ip, prj_per_ion;
    int *pvec, *dvec, *ivec;
    int ilow, jlow, klow, ihi, jhi, khi, map, map2, icount;
    int alloc;
    int *Aix, *Aiy, *Aiz;
    int *Aix2, *Aiy2, *Aiz2;
    int icut, itmp, icenter;
    rmg_double_t vect[3];
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
    my_malloc (pvec, 2*alloc, int);
    dvec = pvec + alloc;

    my_malloc (Aix, get_NX_GRID(), int);
    my_malloc (Aiy, get_NY_GRID(), int);
    my_malloc (Aiz, get_NZ_GRID(), int);

    my_malloc (Aix2, get_FNX_GRID(), int);
    my_malloc (Aiy2, get_FNY_GRID(), int);
    my_malloc (Aiz2, get_FNZ_GRID(), int);

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
        vect[0] -= (sp->nldim / 2) / (rmg_double_t) get_NX_GRID();
        vect[1] -= (sp->nldim / 2) / (rmg_double_t) get_NY_GRID();
        vect[2] -= (sp->nldim / 2) / (rmg_double_t) get_NZ_GRID();

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
                my_malloc (pct.idxflag[ion], idx, int);

                ivec = pct.idxflag[ion];
                for (ix = 0; ix < idx; ix++)
                    ivec[ix] = (int) dvec[ix];
            }

            else
                pct.idxflag[ion] = NULL;



            if (icount)
            {
                my_calloc (pct.nlindex[ion], icount + 128, int);

                ivec = pct.nlindex[ion];
                for (idx = 0; idx < icount; idx++)
                    ivec[idx] = (int) pvec[idx];
            }
            else
                pct.nlindex[ion] = NULL;



            if ((icount * prj_per_ion))
            {

            size_t weight_size = prj_per_ion * get_P0_BASIS() + 128;
#if FDIFF_BETA
                my_calloc (pct.weight_derx[ion], weight_size, rmg_double_t);
                my_calloc (pct.weight_dery[ion], weight_size, rmg_double_t);
                my_calloc (pct.weight_derz[ion], weight_size, rmg_double_t);
#endif
            }
            else
            {
#if FDIFF_BETA
                pct.weight_derx[ion] = NULL;
                pct.weight_dery[ion] = NULL;
                pct.weight_derz[ion] = NULL;
#endif
            }


#if !GAMMA_PT

            /* Allocate memory for the phase array */
            if ((icount * prj_per_ion))
                my_calloc (pct.phaseptr[ion], 2 * icount * ct.num_kpts + 128, rmg_double_t);
            else
                pct.phaseptr[ion] = NULL;

            get_phase (iptr, pct.phaseptr[ion], ip, icount, dvec);

#endif

        }                       /* end if (map) */


        /*Add ion into list of nonlocal ions if it has overlap with given processor */
        if (icount || pct.Qidxptrlen[ion])
        {

            if (pct.num_nonloc_ions >= MAX_NONLOC_IONS)
                error_handler
                    ("Too many nonlocal ions (pct.num_nonloc_ions = %d MAX = %d ), pct.nonloc_ions_list will overflow",
                     pct.num_nonloc_ions, MAX_NONLOC_IONS);

            pct.nonloc_ions_list[pct.num_nonloc_ions] = ion;

            /*Ownership flag for current ion */
            pct.nonloc_ion_ownflag[pct.num_nonloc_ions] = 0;

            /*See if this processor owns current ion */
            if (pct.gridpe ==
                claim_ion (iptr->xtal, get_PX0_GRID(), get_PY0_GRID(), get_PZ0_GRID(), get_NX_GRID(), get_NY_GRID(),
                           get_NZ_GRID()))
            {
                if (pct.num_owned_ions >= MAX_NONLOC_IONS)
                    error_handler ("Too many owned ions, pct.owned_ions_list will overflow");

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
    if( cudaSuccess != cudaMallocHost((void **)&pct.weight, weight_size * sizeof(rmg_double_t) ))
        error_handler("Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.weight[idx] = 0.0;
    if( cudaSuccess != cudaMallocHost((void **)&pct.Bweight, weight_size * sizeof(rmg_double_t) ))
        error_handler("Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.Bweight[idx] = 0.0;

    weight_size = pct.num_tot_proj * pct.num_tot_proj + 128;
    if( cudaSuccess != cudaMallocHost((void **)&pct.M_dnm, weight_size * sizeof(rmg_double_t) ))
        error_handler("Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.M_dnm[idx] = 0.0;

    if( cudaSuccess != cudaMallocHost((void **)&pct.M_qqq, weight_size * sizeof(rmg_double_t) ))
        error_handler("Error: cudaMallocHost failed for: get_nlop \n");
    for(idx = 0;idx < weight_size;idx++) pct.M_dnm[idx] = 0.0;

#else
    my_calloc (pct.weight, weight_size, rmg_double_t);
    my_calloc (pct.Bweight, weight_size, rmg_double_t);
    weight_size = pct.num_tot_proj * pct.num_tot_proj + 128;
    my_calloc (pct.M_dnm, weight_size, rmg_double_t);
    my_calloc (pct.M_qqq, weight_size, rmg_double_t);
#endif



    /*Make sure that ownership of ions is properly established
     * This conditional can be removed if it is found that claim_ions works reliably*/
    if (int_sum_all (pct.num_owned_ions, pct.grid_comm) != ct.num_ions)
        error_handler ("Problem with claimimg ions, %d ions claimed, but num_ions is %d",
                       int_sum_all (pct.num_owned_ions, pct.grid_comm), ct.num_ions);

    printf ("\n PE %d: Number of nonlocal ions is %d", pct.gridpe, pct.num_nonloc_ions);

    for (i = 0; i < pct.num_nonloc_ions; i++)
        printf (" %d", pct.nonloc_ions_list[i]);

    printf ("\n PE %d: Number of claimed ions is %d", pct.gridpe, pct.num_owned_ions);
    for (i = 0; i < pct.num_owned_ions; i++)
        printf (" %d", pct.owned_ions_list[i]);

    
    
    /*Memory for nonlocal projectors */
    if (pct.newsintR_local)
        my_free (pct.newsintR_local);
    if (pct.oldsintR_local)
        my_free (pct.oldsintR_local);
    
    my_calloc (pct.newsintR_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl,
               rmg_double_t);
    my_calloc (pct.oldsintR_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl,
               rmg_double_t);

#if !GAMMA_PT
    if (pct.newsintI_local)
        my_free (pct.newsintI_local);
    if (pct.oldsintI_local)
        my_free (pct.oldsintI_local);
    
    my_calloc (pct.newsintI_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl,
               rmg_double_t);
    my_calloc (pct.oldsintI_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl,
               rmg_double_t);
#endif



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
	Dprintf ("Owner of ion %d nlion %d is PE %d", ion, nlion, owner); 
        
	owned = 0;
        if (pct.gridpe == owner)
            owned = 1;
	Dprintf("owned set to %d for ion %d", owned, ion);



	/*Case when current PE owns an ion, in such a case we need to assemble
	 * list of other PEs that have overlap - those non-owning PEs will be sending data*/
	if (owned)
	{
	    Dprintf("ion %d is owned", ion);
	    
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

		Dprintf("Overlap condition for ion %d and PE %d is %d, map is %d, map2 is %d ", ion, pe, map || map2, map, map2); 

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

		    Dprintf("known_nonowner is %d and nonwner_index is %d for PE %d and ion %d", known_nonowner, nonwner_index, pe, ion);

		    /* Record pe index into communication list, if it has not been recorded yet*/
		    if (!known_nonowner)
		    {
			if (pct.num_owned_pe >= MAX_NONLOC_PROCS)
			    error_handler ("pct.num_nonloc_pes (%d) is too large (max: %d)",
				    pct.num_owned_pe, MAX_NONLOC_PROCS);

			pct.owned_pe_list[pct.num_owned_pe] = pe;

			pct.num_owned_ions_per_pe[pct.num_owned_pe] = 1;

			pct.list_owned_ions_per_pe[pct.num_owned_pe][0] = nlion;

			pct.num_owned_pe++;

			Dprintf("Recording previously unknown PE %d, pct.owned_pe_list[%d]=%d,  pct.num_owned_ions_per_pe[%d]=%d, pct.list_owned_ions_per_pe[%d][0]=%d,  pct.num_owned_pe = %d for ion %d", 
				pe, pct.num_owned_pe, pct.owned_pe_list[pct.num_owned_pe],pct.num_owned_pe, pct.num_owned_ions_per_pe[pct.num_owned_pe], pct.num_owned_pe, pct.list_owned_ions_per_pe[pct.num_owned_pe][0], pct.num_owned_pe, ion);

		    }

		    else
		    {
			/*Number of ions to communicate with current non-owner PE*/
			owned_ions_per_pe = pct.num_owned_ions_per_pe[nonwner_index];

			if (owned_ions_per_pe >= MAX_NONLOC_IONS)
			    error_handler
				("pct.num_owned_ions_per_pe too large (%d, MAX = %d ), for ion %d and PE %d",
				 pct.num_owned_ions_per_pe[nonwner_index], MAX_NONLOC_IONS, ion, pe);

			pct.list_owned_ions_per_pe[nonwner_index][owned_ions_per_pe] = nlion;
			pct.num_owned_ions_per_pe[nonwner_index]++;

			Dprintf("ion %d: PE %d is already known, pct.list_owned_ions_per_pe[%d][%d]=%d, pct.num_owned_ions_per_pe[%d]=%d",
				ion, pe, nonwner_index, owned_ions_per_pe, pct.list_owned_ions_per_pe[nonwner_index][owned_ions_per_pe],
				nonwner_index, pct.num_owned_ions_per_pe[nonwner_index]);
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

	    Dprintf("known_owner is %d and owner_index is %d for ion %d", known_owner, owner_index, ion);


	    if (!known_owner)
	    {
		if (pct.num_owners >= MAX_NONLOC_PROCS)
		    error_handler ("pct.num_owners (%d) is too large (max: %d)", pct.num_owners,
			    MAX_NONLOC_PROCS);

		pct.owners_list[pct.num_owners] = owner;

		pct.num_nonowned_ions_per_pe[pct.num_owners] = 1;

		pct.list_ions_per_owner[pct.num_owners][0] = nlion;

		pct.num_owners++;

		Dprintf("Owner of ion %d PE %d has not been previously recorded: pct.owners_list[%d]=%d, pct.num_nonowned_ions_per_pe[%d]=%d, pct.list_ions_per_owner[%d][0]=%d, pct.num_owners=%d", 
			ion, owner, pct.num_owners, pct.owners_list[pct.num_owners], pct.num_owners, pct.num_nonowned_ions_per_pe[pct.num_owners], pct.num_owners, pct.list_ions_per_owner[pct.num_owners][0], pct.num_owners); 
	    }

	    else
	    {

		nonowned_ions_per_pe = pct.num_nonowned_ions_per_pe[owner_index];

		if (nonowned_ions_per_pe >= MAX_NONLOC_IONS)
		    error_handler
			("pct.num_nonowned_ions_per_pe too large (%d, MAX = %d ), for ion %d and owner PE %d",
			 pct.num_nonowned_ions_per_pe[owner_index], MAX_NONLOC_IONS, ion, owner);

		pct.list_ions_per_owner[owner_index][nonowned_ions_per_pe] = nlion;
		pct.num_nonowned_ions_per_pe[owner_index]++;

		Dprintf("ion %d owned by %d: pct.list_ions_per_owner[%d][%d]=%d, pct.num_nonowned_ions_per_pe[%d]=%d", 
			ion, owner, owner_index, nonowned_ions_per_pe, pct.list_ions_per_owner[owner_index][nonowned_ions_per_pe], owner_index, pct.num_nonowned_ions_per_pe[owner_index]);
	    }

	}

#endif

    }/*for (nlion = 0; nlion < pct.num_nonloc_ions; nlion++)*/


    /* Release temporary memory */
    my_free (Aiz2);
    my_free (Aiy2);
    my_free (Aix2);
    my_free (Aiz);
    my_free (Aiy);
    my_free (Aix);

    my_free (pvec);

}                               /* end get_nlop */



/*Resets pct projector arrays, deallocates memory*/
static void reset_pct_arrays (int num_ions)
{

    int ion;
    for(ion = 0; ion < num_ions; ion++)
    {
        pct.idxptrlen[ion] = 0;

        if (pct.idxflag[ion])
            my_free (pct.idxflag[ion]);

        if (pct.nlindex[ion])
            my_free (pct.nlindex[ion]);
#if FDIFF_BETA
        if (pct.weight_derx[ion])
            my_free (pct.weight_derx[ion]);
        if (pct.weight_dery[ion])
            my_free (pct.weight_dery[ion]);
        if (pct.weight_derz[ion])
            my_free (pct.weight_derz[ion]);
#endif


#if !GAMMA_PT
        if (pct.phaseptr[ion])
            my_free (pct.phaseptr[ion]);
#endif
    }



    if (pct.weight != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.weight);
#else
        my_free (pct.weight);
#endif
    }
    if (pct.Bweight != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.Bweight);
#else
        my_free (pct.Bweight);
#endif
    }

    if (pct.M_dnm != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.M_dnm);
#else
        my_free (pct.M_dnm);
#endif
    }

    if (pct.M_qqq != NULL) {
#if GPU_ENABLED
        cudaFreeHost(pct.M_qqq);
#else
        my_free (pct.M_qqq);
#endif
    }


}

/******/
