/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"


static void reset_pct_arrays (int ion);

void get_nlop (void)
{

    int ion, idx, i, j, pe, pe_bin;
    int ix, iy, iz, ip, pcount;
    int *pvec, *dvec, *ivec;
    int ilow, jlow, klow, ihi, jhi, khi, map, icount;
    int alloc;
    int Aix[NX_GRID], Aiy[NY_GRID], Aiz[NZ_GRID];
    int icut, itmp, icenter;
    REAL vect[3];
    SPECIES *sp;
    ION *iptr;
    REAL nlxcstart, nlycstart, nlzcstart;


    /*Reset number of nonlocal ions */
    pct.num_nonloc_ions = 0;
    pct.num_nonloc_pes = 0;


    /* Grab some memory for temporary storage */
    alloc = ct.max_nlpoints;
    if (alloc < P0_BASIS)
	alloc = P0_BASIS;
    my_malloc (pvec, 2 * alloc, int);
    dvec = pvec + alloc;

    /* Loop over ions */
    for (ion = 0; ion < ct.num_ions; ion++)
    {

	/*Release memory and reset other parameters for given ion */
	reset_pct_arrays (ion);

	/* Generate ion pointer */
	iptr = &ct.ions[ion];


	/* Get species type */
	sp = &ct.sp[iptr->species];


	icenter = sp->nldim / 2;
	icut = (icenter + 1) * (icenter + 1);

	/* Set up Kleinman-Bylander normalization coefficients */
	pcount = 0;
	for (ip = 0; ip < sp->nbeta; ip++)
	{
	    switch (sp->llbeta[ip])
	    {
		case S_STATE:
		    pcount++;
		    break;

		case P_STATE:
		    pcount += 3;
		    break;

		case D_STATE:
		    pcount += 5;
		    break;
		default:
		    error_handler ("Angular momentum state not programmed");
	    }

	}

	pct.prj_per_ion[ion] = pcount;


	/* Determine mapping indices or even if a mapping exists */
	map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
		sp->nldim, PX0_GRID, PY0_GRID, PZ0_GRID,
		ct.psi_nxgrid, ct.psi_nygrid, ct.psi_nzgrid,
		&iptr->nlxcstart, &iptr->nlycstart, &iptr->nlzcstart);

	/*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
	/*xtal vector between ion and left bottom corner of the box */
	vect[0] = iptr->xtal[0] - iptr->nlxcstart;
	vect[1] = iptr->xtal[1] - iptr->nlycstart;
	vect[2] = iptr->xtal[2] - iptr->nlzcstart;

	/*Substract vector between left bottom corner of the box and center of the box */
	vect[0] -= (sp->nldim / 2) / (REAL) ct.psi_nxgrid;
	vect[1] -= (sp->nldim / 2) / (REAL) ct.psi_nygrid;
	vect[2] -= (sp->nldim / 2) / (REAL) ct.psi_nzgrid;

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

			    /*   if(icut >= itmp) { */

			    pvec[icount] =
				PY0_GRID * PZ0_GRID * (Aix[ix] % PX0_GRID) +
				PZ0_GRID * (Aiy[iy] % PY0_GRID) + (Aiz[iz] % PZ0_GRID);

			    dvec[idx] = TRUE;

			    icount++;

			    /*  } */

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



	    if ((icount * pct.prj_per_ion[ion]))
	    {
		size_t weight_size = pct.prj_per_ion[ion] * icount + 128;

		my_calloc (pct.weight[ion], weight_size, REAL);


#if FDIFF_BETA
		my_calloc (pct.weight_derx[ion], weight_size, REAL);
		my_calloc (pct.weight_dery[ion], weight_size, REAL);
		my_calloc (pct.weight_derz[ion], weight_size, REAL);
#endif
	    }
	    else
	    {
		pct.weight[ion] = NULL;
#if FDIFF_BETA
		pct.weight_derx[ion] = NULL;
		pct.weight_dery[ion] = NULL;
		pct.weight_derz[ion] = NULL;
#endif
	    }


#if !GAMMA_PT

	    /* Allocate memory for the phase array */
	    if ((icount * pct.prj_per_ion[ion]))
		my_calloc (pct.phaseptr[ion], 2 * icount * ct.num_kpts + 128, REAL);
	    else
		pct.phaseptr[ion] = NULL;

	    get_phase (iptr, pct.phaseptr[ion], ip, icount, dvec);
#endif

	}                       /* end if (map) */


#if 1
	/*Add ion into list of nonlocal ions if it has overlap with given processor */
	if (icount)
	{
	    if (pct.num_nonloc_ions >= MAX_NONLOC_IONS) error_handler ("Too many nonlocal ions, pct.nonloc_ions_list will overflow");
	    pct.nonloc_ions_list[pct.num_nonloc_ions] = ion;
	    pct.num_nonloc_ions++;
	}
#endif

    }                           /* end for (ion = 0; ion < ct.num_ions; ion++) */

    printf("\n PE %d: Number of nonlocal ions is %d", pct.gridpe, pct.num_nonloc_ions);

    for (i=0; i<pct.num_nonloc_ions; i++)
	printf(" %d", pct.nonloc_ions_list[i]);

    /*Memory for nonlocal projectors*/
    if (pct.newsintR_local) my_free (pct.newsintR_local);
    if (pct.oldsintR_local) my_free (pct.oldsintR_local);
    my_calloc (pct.newsintR_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl, REAL);
    my_calloc (pct.oldsintR_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl, REAL);
    
#if !GAMMA_PT
    if (pct.newsintI_local) my_free (pct.newsintI_local);
    if (pct.oldsintI_local) my_free (pct.oldsintI_local);
    my_calloc (pct.newsintI_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl, REAL)
    my_calloc (pct.oldsintI_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl, REAL);
#endif
    


    /* Loop over all nonlocal ions to obtain the lists necessary for communication*/
    for (i = 0; i < pct.num_nonloc_ions; i++)
    {
	ion = pct.nonloc_ions_list[i];

	//printf("\n\n Considering %dth ion number %d", i, ion); 


	/* Generate ion pointer */
	iptr = &ct.ions[ion];

	/* Get species type */
	sp = &ct.sp[iptr->species];

	icenter = sp->nldim / 2;
	icut = (icenter + 1) * (icenter + 1);

	/*Loop over all processors in caculation to determine if it overlaps with current ion*/
	for (pe=0; pe < NPES; pe++)
	{

	    /* Skip current processor*/
	    if (pe == pct.gridpe) continue;


	    /************* This is a copy of what is done above to find whether an ion has overlap with a given processor
	     ************ The idea is follow the same procedure (for any processor) and determine 
	     variable   icount. If it is non-zero, we have an ovelap*/

	    /* Determine mapping indices or even if a mapping exists */
	    map = get_index (pe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
		    sp->nldim, PX0_GRID, PY0_GRID, PZ0_GRID,
		    ct.psi_nxgrid, ct.psi_nygrid, ct.psi_nzgrid,
		    &nlxcstart, &nlycstart, &nlzcstart);

	    /*Find nlcdrs, vector that gives shift of ion from center of its ionic box */
	    /*xtal vector between ion and left bottom corner of the box */
	    vect[0] = iptr->xtal[0] - nlxcstart;
	    vect[1] = iptr->xtal[1] - nlycstart;
	    vect[2] = iptr->xtal[2] - nlzcstart;

	    /*Substract vector between left bottom corner of the box and center of the box */
	    vect[0] -= (sp->nldim / 2) / (REAL) ct.psi_nxgrid;
	    vect[1] -= (sp->nldim / 2) / (REAL) ct.psi_nygrid;
	    vect[2] -= (sp->nldim / 2) / (REAL) ct.psi_nzgrid;

	    /*The vector we are looking for should be */
	    to_cartesian (vect, iptr->nlcrds);



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

				/*   if(icut >= itmp) { */

				pvec[icount] =
				    PY0_GRID * PZ0_GRID * (Aix[ix] % PX0_GRID) +
				    PZ0_GRID * (Aiy[iy] % PY0_GRID) + (Aiz[iz] % PZ0_GRID);

				dvec[idx] = TRUE;

				icount++;

				/*  } */

			    }
			    idx++;
			}
		    }
		}

		/******* At this point variable icount is determined and what follows is jus recording the information
		 * into the list */

		/*If this condition is satisfied then processor "pe" has overlap with the current nonlocal ion*/
		if (icount)
		{
		    //printf("\n PE %d has overlap!", pe); 

		    /*Loop over list of processors with overlap to see if this a processor that is already recorded*/
		    pe_bin = -1;
		    for (j=0; j< pct.num_nonloc_pes; j++)
			if (pe == pct.nonloc_pe_list[j]) pe_bin = j;

		    /* The case when this is a processor that has not been recorded yet*/
		    if (pe_bin == -1)
		    {
			//printf("\n PE %d has NOT been recorded yet!", pe); 
			if (pct.num_nonloc_pes >= MAX_NONLOC_PROCS) error_handler("Too many PEs to fit into pct.nonloc_pe_list, max is %d", MAX_NONLOC_PROCS);
			pct.nonloc_pe_list[pct.num_nonloc_pes] = pe;

			/*Record atom*/
			pct.nonloc_atom_list_per_pe[pct.num_nonloc_pes][0] = ion;
			//printf("\n Recording ion %d into pct.nonloc_atom_list_per_pe[%d][0]", ion, pct.num_nonloc_pes); 


			/*This keeps count for number of ions that need to be communicated for a given processor*/
			pct.nonloc_atom_count_per_pe[pct.num_nonloc_pes] = 1;

			/*This keeps count for number of processors to communicate with*/
			pct.num_nonloc_pes ++;

		    }

		    else
		    {
			//printf("\n PE %d has been recorded yet!", pe); 
			/*if we are here, processor was already recorded, we need to additional ion to the proper list*/
			pct.nonloc_atom_list_per_pe[pe_bin][pct.nonloc_atom_count_per_pe[pe_bin]] = ion;

			//printf("\n Recording ion %d into pct.nonloc_atom_list_per_pe[%d][%d]", ion, pe_bin, pct.nonloc_atom_count_per_pe[pe_bin]); 
			pct.nonloc_atom_count_per_pe[pe_bin] ++;

		    }

		}


	    }

	}

    }

#if 0
    printf("\n PE %d: Number of nonloc ions %d, number of PEs to communicate with is %d", pct.gridpe, pct.num_nonloc_ions, pct.num_nonloc_pes);

    for (i=0; i<pct.num_nonloc_pes; i++)
    {
	printf("\n Atoms to communicate about with PE %d:", pct.nonloc_pe_list[i]);

	for (j=0; j < pct.nonloc_atom_count_per_pe[i]; j++)
	    printf("  %d", pct.nonloc_atom_list_per_pe[i][j]);

    }
#endif



    /* Release temporary memory */
    my_free (pvec);

}                               /* end get_nlop */








/*Resets pct projector arrays, deallocates memory*/
static void reset_pct_arrays (int ion)
{

    pct.idxptrlen[ion] = 0;

    if (pct.idxflag[ion])
	my_free (pct.idxflag[ion]);

    if (pct.nlindex[ion])
	my_free (pct.nlindex[ion]);


    if (pct.weight[ion])
	my_free (pct.weight[ion]);

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


/******/
