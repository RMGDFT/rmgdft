/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"

/*Call to this function needs to be preceeded by get_QI, since we use pct.Qidxptrlen,
 * which is setup in that function*/

static void reset_pct_arrays (int ion);

void get_nlop (void)
{

    int ion, idx, i, j, pe, pe_bin;
    int ix, iy, iz, ip, prj_per_ion;
    int *pvec, *dvec, *ivec;
    int ilow, jlow, klow, ihi, jhi, khi, map, map2, flag, icount;
    int alloc, pe_index;
    int Aix[NX_GRID], Aiy[NY_GRID], Aiz[NZ_GRID];
    int Aix2[FNX_GRID], Aiy2[FNY_GRID], Aiz2[FNZ_GRID];
    int icut, itmp, icenter;
    REAL vect[3];
    SPECIES *sp;
    ION *iptr;
    REAL nlxcstart, nlycstart, nlzcstart;


    /*Reset number of nonlocal ions */
    pct.num_nonloc_ions = 0;
    pct.num_owned_ions = 0;
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

        prj_per_ion = sp->nh;


        icenter = sp->nldim / 2;
        icut = (icenter + 1) * (icenter + 1);


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

                            if (icut >= itmp) { 

                                pvec[icount] =
                                    PY0_GRID * PZ0_GRID * (Aix[ix] % PX0_GRID) +
                                    PZ0_GRID * (Aiy[iy] % PY0_GRID) + (Aiz[iz] % PZ0_GRID);

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
                size_t weight_size = prj_per_ion * icount + 128;

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
            if ((icount * prj_per_ion))
                my_calloc (pct.phaseptr[ion], 2 * icount * ct.num_kpts + 128, REAL);
            else
                pct.phaseptr[ion] = NULL;

            get_phase (iptr, pct.phaseptr[ion], ip, icount, dvec);
#endif

        }                       /* end if (map) */


        /*Add ion into list of nonlocal ions if it has overlap with given processor */
        if (icount || pct.Qidxptrlen[ion])
        {
            
            if (pct.num_nonloc_ions >= MAX_NONLOC_IONS) 
                error_handler ("Too many nonlocal ions (pct.num_nonloc_ions = %d MAX = %d ), pct.nonloc_ions_list will overflow",pct.num_nonloc_ions, MAX_NONLOC_IONS);

            pct.nonloc_ions_list[pct.num_nonloc_ions] = ion;
            pct.num_nonloc_ions++;

            /*See if this processor owns current ion*/
            if (claim_ion(pct.gridpe, iptr, PX0_GRID, PY0_GRID, PZ0_GRID, ct.psi_nxgrid, ct.psi_nygrid, ct.psi_nzgrid))
            {
                if (pct.num_owned_ions >= MAX_NONLOC_IONS) 
                    error_handler ("Too many owned ions, pct.owned_ions_list will overflow");

                pct.owned_ions_list[pct.num_owned_ions] = ion;
                pct.num_owned_ions++;
                //dprintf("I own ion %d", ion);
            }



        }

    }                           /* end for (ion = 0; ion < ct.num_ions; ion++) */

    /*Make sure that ownership of ions is properly established
     * This conditional can be removed if it is found that claim_ions works reliably*/
    if (real_sum_all(pct.num_owned_ions, pct.grid_comm) != ct.num_ions)
        error_handler("Problem with claimimg ions, %d ions claimed, but num_ions is %d", real_sum_all(pct.num_owned_ions, pct.grid_comm), ct.num_ions);

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
    my_calloc (pct.newsintI_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl, REAL);
    my_calloc (pct.oldsintI_local, ct.num_kpts * pct.num_nonloc_ions * ct.num_states * ct.max_nl, REAL);
#endif



#if 0
    
    /* Loop over all ions to obtain the lists necessary for communication*/
    for (j = 0; j < pct.num_nonloc_ions; j++)
    {

        ion = pct.nonloc_ions_list[j];
        //printf("\n\n Considering ion number %d", ion); 

        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];


        /*Loop over all processors in caculation to determine if it overlaps with current ion*/
        /* This should only loop over processors in a grid */
        for (pe=0; pe < NPES; pe++)
        {
            if (pe == pct.gridpe)
                continue;

            /* Determine if ion has overlap with a given PE becasue of beta functions */
            map = test_overlap(pe, iptr, Aix, Aiy, Aiz, sp->nldim,
                    PX0_GRID, PY0_GRID, PZ0_GRID,
                    NX_GRID, NY_GRID, NZ_GRID);
        
            /* Determine if ion has overlap with a given PE becasue of Q function */
            map2 = test_overlap(pe, iptr, Aix2, Aiy2, Aiz2,
                         sp->qdim, FPX0_GRID, FPY0_GRID, FPZ0_GRID,
                         FNX_GRID, FNY_GRID, FNZ_GRID);

            printf("\n test_overlap returned %d and %d", map, map2);
                         

            if (map || map2)
            {
            
                /*Loop over list of processors that we already have to see if "pe" has alredy been recorded*/
                flag = 0;
                for (i=0; i < pct.num_nonloc_pes; i++)
                {
                    if (pct.nonloc_pe_list[i] == pe)
                    {
                          flag ++;
                          pe_index = i;
                          break;
                    }
                }

                /*Record pe index, if not known yet*/ 
                if (!flag)
                {
                    if (pct.num_nonloc_pes >= MAX_NONLOC_PROCS) 
                        error_handler("pct.num_nonloc_pes (%d) is too large (max: %d)", pct.nonloc_pe_list_ions, MAX_NONLOC_PROCS);

                    pct.nonloc_pe_list[pct.num_nonloc_pes] = pe;
                        
                    pct.nonloc_pe_list_ions[pct.num_nonloc_pes][0] = ion;

                    pct.nonloc_pe_num_ions[pct.num_nonloc_pes] = 1;
                    
                    pct.num_nonloc_pes ++;
                }

                else
                {
                    if (pct.nonloc_pe_num_ions[pe_index] >= MAX_NONLOC_IONS) 
                        error_handler ("pct.nonloc_pe_num_ions too large (%d, MAX = %d ), for ion %d and PE %d",pct.nonloc_pe_num_ions[pe_index], MAX_NONLOC_IONS, ion, pe);
                        
                    pct.nonloc_pe_list_ions[pe_index][pct.nonloc_pe_num_ions[pe_index]] = ion;
                    pct.nonloc_pe_num_ions[pe_index] ++;
                }

            }


        }

    }

#endif

#if 0
    printf("\n PE %d: Number of nonloc ions %d, number of PEs to communicate with is %d", pct.gridpe, pct.num_nonloc_ions, pct.num_nonloc_pes);

    for (i=0; i<pct.num_nonloc_pes; i++)
    {
        printf("\n Atoms to communicate about with PE %d:", pct.nonloc_pe_list[i]);

        for (j=0; j < pct.nonloc_pe_num_ions[i]; j++)
            printf("  %d",  pct.nonloc_pe_list_ions[i][j]);

    }

    dprintf("\n PE %d: Number of nonloc ions %d", pct.gridpe, pct.num_nonloc_ions);
    for (i=0; i<pct.num_nonloc_ions; i++)
    {
        dprintf("\n PE %d, nonlocal ion %d: %d", pct.gridpe, i, pct.nonloc_ions_list[i]);
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
