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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"

/*Initializes atomic wavefunctions (multiplies radial part read from pseudpotential file 
 * with spherical harmonics)*/


#if 0
void mulliken (STATE * states)
{

    int ion, idx, species, state, dim, icenter, i, j, icut, map, ilow, ihi, jlow, jhi, klow, khi;
    int ix, iy, iz, grid_index, count, l_state, m, index, oindex, st, nindex;
    double *awave, *overlap, xcstart, ycstart, zcstart, xc, yc, zc, x[3], r, cx[3], radial_wave,
        dot_product;
    double a, b, rdiff, rleft, rright, *norm_factor, max1, max2, max3;
    int basis =get_P0_BASIS();
    int tot_atomic_states;
    int Aix[NX_GRID], Aiy[NY_GRID], Aiz[NZ_GRID];
    ION *iptr;
    SPECIES *sp;
    int rindex, awave_max, ii;


    /*First, let us calculate total number of atomic waves (sum over all ions */
    tot_atomic_states = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {
        iptr = &ct.ions[ion];
        species = iptr->species;
        sp = &ct.sp[species];

        tot_atomic_states += sp->sum_atomic_waves;
    }

    /*Allocate memory to hold overlap between atomic waves and wavefunctions */
    my_malloc (overlap, tot_atomic_states * ct.num_states, double);

    awave_max = 16 * basis;
    my_malloc (awave, awave_max, double);

    my_malloc (norm_factor, tot_atomic_states, double);





    oindex = 0;
    nindex = 0;
    for (ion = 0; ion < ct.num_ions; ion++)
    {

        /*Reset awave array */
        for (idx = 0; idx < awave_max; idx++)
            awave[idx] = 0.0;

        iptr = &ct.ions[ion];
        species = iptr->species;
        sp = &ct.sp[species];
        dim = sp->mill_dim;

        icenter = sp->mill_dim / 2;
        icut = (icenter + 1) * (icenter + 1);

        a = sp->aa;
        b = sp->bb;

        /*if (pct.gridpe == 0)
           printf("\n Ion:%d a:%f b:%f  num_atomic_waves:%d  sum_atomic_waves:%d rg_points:%d", ion, a, b, sp->num_atomic_waves, sp->sum_atomic_waves, sp->rg_points); */


        count = 0;

        map = get_index (pct.gridpe, iptr, Aix, Aiy, Aiz, &ilow, &ihi, &jlow, &jhi, &klow, &khi,
                         dim, pct.PX0_GRID, pct.PY0_GRID, pct.PZ0_GRID,
                         get_NX_GRID(), get_NY_GRID(), get_NZ_GRID(), &xcstart, &ycstart, &zcstart);



        if (map)
        {

            xc = xcstart;
            for (ix = 0; ix < dim; ix++)
            {

                yc = ycstart;
                for (iy = 0; iy < dim; iy++)
                {

                    zc = zcstart;
                    for (iz = 0; iz < dim; iz++)
                    {

                        if (((Aix[ix] >= ilow) && (Aix[ix] <= ihi)) &&
                            ((Aiy[iy] >= jlow) && (Aiy[iy] <= jhi)) &&
                            ((Aiz[iz] >= klow) && (Aiz[iz] <= khi)))
                        {

                            /*Maybe spherical cutoff here ?? */

                            grid_index = pct.PY0_GRID * pct.PZ0_GRID * (Aix[ix] % pct.PX0_GRID) +
                                pct.PZ0_GRID * (Aiy[iy] % pct.PY0_GRID) + (Aiz[iz] % PZ0_GRID);

                            x[0] = xc - iptr->xtal[0];
                            x[1] = yc - iptr->xtal[1];
                            x[2] = zc - iptr->xtal[2];

                            r = metric (x);
                            to_cartesian (x, cx);



                            /*Calculate index in radial grid that give r corresponds to */
                            rindex = (int) ((log ((r + a) / a) / b) + 1.0);

                            /*Index in array read from psedopotential is actually shifted by 2 */
                            rindex -= 2;



                            /*Just a check that rindex was determined correctly */
                            if (rindex >= 0)
                            {

                                if (r < sp->r[rindex] || r > sp->r[rindex + 1])
                                {
                                    rleft = sp->aa * (exp ((rindex + 2) * sp->bb) - 1.0);
                                    rright = sp->aa * (exp ((rindex + 3) * sp->bb) - 1.0);

                                    printf
                                        ("\nPE%d: Possibly wrong index:  ion %d, r  %f sp->r[rindex] %f sp->r[rindex+1] %f rleft %f rright %f",
                                         pct.gridpe, ion, r, sp->r[rindex], sp->r[rindex + 1],
                                         rleft, rright);
                                    //error_handler("wrong rindex");
                                }
                            }



                            count = 0;
                            for (i = 0; i < sp->num_atomic_waves; i++)
                            {

                                /*Get radial part of atomic wave function using linear interpolation */
                                if (rindex >= 0)
                                {
                                    rdiff = r - sp->r[rindex];

                                    radial_wave =
                                        (sp->atomic_wave[i][rindex + 1] -
                                         sp->atomic_wave[i][rindex]) * rdiff / (sp->r[rindex + 1] -
                                                                                sp->r[rindex]) +
                                        sp->atomic_wave[i][rindex];
                                }

                                /*For very small distances just use the smallest known value of atomic wave function */
                                else
                                    radial_wave = sp->atomic_wave[i][0];


                                /*Check only */
                                /*if (rindex >= 0)
                                   printf("\n r:%f r[rindex]:%f r[rindex+1]:%f, wave: %f %f %f",
                                   r, sp->r[rindex], sp->r[rindex+1], radial_wave, sp->atomic_wave[i][rindex], sp->atomic_wave[i][rindex+1]); */


                                l_state = sp->lstate_atomic_wave[i];


                                /*Starting index for ylm function
                                 * The functions are indexed (ylm)
                                 * 0:s, 1:px, 2:py, 3:pz, 4:dxx, etc.*/
                                index = l_state * l_state;


                                /*Calculate spherical harmonics and add them 
                                 * to radial part get full atomic wave functions*/
                                for (m = index; m < index + 2 * l_state + 1; m++)
                                {

                                    //awave[count*basis+grid_index] = pow(r,l_state) * radial_wave * ylm(m, cx);
                                    awave[count * basis + grid_index] = radial_wave * ylm (m, cx);
                                    count++;
                                }




                            }   /*end for (i=0; i<sp->num_atomic_waves; i++) */



                        }


                        zc += get_hzgrid();

                    }           /* end for iz */
                    yc += get_hygrid();

                }               /* end for iy */
                xc += get_hxgrid();

            }                   /* end for ix */

        }                       /*end if map */



        /*This should be satisfied, otherwise we are doing something wrong */
        if (count && (sp->sum_atomic_waves != count))
        {
            printf ("\n PE %d: Ion %d, species %d, sum_state_atomic_wave is %d, count is %d",
                    pct.gridpe, ion, species, sp->sum_atomic_waves, count);

            error_handler ("count should be equal to sum_state_atomic_wave !!!");
        }






        /*This checks orthogonality and normalization of atomic wave functions */
        for (i = 0; i < sp->sum_atomic_waves; i++)
        {


#if 0
            /* Apply non-local operator to psi and store in work2 */
            app_nl_eig (&awave[i * basis], NULL, work2, NULL, sp->istate, FALSE, kidx, 0);

            /* Apply the S operator acting on psi and store in work3 */
            app_ns_eig (&awave[i * basis], NULL, work3, NULL, sp->istate, sp->kidx, 0);

            for (idx = 0; idx < basis; idx++)
                awave[i * basis + idx] = work3;
#endif


            /*Check normalization of atomic waves */
            dot_product = 0.0;

            for (idx = 0; idx < basis; idx++)
                dot_product += awave[i * basis + idx] * awave[i * basis + idx];


            dot_product = get_vel() * real_sum_all (dot_product, pct.grid_comm);

            norm_factor[nindex] = dot_product;

            if (pct.gridpe == 0)
                printf ("\n Ion %d:   atomic state %d is normalized to %e", ion, i, dot_product);

            /*for(idx=0; idx < basis; idx++) 
               awave[i*basis+idx] /= norm_factor; */

            if (i >= 1)
                for (j = i - 1; j >= 0; j--)
                {
                    /*Check orthogonalization of atomic waves */
                    dot_product = 0.0;

                    for (idx = 0; idx < basis; idx++)
                        dot_product += awave[i * basis + idx] * awave[j * basis + idx];

                    dot_product = get_vel() * real_sum_all (dot_product, pct.grid_comm);

                    if (pct.gridpe == 0)
                        printf ("\n Ion %d:  Dot product between  atomic states %d and %d  is %e",
                                ion, i, j, dot_product);
                }


            nindex++;
        }                       /*end for(i=0; i<sp->sum_atomic_waves; i++) */




        for (i = 0; i < sp->sum_atomic_waves; i++)
        {

            if (pct.gridpe == 0)
                printf ("\n");

            for (state = 0; state < ct.num_states; state++)
            {

                /*Now we can do dot product of wave functions with atomic wave functions */
                dot_product = 0.0;

                for (idx = 0; idx < basis; idx++)
                    dot_product += awave[i * basis + idx] * states[state].psiR[idx];



                overlap[oindex] = get_vel() * real_sum_all (dot_product, pct.grid_comm);

                /*if (pct.gridpe == 0)
                   printf("\n Ion %d: Dot product between state %d and atomic state %d is %e", ion, state, i, overlap[oindex]); */

                oindex++;

            }                   /*end for (state=0; state < ct.num_states; state++) */


        }                       /*end for(i=0; i<sp->sum_atomic_waves; i++) */







    }                           /*end for (ion = 0; ion < ct.num_ions; ion++) */


    /*Now let us print table with results */
    if (pct.gridpe == 0)
    {
        printf ("\n\n");



        /*Header first */
        printf ("\n%s     ", "Atom Orbitals");

        for (ion = 0; ion < ct.num_ions; ion++)
        {

            iptr = &ct.ions[ion];
            species = iptr->species;
            sp = &ct.sp[species];

            for (i = 0; i < sp->sum_atomic_waves; i++)
                printf ("%s% -4d  ", "Ion", ion);
        }



        /*Second line */
        printf ("\n%s   ", "------------");

        for (ion = 0; ion < ct.num_ions; ion++)
        {

            iptr = &ct.ions[ion];
            species = iptr->species;
            sp = &ct.sp[species];

            for (i = 0; i < sp->num_atomic_waves; i++)
            {
                /*Quantum number l */

                l_state = sp->lstate_atomic_wave[i];

                for (m = 0; m < 2 * l_state + 1; m++)
                {
                    switch (l_state)
                    {
                        /*s state */
                    case 0:
                        printf ("%7s  ", "s");
                        break;


                        /*px, py and pz states */
                    case 1:
                        switch (m)
                        {
                        case 0:
                            printf ("%7s  ", "px");
                            break;

                        case 1:
                            printf ("%7s  ", "pz");
                            break;

                        case 2:
                            printf ("%7s  ", "py");
                            break;

                        default:
                            error_handler (" m should only be from 0 to 2 for p-state");
                        }       /*end switch(m) */
                        break;


                        /*d states */
                    case 2:
                        switch (m)
                        {
                        case 0:
                            printf ("%7s  ", "dxy");
                            break;

                        case 1:
                            printf ("%7s  ", "dxz");
                            break;

                        case 2:
                            printf ("%7s  ", "dzz");
                            break;

                        case 3:
                            printf ("%7s  ", "dyz");
                            break;


                        case 4:
                            printf ("%7s  ", "dxx-yy");
                            break;

                        default:
                            error_handler (" m should only be from 0 to 4 for d-state");

                        }       /*end switch(m) */
                        break;


                        /*Commenting out due to space constraints, Fz(xx-yy) takes 9 characters,
                         * seems too much since F is rarely used*/
                        /*F-state */
#if 0
                    case 2:
                        switch (m)
                        {
                        case 0:
                            printf ("%7s  ", "Fxxx");
                            break;

                        case 1:
                            printf ("%7s  ", "Fyyy");
                            break;

                        case 2:
                            printf ("%7s  ", "Fxyz");
                            break;

                        case 3:
                            printf ("%7s  ", "Fzzz");
                            break;


                        case 4:
                            printf ("%7s  ", "Fz(xx-yy)");
                            break;

                        case 5:
                            printf ("%7s  ", "Fy(zz-xx)");
                            break;


                        case 6:
                            printf ("%7s  ", "Fx(yy-zz)");
                            break;

                        default:
                            error_handler (" m should only be from 0 to 6 for f-state");

                        }       /*end switch(m) */
                        break;
#endif

                    default:
                        error_handler ("Higher l not programmed");

                    }           /*switch (l_state) */


                }               /*end for(m=0; m< 2*l_state; m++) */

            }                   /*end for(i=0; i<sp->num_atomic_waves; i++) */

        }                       /*end for(ion=0; ion<ct.num_ions; ion++) */





        /*Third line */
        printf ("\n%s     ", "Mol Orbitals");

        /*Print norm for each state */
        for (i = 0; i < tot_atomic_states; i++)
            printf (" N:% -3.2f ", norm_factor[i]);


        /*Header should now be completed, now let us write results of overlap */


        for (st = 0; st < ct.num_states; st++)
        {

            /*First state number and its occupancy */
            printf ("\n% 4d [% 3.2f]     ", st, ct.kp[0].kstate[st].occupation[0]);

            max1 = 0.0;
            max2 = 0.0;
            max3 = 0.0;

            for (i = 0; i < tot_atomic_states; i++)
            {

                if (fabs (overlap[st + i * ct.num_states]) > fabs (max1))
                {
                    max3 = max2;
                    max2 = max1;
                    max1 = overlap[st + i * ct.num_states];
                }
                else
                {
                    if (fabs (overlap[st + i * ct.num_states]) > fabs (max2))
                    {
                        max3 = max2;
                        max2 = overlap[st + i * ct.num_states];
                    }
                    else
                    {
                        if (fabs (overlap[st + i * ct.num_states]) > fabs (max3))
                            max3 = overlap[st + i * ct.num_states];
                    }

                }



                printf (" % 5.4f ", overlap[st + i * ct.num_states]);

            }



            /*Print 3 maxima for each state */
            printf ("    Maxima: % 5.4f % 5.4f % 5.4f", max1, max2, max3);




        }


        printf ("\n\n\n");






        /*Second printout, just print important overlaps */
        printf ("\n ----------- Overlaps over 0.1 ------------");


        for (st = 0; st < ct.num_states; st++)
        {

            /*First state number and its occupancy */
            printf ("\n% 4d [% 3.2f]     ", st, ct.kp[0].kstate[st].occupation[0]);

            for (i = 0; i < tot_atomic_states; i++)
            {

                if (fabs (overlap[st + i * ct.num_states]) > 0.1)
                {
                    /*We need to figure out ion that corresponds to current i */
                    ii = -1;
                    ion = 0;
                    do
                    {
                        ii += ct.sp[ct.ions[ion].species].sum_atomic_waves;
                        ion++;
                    }
                    while (ii < i);

                    /*ion overcounts by one */
                    ion--;

                    /*ii will be index of atomic orbital in ion "ion" */
                    ii = ii - ct.sp[ct.ions[ion].species].sum_atomic_waves + 1;
                    ii = i - ii;




                    printf (" Overlap: % 5.4f  Ion %d[%s] State %s   ",
                            overlap[st + i * ct.num_states], ion,
                            ct.sp[ct.ions[ion].species].pseudo_symbol,
                            ct.sp[ct.ions[ion].species].atomic_wave_symbol[ii]);

                }

            }

        }
















    }                           /*end if (pct.gridpe == 0) */


    my_free (awave);
    my_free (overlap);
    my_free (norm_factor);


}
#endif
