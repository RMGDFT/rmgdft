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

#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"

template void LcaoGetPsi(State<double> *);
template void LcaoGetPsi(State<std::complex<double> > *);

template <typename StateType>
void LcaoGetPsi (State<StateType> * states)
{

    int ion, idx, ip, l, m, state_count, st, P0_BASIS;
    int PX0_GRID, PY0_GRID, PZ0_GRID;
    int PX_OFFSET, PY_OFFSET, PZ_OFFSET;

    SPECIES *sp;
    ION *iptr;
    long idum;
    double coeff;

    PX0_GRID = get_PX0_GRID();
    PY0_GRID = get_PY0_GRID();
    PZ0_GRID = get_PZ0_GRID();
    PX_OFFSET = get_PX_OFFSET();
    PY_OFFSET = get_PY_OFFSET();
    PZ_OFFSET = get_PZ_OFFSET();


    /*The first index is due to k-point*/
    P0_BASIS = get_P0_BASIS();

    for (st = 0; st < ct.num_states; st++)
    {
        for (idx = 0; idx < P0_BASIS; idx++)
        {
            states[st].psi[idx] = 0.0;
        }
    }

    idum = 1314; 
    rand0 (&idum);

    /* Loop over ions */
    state_count = 0;

    for (ion = 0; ion < ct.num_ions; ion++)
    {
        /* Generate ion pointer */
        iptr = &ct.ions[ion];

        /* Get species type */
        sp = &ct.sp[iptr->species];

        /*Make sure that the wavefunctions have been read*/
        if (!sp->num_atomic_waves) {
            rmg_printf("No initial wavefunctions for ion %d, most likely the PP file does not have them", ion);
            rmg_error_handler(__FILE__,__LINE__,"Terminating.");
        }

        /*Loop over atomic wavefunctions for given ion*/
        for (ip = 0; ip < sp->num_atomic_waves; ip++)
        {
            l = sp->atomic_wave_l[ip];

            /*Loop over all m values for given l and get wavefunctions */
            for (m=0; m < 2*l+1; m++)
            {

                state_count++;

            }

        }

    }

    if(state_count <= ct.num_states)
    {
        coeff = 1.0;
        st = 0;
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            /* Generate ion pointer */
            iptr = &ct.ions[ion];

            /* Get species type */
            sp = &ct.sp[iptr->species];

            /*Make sure that the wavefunctions have been read*/
            if (!sp->num_atomic_waves) {
                rmg_printf("No initial wavefunctions for ion %d, most likely the PP file does not have them", ion);
                rmg_error_handler(__FILE__,__LINE__,"Terminating.");
            }

            /*Loop over atomic wavefunctions for given ion*/
            for (ip = 0; ip < sp->num_atomic_waves; ip++)
            {
                l = sp->atomic_wave_l[ip];

                /*Loop over all m values for given l and get wavefunctions */
                for (m=0; m < 2*l+1; m++)
                {
                    LcaoGetAwave(states[st].psi, iptr, ip, l, m, coeff);
                    st++;
                }
            }

        }
    }
    else
    {

        //StateType *aidum = new StateType[ct.num_states];
        //StateType *apsi = new StateType [P0_BASIS];
        long *aidum = new long[ct.num_states];
        double *apsi = new double [P0_BASIS];
        for(int st = 0;st < ct.num_states;st++) {
            aidum[st] = idum = st + 3314;
        }

        coeff = 1.0;
        for (ion = 0; ion < ct.num_ions; ion++)
        {
            /* Generate ion pointer */
            iptr = &ct.ions[ion];

            /* Get species type */
            sp = &ct.sp[iptr->species];

            /*Make sure that the wavefunctions have been read*/
            if (!sp->num_atomic_waves) {
                rmg_printf("No initial wavefunctions for ion %d, most likely the PP file does not have them", ion);
                rmg_error_handler(__FILE__,__LINE__,"Terminating.");
            }

            /*Loop over atomic wavefunctions for given ion*/
            for (ip = 0; ip < sp->num_atomic_waves; ip++)
            {
                l = sp->atomic_wave_l[ip];

                /*Loop over all m values for given l and get wavefunctions */
                for (m=0; m < 2*l+1; m++)
                {
                    for(idx = 0;idx < P0_BASIS;idx++)  apsi[idx] = 0.0;
                    LcaoGetAwave(apsi, iptr, ip, l, m, coeff);
                    for(int st = 0;st < ct.num_states;st++) 
                    {
                        double tem1 = rand0(&aidum[st]);
                        for(idx = 0;idx < P0_BASIS;idx++) 
                        {
                            states[st].psi[idx] += tem1 * apsi[idx];
                        }

                    }
                }

            }
        }

        delete [] apsi;
        delete [] aidum;

    }

    /*Initialize any additional states to random start*/
    if ( ct.num_states > state_count)
    {
        int ix, iy, iz;
        int xoff, yoff, zoff;
        State<StateType> *state_p;

        double *xrand = new double[2 * get_NX_GRID()];
        double *yrand = new double[2 * get_NY_GRID()];
        double *zrand = new double[2 * get_NZ_GRID()];

        pe2xyz (pct.gridpe, &ix, &iy, &iz);
        xoff = PX_OFFSET;
        yoff = PY_OFFSET;
        zoff = PZ_OFFSET;

        /* Initialize the random number generator */
        idum = 3356;
        rand0 (&idum);


        for (st = state_count; st < ct.num_states; st++)
        {

            /* Generate x, y, z random number sequences */
            for (idx = 0; idx < get_NX_GRID(); idx++)
                xrand[idx] = rand0 (&idum) - 0.5;
            for (idx = 0; idx < get_NY_GRID(); idx++)
                yrand[idx] = rand0 (&idum) - 0.5;
            for (idx = 0; idx < get_NZ_GRID(); idx++)
                zrand[idx] = rand0 (&idum) - 0.5;

            state_p = &states[st];

            idx = 0;
            for (ix = 0; ix < PX0_GRID; ix++)
            {

                for (iy = 0; iy < PY0_GRID; iy++)
                {

                    for (iz = 0; iz < PZ0_GRID; iz++)
                    {

                        state_p->psi[idx] = xrand[xoff + ix] * yrand[yoff + iy] * zrand[zoff + iz];
                        state_p->psi[idx] = state_p->psi[idx] * state_p->psi[idx];
                        idx++;

                    }               /* end for */
                }                   /* end for */
            }                       /* end for */


        }                           /* end for */

        delete [] zrand;
        delete [] yrand;
        delete [] xrand;

    }


}
/******/
