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


#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "FiniteDiff.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "InputKey.h"
#include "Kpoint.h"
#include "RmgException.h"
#include "blas.h"
#include <complex>
#include <omp.h>
#include "../Headers/prototypes.h"

extern "C" void zaxpy(int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);

template Kpoint<double>::Kpoint(double*, double, int, MPI_Comm, BaseGrid *, TradeImages *, Lattice *, std::unordered_map<std::string, InputKey *>& ControlMap);
template Kpoint<std::complex <double> >::Kpoint(double*, double, int, MPI_Comm, BaseGrid *, TradeImages *, Lattice *, std::unordered_map<std::string, InputKey *>& ControlMap);
template void Kpoint<double>::sort_orbitals(void);
template void Kpoint<std::complex <double> >::sort_orbitals(void);
template void Kpoint<double>::set_pool(double *pool);
template void Kpoint<std::complex <double> >::set_pool(std::complex<double> *pool);
template void Kpoint<double>::init_states(void);
template void Kpoint<std::complex <double> >::init_states(void);
template int Kpoint<double>::get_nstates(void);
template int Kpoint<std::complex <double> >::get_nstates(void);
template void Kpoint<double>::random_init(void);
template void Kpoint<std::complex <double> >::random_init(void);
template void Kpoint<double>::orthogonalize(double *tpsi);
template void Kpoint<std::complex <double> >::orthogonalize(std::complex <double> *tpsi);
template void Kpoint<double>::write_occ(void);
template void Kpoint<std::complex <double> >::write_occ(void);

template <class KpointType> Kpoint<KpointType>::Kpoint(double *kkpt, double kkweight, int kindex, MPI_Comm newcomm, BaseGrid *newG, TradeImages *newT, Lattice *newL, std::unordered_map<std::string, InputKey *>& ControlMap) : ControlMap(ControlMap)
{

    this->kpt[0] = kkpt[0];
    this->kpt[1] = kkpt[1];
    this->kpt[2] = kkpt[2];
    this->comm = newcomm;
    this->kidx = kindex;
    this->kweight = kkweight;
    this->nl_weight = NULL;
    this->nl_Bweight = NULL;

    this->G = newG;
    this->T = newT;
    this->L = newL;
    this->pbasis = this->G->get_P0_BASIS(1);

    double v1, v2, v3;
    v1 = twoPI * this->kpt[0] / Rmg_L.get_xside();
    v2 = twoPI * this->kpt[1] / Rmg_L.get_yside();
    v3 = twoPI * this->kpt[2] / Rmg_L.get_zside();

    this->kvec[0] = v1;
    this->kvec[1] = v2;
    this->kvec[2] = v3;
    this->kmag = v1 * v1 + v2 * v2 + v3 * v3;

    this->init_states();

}

#define MAX_NOCC 10
template <class KpointType> void Kpoint<KpointType>::init_states(void)
{

    int ii, is, ns, idx, i, j;
    int count_states[2]={0,0}, nocc[2]={0,0}, num_states_spf[2], nspin = (ct.spin_flag + 1);
    char *tbuf[2];

    struct
    {
        int n;
        double occ;
    } occ[nspin * MAX_NOCC];

    int repeat_occ;


    /* calculate total number of electrons in the system from pseudopotential information */
    ct.ionic_charge = 0.0;

    for (ii = 0; ii < ct.num_ions; ii++)
        ct.ionic_charge += ct.sp[ct.ions[ii].species].zvalence;


    if (ct.spin_flag)
    {
        repeat_occ =( (strcmp(ct.occupation_str_spin_up, "") != 0) && (strcmp(ct.occupation_str_spin_down, "")!= 0) );
        if( (strcmp(ct.occupation_str_spin_up, "") != 0) + (strcmp(ct.occupation_str_spin_down, "")!= 0) == 1 )
                rmg_error_handler (__FILE__, __LINE__, "Fixed occupation for both spin up and down must be specified !!!");
        tbuf[0] = ct.occupation_str_spin_up;
        tbuf[1] = ct.occupation_str_spin_down;
        num_states_spf[0] = 0;
        num_states_spf[1] = 0;
    }
    else
    {
        repeat_occ = (strcmp (ct.occupation_str, "") != 0);
        num_states_spf[0] = 0;
        tbuf[0] = ct.occupation_str;
    }

    /* get repeat count occupation info from control file*/
    if (repeat_occ)
    {
        int n;
        ct.nel = 0;
        ct.nel_up = 0;
        ct.nel_down = 0;

        for (idx = 0; idx < nspin; idx++)
        {
                /* count the fixed occupations of states */
                while ((n = strtol (tbuf[idx], &tbuf[idx], 10)) > 0)
                {
                        count_states[idx] += n;
                        if (nocc[idx] == MAX_NOCC)
                                rmg_error_handler (__FILE__, __LINE__, "Too many blocks in repeat count for state occupations");
                                /* two block example  3 1.0  4 0.0*/
                        occ[nocc[idx] + MAX_NOCC * idx].n = n;
                        occ[nocc[idx] + MAX_NOCC * idx].occ = strtod (tbuf[idx], &tbuf[idx]);
                        ct.nel += n * occ[nocc[idx] + MAX_NOCC * idx].occ;
                        if(idx == 0)ct.nel_up += n * occ[nocc[idx] + MAX_NOCC * idx].occ;
                        if(idx == 1)ct.nel_down += n * occ[nocc[idx] + MAX_NOCC * idx].occ;
                        nocc[idx]++;
                }

                num_states_spf[idx] = count_states[idx];

        }

        if ( (nspin == 2) && (num_states_spf[0] != num_states_spf[1]) )
        {
                rmg_printf("number of states for spin up: %d, number of states for spin down %d\n", num_states_spf[0], num_states_spf[1]);
                rmg_error_handler(__FILE__, __LINE__, "num_of_states_spin_up not equal to num_states_spin_down, you are wasting memory address for extra STATE structures !");
        }

        ct.background_charge = ct.nel - ct.ionic_charge;

    }
    else     /* in case no fixed occupations available, calculate number of states */
    {
        ct.nel = ct.ionic_charge + ct.background_charge;
        for (idx = 0; idx < nspin; idx++) {
            num_states_spf[idx] = (int) ceil(0.5 * ct.nel) + ct.num_unocc_states;
            if(idx == 0) ct.nel_up = 0.5*(ct.nel + 1.0);
            if(idx == 1) ct.nel_down = 0.5*(ct.nel - 1.0);
        }
    }

    /* re-assign the number of states for global variables */
    if(ct.num_states <= 0 ) ct.num_states = num_states_spf[0];

    // When LCAO init is selected we may use more orbitals during the initialization
    // than during the rest of the run so we need to count the number of atomic orbitals
    // here and allocate state structures for the largest possible number of states
    ct.total_atomic_orbitals = CountAtomicOrbitals();
    if (Verify ("start_mode","LCAO Start", ControlMap) || (ct.forceflag == BAND_STRUCTURE)) {
        ct.init_states = ct.total_atomic_orbitals + ct.extra_random_lcao_states;
        if(ct.init_states < ct.num_states) {
            ct.init_states = ct.num_states;
        }
    }
    else {
        ct.init_states = ct.num_states;
    }

    if (Verify ("start_mode", "Modified LCAO Start", ControlMap)) {
        ct.init_states = ct.num_states + ct.extra_random_lcao_states;
    }

    ct.run_states = ct.num_states;

    if(ct.state_block_size > ct.num_states ) ct.state_block_size = ct.num_states;
    // Now figure out some buffer sizes
    ct.max_states = ct.run_states + 3 * ct.state_block_size;
    if (Verify ("kohn_sham_solver", "davidson", ControlMap)) ct.max_states = std::max(ct.max_states, 4*ct.run_states);
    if (Verify ("start_mode","LCAO Start", ControlMap)) ct.max_states = std::max(ct.max_states, ct.init_states);
    if (Verify ("start_mode","Modified LCAO Start", ControlMap)) ct.max_states = std::max(ct.max_states, ct.init_states);

    /* Allocate memory for the state structures */
    this->Kstates = new State<KpointType>[ct.max_states];
    this->nstates = ct.init_states;

    // Set the size of the wavefunction array to allocate. This needs to be 4 times
    // the number of run_states for force calculations right now but will be updated later
    //ct.alloc_states = std::max(4*ct.run_states, ct.max_states);
    ct.alloc_states = std::max(2*ct.init_states + 3 * ct.state_block_size, ct.max_states + 3 * ct.state_block_size);


//    if (verify ("calculation_mode", "Band Structure Only"))
//        nk = 1;
//    else
//        nk = ct.num_kpts;
//    my_malloc (states, ct.num_states * nk, STATE);

    if (nspin == 2)
    {
        ct.num_states_up = num_states_spf[0];
        ct.num_states_down = num_states_spf[1];
    }

    /* set the initial occupations of the states */
    if (repeat_occ)
    {

        if (nspin ==1)
                pct.spinpe = 0;

        for (idx = 0; idx < nspin; idx++)
        {
                ns = 0;
                j = (idx + pct.spinpe) % 2;
                for (i = 0; i < nocc[j]; i++)
                {
                        for (is = ns; is < ns + occ[i + j * MAX_NOCC].n; is++)
                            this->Kstates[is].occupation[idx] = occ[i + j * MAX_NOCC].occ;

                        ns += occ[i + j * MAX_NOCC].n;
                }
        }
    }

    else
    {
        double ne[nspin], oc;

        for (idx = 0; idx < nspin; idx++)
                ne[idx] = ct.nel / ((double) nspin);

        for (idx = 0; idx < nspin; idx++)
        {
                for (is = 0; is < ct.init_states; is++)
                {
                        oc = 0.0;
                        if ( ne[idx] >= (3.0 - nspin) )
                                oc = (3.0 - nspin);
                        else if (ne[idx] >= 0.0)
                                oc = ne[idx];
                        ne[idx] -= oc;
                        this->Kstates[is].occupation[idx] = oc;
                }
        }

    }

}

template <class KpointType> void Kpoint<KpointType>::set_pool(KpointType *pool)
{
    KpointType *tptr;
    int state;

    this->orbital_storage = pool;

    tptr = pool;

    for(state = 0;state < ct.max_states;state++) {
        Kstates[state].set_storage(tptr); 
        Kstates[state].Kptr = this;
        Kstates[state].istate = state;
        tptr += this->pbasis;
    }

}
template <class KpointType> int Kpoint<KpointType>::get_nstates(void)
{
    return this->nstates;
}

template <class KpointType> int Kpoint<KpointType>::get_index(void)
{
    return this->kidx;
}

template <class KpointType> void Kpoint<KpointType>::sort_orbitals(void)
{

    int state;
    double t1;
    State<KpointType> *sp, *sp1;

    for(state = 0;state < this->nstates - 1;state++)
    {
        sp = &this->Kstates[state];
        sp1 = &this->Kstates[state+1];
        if (sp->eig[0] > sp1->eig[0])
        {

            if (((sp->occupation[0] > 0.1) && (sp1->occupation[0] > 0.1))
                || ((sp->occupation[0] < 0.1) && (sp1->occupation[0] < 0.1)))
            {

                t1 = sp->eig[0];
                sp->eig[0] = sp1->eig[0];
                sp1->eig[0] = t1;

                t1 = sp->oldeig[0];
                sp->oldeig[0] = sp1->oldeig[0];
                sp1->oldeig[0] = t1;

                t1 = sp->occupation[0];
                sp->occupation[0] = sp1->occupation[0];
                sp1->occupation[0] = t1;


            }                   /* end if */

        }                       /* end if */

    } 

}


// Generates random initial wavefunctions for this Kpoint
template <class KpointType> void Kpoint<KpointType>::random_init(void)
{

    KpointType ONE_t(1.0);
    int PX0_GRID = Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    double *tmp_psiR = new double[pbasis];
    double *tmp_psiI = new double[pbasis];

    double *xrand = new double[2 * Rmg_G->get_NX_GRID(1)];
    double *yrand = new double[2 * Rmg_G->get_NY_GRID(1)];
    double *zrand = new double[2 * Rmg_G->get_NZ_GRID(1)];

    int factor = 2;
    if(ct.is_gamma) factor = 1;

    // Set state 0 to a constant 
    for(int idx = 0;idx < pbasis;idx++) {

        this->Kstates[0].psi[idx] = ONE_t;

    }
    if(typeid(KpointType) == typeid(std::complex<double>)) {
        for(int idx = 0;idx < pbasis;idx++) {
            double *a = (double *)&this->Kstates[0].psi[idx];
            a[1] = 1.0;
        }
    }

    /* If random start and Fermi occupation, start with
       each state equally occupied  */
    if (ct.occ_flag && (ct.runflag != RESTART))
    {
        /* Set occupation for the first state */
        for (int idx = 0; idx < (ct.spin_flag+1); idx++) {
            this->Kstates[0].occupation[idx] = ct.nel / ((ct.spin_flag+1) * this->nstates);
        }
    }

    long int idum = 3356;
    int xoff = Rmg_G->get_PX_OFFSET(1);
    int yoff = Rmg_G->get_PY_OFFSET(1);
    int zoff = Rmg_G->get_PZ_OFFSET(1);

    /* Initialize the random number generator */
    rand0 (&idum);

    for (int state = 1; state < this->nstates; state++)
    {


        /* Generate x, y, z random number sequences */
        for (int idx = 0; idx < factor*Rmg_G->get_NX_GRID(1); idx++)
            xrand[idx] = rand0 (&idum) - 0.5;
        for (int idx = 0; idx < factor*Rmg_G->get_NY_GRID(1); idx++)
            yrand[idx] = rand0 (&idum) - 0.5;
        for (int idx = 0; idx < factor*Rmg_G->get_NZ_GRID(1); idx++)
            zrand[idx] = rand0 (&idum) - 0.5;


        /* If random start and Fermi occupation, start with
           each state equally occupied  */

        if (ct.occ_flag && (ct.runflag != RESTART))
        {
            for (int idx = 0; idx < (ct.spin_flag+1); idx++) {
                this->Kstates[state].occupation[idx] = ct.nel / ((ct.spin_flag+1) * this->nstates);
            }
        }



        int idx = 0;
        for (int ix = 0; ix < PX0_GRID; ix++)
        {

            for (int iy = 0; iy < PY0_GRID; iy++)
            {

                for (int iz = 0; iz < PZ0_GRID; iz++)
                {

                    
                    tmp_psiR[idx] = xrand[xoff + ix] * 
                                    yrand[yoff + iy] * 
                                    zrand[zoff + iz];
                    tmp_psiR[idx] = tmp_psiR[idx] * tmp_psiR[idx];


                    if(!ct.is_gamma) {

                        tmp_psiI[idx] = xrand[Rmg_G->get_NX_GRID(1) + xoff + ix] * 
                                        yrand[Rmg_G->get_NY_GRID(1) + yoff + iy] * 
                                        zrand[Rmg_G->get_NZ_GRID(1) + zoff + iz];
                        tmp_psiI[idx] = tmp_psiI[idx] * tmp_psiI[idx];

                    }

                    idx++;

                }               /* end for */
            }                   /* end for */
        }                       /* end for */

        // Copy data from tmp_psi into orbital storage
        for(idx = 0;idx < pbasis;idx++) {
            this->Kstates[state].psi[idx] = tmp_psiR[idx];
        }
        if(typeid(KpointType) == typeid(std::complex<double>)) {
            for(idx = 0;idx < pbasis;idx++) {
                double *a = (double *)&this->Kstates[state].psi[idx];
                if(!ct.is_gamma)
                    a[1] = tmp_psiI[idx];
                //a[1] = 0.0;

            }

        }

        // Hit the orbital with the right hand mehrstellen operator which should smooth it a bit
        CPP_app_cir_driver (this->L, this->T, this->Kstates[state].psi, this->Kstates[state].psi, PX0_GRID, PY0_GRID, PZ0_GRID, APP_CI_FOURTH);

    }                           /* end for */


    
    delete [] zrand;
    delete [] yrand;
    delete [] xrand;
    delete [] tmp_psiI;
    delete [] tmp_psiR;

}

template <class KpointType> void Kpoint<KpointType>::orthogonalize(double *tpsi)
{

    RmgTimer RT("Orthogonalization");

    double vel = (double) (this->G->get_NX_GRID(1) * this->G->get_NY_GRID(1) * this->G->get_NZ_GRID(1));
    vel = this->L->get_omega() / vel;


    if(ct.norm_conserving_pp) {

        int st, st1, length, idx, omp_tid;
        double zero = 0.0;
        double one = 1.0;
        double *sarr;
        char *transt = "t";
        char *uplo = "l";

        double *tarr = new double[this->nstates];
        double *global_matrix = new double[this->nstates * this->nstates];

        dsyrk( uplo, transt, &this->nstates, &this->pbasis, &one, this->orbital_storage, &this->pbasis,
                    &zero, global_matrix, &this->nstates);

        /* get the global part */
        length = this->nstates * this->nstates;
        MPI_Allreduce(MPI_IN_PLACE, global_matrix, length, MPI_DOUBLE, MPI_SUM, this->comm);


        /* compute the cholesky factor of the overlap matrix */
        int info;
        dpotrf(uplo, &this->nstates, global_matrix, &this->nstates, &info);
        if (info != 0)
            throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << ". Matrix not positive definite or argument error. Terminating";


        // Get inverse of diagonal elements
        for(st = 0;st < this->nstates;st++) tarr[st] = 1.0 / global_matrix[st + this->nstates * st];


        // This code may look crazy but there is a method to the madness. We copy a slice
        // of the wavefunction array consisting of the values for all orbitals of a given
        // basis point into a temporary array. Then we do the updates on each slice and
        // parallelize over slices with OpenMP. This produces good cache behavior
        // and excellent parformance on the XK6.

        double *darr;
        #pragma omp parallel private(idx,st,st1,omp_tid,sarr)
        {
               omp_tid = omp_get_thread_num();
               if(omp_tid == 0) darr = new double[this->nstates * omp_get_num_threads()];
        #pragma omp barrier

        #pragma omp for schedule(static, 1) nowait
            for(idx = 0;idx < this->pbasis;idx++) {

                sarr = &darr[omp_tid*this->nstates];

                for (st = 0; st < this->nstates; st++) sarr[st] = this->orbital_storage[st*this->pbasis + idx];

                for (st = 0; st < this->nstates; st++) {

                    sarr[st] *= tarr[st];

                    for (st1 = st+1; st1 < this->nstates; st1++) {
                        sarr[st1] -= global_matrix[st1 + this->nstates*st] * sarr[st];
                    }

                }

                for (st = 0; st < this->nstates; st++) this->orbital_storage[st*this->pbasis + idx] = sarr[st];

            }
        }

        delete [] darr;

        double tmp = 1.0 / sqrt(vel);
        idx = this->nstates * this->pbasis;
        for(int idx = 0;idx < this->nstates * this->pbasis;idx++) {
            this->orbital_storage[idx] *= tmp;
        }

        delete [] global_matrix;
        delete [] tarr;

    }
    else {

        int incx=1;
        double *cR = new double[this->nstates];

        for(int ist1 = 0;ist1 < this->nstates;ist1++) {


            // Normalize this orbital
            this->Kstates[ist1].normalize(this->Kstates[ist1].psi, ist1);

            /* This will calculate cR coefficients */
            for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {

                int sidx1 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist1 * pct.num_nonloc_ions * ct.max_nl;
                int sidx2 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist2 * pct.num_nonloc_ions * ct.max_nl;
                double sumpsiR = 0.0;
                double sumbetaR = 0.0;

                int nidx = -1;
                for (int ion = 0; ion < pct.num_owned_ions; ion++)
                {
                    int oion = pct.owned_ions_list[ion];

                    ION *iptr = &ct.ions[oion];
                    SPECIES *sp = &ct.sp[iptr->species];

                    int nh = sp->nh;

                    /* Figure out index of owned ion in nonloc_ions_list array*/
                    do {

                        nidx++;
                        if (nidx >= pct.num_nonloc_ions)
                            rmg_error_handler(__FILE__,__LINE__,"Could not find matching entry in pct.nonloc_ions_list for owned ion");

                    } while (pct.nonloc_ions_list[nidx] != oion);

                    double *qqq = pct.qqq[oion];

                    /* get<beta|psi1> and <beta|psi2> */
                    double *sint1R = &pct.newsintR_local[sidx1 + nidx * ct.max_nl];
                    double *sint2R = &pct.newsintR_local[sidx2 + nidx * ct.max_nl];


                    for (int i = 0; i < nh; i++)
                    {
                        int inh = i * nh;
                        double sri = sint1R[i];

                        for (int j = 0; j < nh; j++)
                        {
                            sumbetaR += qqq[inh + j] * sri * sint2R[j];
                        }                   /*end for j */
                    }                       /*end for i */
                }                           /*end for ion */

                for (int idx = 0; idx < this->pbasis; idx++)
                {
                    sumpsiR = sumpsiR + std::real(this->Kstates[ist2].psi[idx] * std::conj(this->Kstates[ist1].psi[idx]));
                }

                cR[ist2] = vel * sumpsiR + sumbetaR;

            }
            int length = this->nstates - (ist1 + 1);
            /*Sum coefficients over all processors */
            if (length)
            {
                global_sums (&cR[ist1 + 1], &length, pct.grid_comm);
            }
            /*Update wavefunctions */
            for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {
                //update_waves (st1, &st[ist2], ist1, ist2, this->kidx, cR[ist2], cI[ist2]);
  
                KpointType cA(cR[ist2]);
                for(int idx = 0;idx < this->pbasis;idx++) {
                    this->Kstates[ist2].psi[idx] = this->Kstates[ist2].psi[idx] 
                                                  - cA * this->Kstates[ist1].psi[idx]; 
                }
                /* update localized <beta|psi2> */
                for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
                {

                    int lsidx1 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist1 * pct.num_nonloc_ions * ct.max_nl;
                    int lsidx2 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist2 * pct.num_nonloc_ions * ct.max_nl;

                    double *ptr1R = &pct.newsintR_local[lsidx1 + ion * ct.max_nl];
                    double *ptr2R = &pct.newsintR_local[lsidx2 + ion * ct.max_nl];

                    QMD_daxpy (ct.max_nl, -cR[ist2], ptr1R, incx, ptr2R, incx);

                }

            }

        }
        delete [] cR;
    }
    
}

template <class KpointType> void Kpoint<KpointType>::orthogonalize(std::complex<double> *tpsi)
{

   RmgTimer RT("Orthogonalization");

   double vel = (double) (this->G->get_NX_GRID(1) * this->G->get_NY_GRID(1) * this->G->get_NZ_GRID(1));
   vel = this->L->get_omega() / vel;



   if(ct.norm_conserving_pp) {

       int ione = 1;
       std::complex<double> *dr = new std::complex<double>[this->nstates];

       // compute the lower-triangular part of the overlap matrix
       for(int st = 0;st < this->nstates;st++) {

           dr[st] = std::complex<double>(0.0,0.0); 

           // Normalize this orbital
           this->Kstates[st].normalize(this->Kstates[st].psi, st);

           // compute the projection along the remaining vectors
           for(int st1 = st + 1;st1 < this->nstates;st1++) {
               dr[st1] = std::complex<double>(0.0,0.0); 
               for(int idx = 0;idx < this->pbasis;idx++) {
                   dr[st1] = dr[st1] + std::conj(this->orbital_storage[st*this->pbasis + idx]) * this->orbital_storage[st1*this->pbasis + idx];
               }
           }            

           int length = 2 * this->nstates;
           MPI_Allreduce(MPI_IN_PLACE, dr, length, MPI_DOUBLE, MPI_SUM, this->comm);

           std::complex<double> ct1(-vel, 0.0);
           for(int st2=0;st2 < this->nstates;st2++) {
               dr[st2] = ct1 * dr[st2];
           }
           
           for(int st1 = st + 1;st1 < this->nstates;st1++) {
               zaxpy(&this->pbasis, &dr[st1], &this->orbital_storage[st*this->pbasis], &ione, &this->orbital_storage[st1*this->pbasis], &ione);
           }            
           
       }

       delete [] dr;

   }
   else {

      double *cR = new double[this->nstates];
      double *cI = new double[this->nstates];

      for(int ist1 = 0;ist1 < this->nstates;ist1++) {


          // Normalize this orbital
          this->Kstates[ist1].normalize(this->Kstates[ist1].psi, ist1);

          /*This will calculate cR and cI coefficients */
          for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {
              //ortho_get_coeff (st1, &st[ist2], ist1, ist2, this->kidx, &cR[ist2], &cI[ist2]);

              int sidx1 = ist1 * ct.max_nl;
              int sidx2 = ist2 * ct.max_nl;
              double sumpsiR = 0.0;
              double sumpsiI = 0.0;
              double sumbetaR = 0.0;
              double sumbetaI = 0.0;

              int nidx = -1;
              for (int ion = 0; ion < pct.num_owned_ions; ion++)
              {
                  int oion = pct.owned_ions_list[ion];

                  ION *iptr = &ct.ions[oion];
                  SPECIES *sp = &ct.sp[iptr->species];

                  int nh = sp->nh;

                  /* Figure out index of owned ion in nonloc_ions_list array*/
                  do {

                      nidx++;
                      if (nidx >= pct.num_nonloc_ions)
                          rmg_error_handler(__FILE__,__LINE__,"Could not find matching entry in pct.nonloc_ions_list for owned ion");

                  } while (pct.nonloc_ions_list[nidx] != oion);

                  double *qqq = pct.qqq[oion];

                  /* get<beta|psi1> and <beta|psi2> */
                  KpointType *sint1 = &this->newsint_local[sidx1 + nidx * ct.max_nl];
                  KpointType *sint2 = &this->newsint_local[sidx2 + nidx * ct.max_nl];


                  for (int i = 0; i < nh; i++)
                  {
                      int inh = i * nh;
                      double sri = std::real(sint1[i]);
                      double sii = std::imag(sint1[i]);

                      for (int j = 0; j < nh; j++)
                      {
                          sumbetaR += qqq[inh + j] * (sri * std::real(sint2[j]) + sii * std::imag(sint2[j]));
                          sumbetaI += qqq[inh + j] * (sri * std::imag(sint2[j]) - sii * std::real(sint2[j]));
                      }                   /*end for j */
                  }                       /*end for i */
              }                           /*end for ion */

              for (int idx = 0; idx < this->pbasis; idx++)
              {
                  //sumpsiR += (tmp_psi2R[idx] * tmp_psi1R[idx] + tmp_psi2I[idx] * tmp_psi1I[idx]);
                  //sumpsiI += (tmp_psi2I[idx] * tmp_psi1R[idx] - tmp_psi2R[idx] * tmp_psi1I[idx]);
                  sumpsiR = sumpsiR + std::real(this->Kstates[ist2].psi[idx] * std::conj(this->Kstates[ist1].psi[idx]));
                  sumpsiI = sumpsiI + std::imag(this->Kstates[ist2].psi[idx] * std::conj(this->Kstates[ist1].psi[idx]));
              }

              cR[ist2] = vel * sumpsiR + sumbetaR;
              cI[ist2] = vel * sumpsiI + sumbetaI;

          }
          int length = this->nstates - (ist1 + 1);
          /*Sum coefficients over all processors */
          if (length)
          {
              global_sums (&cR[ist1 + 1], &length, pct.grid_comm);
              global_sums (&cI[ist1 + 1], &length, pct.grid_comm);
          }
          /*Update wavefunctions */
          for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {
              //update_waves (st1, &st[ist2], ist1, ist2, this->kidx, cR[ist2], cI[ist2]);

              KpointType cA(cR[ist2], cI[ist2]);
              for(int idx = 0;idx < this->pbasis;idx++) {
                  this->Kstates[ist2].psi[idx] = this->Kstates[ist2].psi[idx] 
                                                - cA * this->Kstates[ist1].psi[idx]; 
              }
              /* update localized <beta|psi2> */
              for (int ion = 0; ion < pct.num_nonloc_ions; ion++)
              {

                  int lsidx1 = ist1 * ct.max_nl;
                  int lsidx2 = ist2 * ct.max_nl;

                  KpointType *ptr1 = &this->newsint_local[lsidx1 + ion * ct.max_nl];
                  KpointType *ptr2 = &this->newsint_local[lsidx2 + ion * ct.max_nl];

                  for(int inh=0;inh < ct.max_nl;inh++) {
                      ptr2[inh] = ptr2[inh] - cA * ptr1[inh];
                  }

              }

          }

      }
      delete [] cI;
      delete [] cR;
   }

}


template <class KpointType> void Kpoint<KpointType>::write_occ(void)
{

    int i, idx, nspin = (ct.spin_flag + 1);

    switch (ct.occ_flag)
    {
        case OCC_NONE:
            break;
        case OCC_FD:
            printf ("\nFERMI-DIRAC OCCUPATION WITH PARAMETERS:");
            printf ("\n  TEMP   = %14.8f", ct.occ_width);
            printf ("\n  MIXING = %14.8f", ct.occ_mix);
            break;
        case OCC_GS:
            printf ("\nGAUSSIAN OCCUPATION WITH PARAMETERS:");
            printf ("\n  TEMP   = %14.8f", ct.occ_width);
            printf ("\n  MIXING = %14.8f", ct.occ_mix);
            break;
        case OCC_EF:
            printf ("\nERROR_FUNCTION OCCUPATION WITH PARAMETERS:");
            printf ("\n  TEMP   = %14.8f", ct.occ_width);
            printf ("\n  MIXING = %14.8f", ct.occ_mix);
            break;
        default:
            rmg_error_handler (__FILE__, __LINE__, "unknown filling procedure");
    }


    for (idx = 0; idx < nspin; idx++)
    {
        if (nspin == 1)
                rmg_printf ("\n\n  STATE OCCUPATIONS :\n");
        else if ((nspin == 2) && (idx == 0))
                rmg_printf ("\n\n  STATE OCCUPATIONS FOR SPIN UP:\n");
        else if ((nspin == 2) && (idx == 1))
                rmg_printf ("\n\n  STATE OCCUPATIONS FOR SPIN DOWN:\n");

        for (i = 0; i < ct.num_states; i++)
                rmg_printf (" %7.2f%s", this->Kstates[i].occupation[idx], ((i % 10 == 9) ? "\n" : ""));

        rmg_printf ("\n\n");

    }

}
