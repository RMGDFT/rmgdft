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
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>
#include <complex>
#include <omp.h>
#include <boost/filesystem.hpp>
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
#include "ErrorFuncs.h"
#include "GpuAlloc.h"

extern "C" void zaxpy(int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);

template Kpoint<double>::Kpoint(KSTRUCT &kpin, int, MPI_Comm, BaseGrid *, TradeImages *, Lattice *, std::unordered_map<std::string, InputKey *>& ControlMap);
template Kpoint<std::complex <double> >::Kpoint(KSTRUCT &kpin, int, MPI_Comm, BaseGrid *, TradeImages *, Lattice *, std::unordered_map<std::string, InputKey *>& ControlMap);
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
template void Kpoint<double>::get_nlop(int type);
template void Kpoint<std::complex <double> >::get_nlop(int type);
template void Kpoint<double>::reset_beta_arrays(void);
template void Kpoint<std::complex <double> >::reset_beta_arrays(void);
template void Kpoint<double>::reset_orbital_arrays(void);
template void Kpoint<std::complex <double> >::reset_orbital_arrays(void);
template void Kpoint<double>::get_orbitals(double *orbitals);
template void Kpoint<std::complex <double> >::get_orbitals(std::complex<double> *orbitals);
template void Kpoint<double>::get_ion_orbitals(ION *iptr, double *orbitals);
template void Kpoint<std::complex <double> >::get_ion_orbitals(ION *iptr, std::complex<double> *orbitals);
template void Kpoint<double>::get_ldaUop(int);
template void Kpoint<std::complex <double> >::get_ldaUop(int);

template void Kpoint<double>::DeleteNvmeArrays(void);
template void Kpoint<std::complex <double> >::DeleteNvmeArrays(void);
template void Kpoint<double>::ClearPotentialAcceleration(void);
template void Kpoint<std::complex <double> >::ClearPotentialAcceleration(void);


template <class KpointType> Kpoint<KpointType>::Kpoint(KSTRUCT &kpin, int kindex, MPI_Comm newcomm, BaseGrid *newG, TradeImages *newT, Lattice *newL, std::unordered_map<std::string, InputKey *>& ControlMap) : kp(kpin), ControlMap(ControlMap)
{

    this->grid_comm = newcomm;
    this->kidx = kindex;
    this->nl_weight = NULL;
#if HIP_ENABLED || CUDA_ENABLED || SYCL_ENABLED
    this->nl_weight_gpu = NULL;
#endif
    this->orbital_weight = NULL;
    this->BetaProjector = NULL;
    this->OrbitalProjector = NULL;
    this->ldaU = NULL;
    this->newsint_local = NULL;
    this->orbitalsint_local = NULL;
    this->nvme_weight_fd = -1;
    this->nvme_ldaU_fd = -1;
    this->G = newG;
    this->T = newT;
    this->L = newL;
    this->pbasis = this->G->get_P0_BASIS(1);
    this->pbasis_noncoll = ct.noncoll_factor * this->pbasis;

    // This is a boost pool per thread allocator used in MgEigState. It makes it easy
    // to change the type of the underlying memory (cpu, cuda, hip) by modifying
    // rmg_user_allocator as well as the fixed pool performance advantage. Initial
    // size is 24 blocks with max size of 32.
    if(kalloc.size() == 0)
    {
        kalloc.resize(ct.MG_THREADS_PER_NODE);
        int dimx = this->G->get_PX0_GRID(1) * pct.coalesce_factor;
        int dimy = this->G->get_PY0_GRID(1);
        int dimz = this->G->get_PZ0_GRID(1);
        int sbasis = (dimx + 8) * (dimy + 8) * (dimz + 8);
        int sbasis_noncoll = sbasis * ct.noncoll_factor;

        for(int tid=0;tid < ct.MG_THREADS_PER_NODE;tid++)
            kalloc[tid] = new boost::pool<rmg_user_allocator>(sbasis_noncoll*sizeof(KpointType), 24, 32);

    }
    this->init_states();

}

#define MAX_NOCC 10
template <class KpointType> void Kpoint<KpointType>::init_states(void)
{

    int is, ns, idx, i, j;
    int count_states[2]={0,0}, nocc[2]={0,0}, num_states_spf[2], nspin_occ;
    char *tbuf[2];

    nspin_occ = 1;
    if(ct.nspin == 2) nspin_occ = 2;
    struct
    {
        int n;
        double occ;
    } occ[nspin_occ * MAX_NOCC];

    int repeat_occ;


    /* calculate total number of electrons in the system from pseudopotential information */
    ct.ionic_charge = 0.0;

    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
        ct.ionic_charge += Species[Atoms[ion].species].zvalence;


    if (ct.nspin == 2)
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

        for (idx = 0; idx < nspin_occ; idx++)
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

        if ( (nspin_occ == 2) && (num_states_spf[0] != num_states_spf[1]) )
        {
            rmg_printf("number of states for spin up: %d, number of states for spin down %d\n", num_states_spf[0], num_states_spf[1]);
            rmg_error_handler(__FILE__, __LINE__, "num_of_states_spin_up not equal to num_states_spin_down, you are wasting memory address for extra STATE structures !");
        }

        ct.background_charge = ct.nel - ct.ionic_charge;

    }
    else     /* in case no fixed occupations available, calculate number of states */
    {
        ct.nel = ct.ionic_charge + ct.background_charge;
        if(ct.num_unocc_states < 0)
        {
            // num of unocc states set to 10 or 10% which is larger
            ct.num_unocc_states = std::max(10, int(ct.nel * 0.5 * 0.1) );
        }



        for (idx = 0; idx < nspin_occ; idx++) {
            num_states_spf[idx] = (int) ceil(0.5 * ct.nel * ct.noncoll_factor) + ct.num_unocc_states;
            if(idx == 0) ct.nel_up = 0.5*(ct.nel + 1.0);
            if(idx == 1) ct.nel_down = 0.5*(ct.nel - 1.0);
        }

    }

    if(nspin_occ == 1)
    {
        ct.nel_up = 0.5 * ct.nel;
        ct.nel_down = 0.5 * ct.nel;
    }

    /* re-assign the number of states for global variables */
    if(ct.num_states <= 0 ) ct.num_states = num_states_spf[0];

    if(pct.coalesce_factor > 1)
    {
        // Adjust num_states to be an integral multiple of pct.coalesce_factor. Whatever value is used for active threads here
        // must match the value used in MgridSubspace
        int active_threads = ct.MG_THREADS_PER_NODE;
        if(ct.mpi_queue_mode) active_threads--;
        if(active_threads < 1) active_threads = 1;
        int states_div = (ct.num_states / (active_threads*pct.coalesce_factor)) * active_threads*pct.coalesce_factor;
        int states_rem = ct.num_states % (active_threads*pct.coalesce_factor);
        if(states_rem)
        {
            while(states_div < ct.num_states) states_div += pct.coalesce_factor;
            if((pct.gridpe == 0) && (states_div != ct.num_states))
                std::cout << "Notice: Numstates adjusted from " << ct.num_states << " to " << states_div <<  " to satisfy coalesce conditions." << std::endl;
            ct.num_states = states_div;
        }
    }

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
    if (Verify ("kohn_sham_solver", "davidson", ControlMap)) 
    {
        // Check if user specified autosetting of davidson multiplier
        if(ct.davidx <= 1)
        {
            ct.davidx = 4;
            if(ct.num_states > 600) ct.davidx = 3;
            if(ct.num_states > 900) ct.davidx = 2;
        }
        ct.max_states = std::max(ct.max_states, ct.davidx*ct.run_states);
    }
    if (Verify ("start_mode","LCAO Start", ControlMap)) ct.max_states = std::max(ct.max_states, 2*ct.init_states);
    if (Verify ("start_mode","Modified LCAO Start", ControlMap)) ct.max_states = std::max(ct.max_states, ct.init_states);
    if(ct.forceflag == BAND_STRUCTURE) ct.max_states = std::max(ct.max_states, 3*ct.num_states);
    ct.max_states = std::max(ct.max_states, 2*ct.num_states);


    /* Allocate memory for the state structures */
    this->Kstates = new State<KpointType>[ct.max_states];
    this->nstates = ct.init_states;

    // Set the size of the wavefunction array to allocate. This needs to be 4 times
    // the number of run_states for force calculations right now but will be updated later
    //ct.alloc_states = std::max(4*ct.run_states, ct.max_states);
    ct.alloc_states = ct.run_states + 3 * ct.state_block_size;
    ct.alloc_states = std::max(2*ct.init_states, ct.alloc_states);
    if(ct.alloc_states < ct.max_states) ct.alloc_states = ct.max_states;

    //    if (verify ("calculation_mode", "Band Structure Only"))
    //        nk = 1;
    //    else
    //        nk = ct.num_kpts;
    //    my_malloc (states, ct.num_states * nk, STATE);

    if (nspin_occ == 2)
    {
        ct.num_states_up = num_states_spf[0];
        ct.num_states_down = num_states_spf[1];
    }

    /* set the initial occupations of the states */
    if (repeat_occ)
    {

        for (idx = 0; idx < nspin_occ; idx++)
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
        double ne[2], oc;

        ne[0] = ct.nel;
        if(nspin_occ == 2)
        {
            ne[0] = ct.nel/2.0;
            ne[1] = ct.nel/2.0;


            for (idx = 0; idx < nspin_occ; idx++)
            {
                for (is = 0; is < ct.init_states; is++)
                {
                    if(ne[idx] >= 1.0) oc = 1.0;
                    else if(ne[idx] > 0.0) oc = ne[idx];
                    else oc = 0;

                    ne[idx] -= oc;
                    this->Kstates[is].occupation[idx] = oc;
                }
            }
        }

        else
        {
            int num_fullocc = (int)(ct.nel /2.0 * ct.noncoll_factor);
            for (is = 0; is < ct.init_states; is++)
            {
                this->Kstates[is].occupation[0] = 0.0;
            }
            for (is = 0; is < num_fullocc; is++)
            {
                this->Kstates[is].occupation[0] = 2.0/ct.noncoll_factor;
            }
            this->Kstates[num_fullocc+1].occupation[0] = ct.nel - num_fullocc * 2.0 / ct.noncoll_factor;

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
        tptr += this->pbasis * ct.noncoll_factor;
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


// Generates random initial wavefunctions for this Kpoint
template <class KpointType> void Kpoint<KpointType>::random_init(void)
{

    KpointType ONE_t(1.0);
    int PX0_GRID = Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    double *tmp_psiR = new double[ct.noncoll_factor *pbasis];
    double *tmp_psiI = new double[ct.noncoll_factor *pbasis];

    double *xrand = new double[ct.noncoll_factor *2 * Rmg_G->get_NX_GRID(1)];
    double *yrand = new double[ct.noncoll_factor *2 * Rmg_G->get_NY_GRID(1)];
    double *zrand = new double[ct.noncoll_factor *2 * Rmg_G->get_NZ_GRID(1)];

    int factor = 2;
    if(ct.is_gamma) factor = 1;

    factor = factor * ct.noncoll_factor;
    // Set state 0 to a constant 
    for(int idx = 0;idx < ct.noncoll_factor *pbasis;idx++) {

        this->Kstates[0].psi[idx] = ONE_t;

    }
    if(typeid(KpointType) == typeid(std::complex<double>)) {
        for(int idx = 0;idx <ct.noncoll_factor * pbasis;idx++) {
            double *a = (double *)&this->Kstates[0].psi[idx];
            a[1] = 1.0;
        }
    }

    /* If random start and Fermi occupation, start with
       each state equally occupied  */
    if (ct.occ_flag && (ct.runflag != RESTART))
    {
        /* Set occupation for the first state */
        this->Kstates[0].occupation[0] = ct.nel / (this->nstates);
        this->Kstates[0].occupation[1] = ct.nel / (this->nstates);
        if(ct.nspin == 2) 
        {
            this->Kstates[0].occupation[0] /=2.0;
            this->Kstates[0].occupation[1] /=2.0;
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
            this->Kstates[state].occupation[0] = ct.nel / (this->nstates);
            this->Kstates[state].occupation[1] = ct.nel / (this->nstates);
            if(ct.nspin == 2) 
            {
                this->Kstates[state].occupation[0] /=2.0;
                this->Kstates[state].occupation[1] /=2.0;
            }
        }



        for (int ix = 0; ix < PX0_GRID; ix++)
        {

            for (int iy = 0; iy < PY0_GRID; iy++)
            {

                for (int iz = 0; iz < PZ0_GRID; iz++)
                {


                    for(int inoncoll = 0; inoncoll < ct.noncoll_factor; inoncoll++)
                    {
                        int idxoff = inoncoll * pbasis;
                        int idx = idxoff + ix * PY0_GRID * PZ0_GRID + iy * PZ0_GRID + iz;

                        int ixoff = inoncoll * 2 * Rmg_G->get_NX_GRID(1);
                        int iyoff = inoncoll * 2 * Rmg_G->get_NY_GRID(1);
                        int izoff = inoncoll * 2 * Rmg_G->get_NZ_GRID(1);

                        tmp_psiR[idx] = xrand[ixoff +  xoff + ix] * 
                            yrand[iyoff + yoff + iy] * 
                            zrand[izoff + zoff + iz];
                        tmp_psiR[idx] = tmp_psiR[idx] * tmp_psiR[idx];


                        if(!ct.is_gamma) {
                            ixoff = (inoncoll * 2 + 1) * Rmg_G->get_NX_GRID(1);
                            iyoff = (inoncoll * 2 + 1) * Rmg_G->get_NY_GRID(1);
                            izoff = (inoncoll * 2 + 1) * Rmg_G->get_NZ_GRID(1);

                            tmp_psiI[idx] = xrand[ixoff + xoff + ix] * 
                                yrand[iyoff + yoff + iy] * 
                                zrand[izoff + zoff + iz];
                            tmp_psiI[idx] = tmp_psiI[idx] * tmp_psiI[idx];

                        }



                    }               /* end for */
                }                   /* end for */
            }                       /* end for */

        }
        // Copy data from tmp_psi into orbital storage
        for(int idx = 0;idx < ct.noncoll_factor * pbasis;idx++) {
            this->Kstates[state].psi[idx] = tmp_psiR[idx];
        }
        if(typeid(KpointType) == typeid(std::complex<double>)) {
            for(int idx = 0;idx < ct.noncoll_factor * pbasis;idx++) {
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
    int num_nonloc_ions = this->BetaProjector->get_num_nonloc_ions();
    int num_owned_ions = this->BetaProjector->get_num_owned_ions();
    int *owned_ions_list = this->BetaProjector->get_owned_ions_list();
    int *nonloc_ions_list = this->BetaProjector->get_nonloc_ions_list();


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
        MPI_Allreduce(MPI_IN_PLACE, global_matrix, length, MPI_DOUBLE, MPI_SUM, grid_comm);


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

                int sidx1 = ist1 * num_nonloc_ions * ct.max_nl;
                int sidx2 = ist2 * num_nonloc_ions * ct.max_nl;
                double sumpsiR = 0.0;
                double sumbetaR = 0.0;

                int nidx = -1;
                for (int ion = 0; ion < num_owned_ions; ion++)
                {
                    int oion = owned_ions_list[ion];

                    ION *iptr = &Atoms[oion];
                    SPECIES *sp = &Species[iptr->species];

                    int nh = sp->nh;

                    /* Figure out index of owned ion in nonloc_ions_list array*/
                    do {

                        nidx++;
                        if (nidx >= num_nonloc_ions)
                            rmg_error_handler(__FILE__,__LINE__,"Could not find matching entry in nonloc_ions_list for owned ion");

                    } while (nonloc_ions_list[nidx] != oion);

                    double *qqq = Atoms[oion].qqq;

                    /* get<beta|psi1> and <beta|psi2> */
                    double *sint1R = &this->newsint_local[sidx1 + nidx * ct.max_nl];
                    double *sint2R = &this->newsint_local[sidx2 + nidx * ct.max_nl];


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
                global_sums (&cR[ist1 + 1], &length, grid_comm);
            }
            /*Update wavefunctions */
            for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {

                KpointType cA(cR[ist2]);
                for(int idx = 0;idx < this->pbasis;idx++) {
                    this->Kstates[ist2].psi[idx] = this->Kstates[ist2].psi[idx] 
                        - cA * this->Kstates[ist1].psi[idx]; 
                }
                /* update localized <beta|psi2> */
                for (int ion = 0; ion < num_nonloc_ions; ion++)
                {

                    int lsidx1 = ist1 * num_nonloc_ions * ct.max_nl;
                    int lsidx2 = ist2 * num_nonloc_ions * ct.max_nl;

                    double *ptr1R = &this->newsint_local[lsidx1 + ion * ct.max_nl];
                    double *ptr2R = &this->newsint_local[lsidx2 + ion * ct.max_nl];

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
    int num_nonloc_ions = this->BetaProjector->get_num_nonloc_ions();
    int num_owned_ions = this->BetaProjector->get_num_owned_ions();
    int *owned_ions_list = this->BetaProjector->get_owned_ions_list();
    int *nonloc_ions_list = this->BetaProjector->get_nonloc_ions_list();



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
            MPI_Allreduce(MPI_IN_PLACE, dr, length, MPI_DOUBLE, MPI_SUM, grid_comm);

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

                int sidx1 = ist1 * ct.max_nl;
                int sidx2 = ist2 * ct.max_nl;
                double sumpsiR = 0.0;
                double sumpsiI = 0.0;
                double sumbetaR = 0.0;
                double sumbetaI = 0.0;

                int nidx = -1;
                for (int ion = 0; ion < num_owned_ions; ion++)
                {
                    int oion = owned_ions_list[ion];

                    ION *iptr = &Atoms[oion];
                    SPECIES *sp = &Species[iptr->species];

                    int nh = sp->nh;

                    /* Figure out index of owned ion in nonloc_ions_list array*/
                    do {

                        nidx++;
                        if (nidx >= num_nonloc_ions)
                            rmg_error_handler(__FILE__,__LINE__,"Could not find matching entry in pct.nonloc_ions_list for owned ion");

                    } while (nonloc_ions_list[nidx] != oion);

                    double *qqq = Atoms[oion].qqq;

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
                global_sums (&cR[ist1 + 1], &length, grid_comm);
                global_sums (&cI[ist1 + 1], &length, grid_comm);
            }
            /*Update wavefunctions */
            for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {

                KpointType cA(cR[ist2], cI[ist2]);
                for(int idx = 0;idx < this->pbasis;idx++) {
                    this->Kstates[ist2].psi[idx] = this->Kstates[ist2].psi[idx] 
                        - cA * this->Kstates[ist1].psi[idx]; 
                }
                /* update localized <beta|psi2> */
                for (int ion = 0; ion < num_nonloc_ions; ion++)
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


// Used to generate a set of delocalized atomic orbitals for a specific ion.
template <class KpointType> void Kpoint<KpointType>::get_ion_orbitals(ION *iptr, KpointType *orbitals)
{

    std::complex<double> I_t(0.0, 1.0);

    int pbasis = this->pbasis;
    int nlxdim = get_NX_GRID();
    int nlydim = get_NY_GRID();
    int nlzdim = get_NZ_GRID();

    /* Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);

    if ((beptr == NULL) || (gbptr == NULL))
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];

    /* Get species type */
    SPECIES *sp = &Species[iptr->species];

    double vect[3], nlcrds[3];

    /* Find nlcdrs, vector that gives shift of ion from center of its ionic box */
    /* for delocalized case it's just half the cell dimensions */
    vect[0] = iptr->xtal[0] - 0.5;
    vect[1] = iptr->xtal[1] - 0.5;
    vect[2] = iptr->xtal[2] - 0.5;

    /*The vector we are looking for should be */
    to_cartesian (vect, nlcrds);

    /*Calculate the phase factor */
    FindPhaseKpoint (this->kp.kvec, nlxdim, nlydim, nlzdim, nlcrds, fftw_phase, false);

    /*Temporary pointer to the already calculated forward transform */
    fptr = (std::complex<double> *)&sp->forward_orbital[this->kidx * sp->num_orbitals * pbasis];

    KpointType *orbit = orbitals;
    /* Loop over atomic orbitals */
    for (int ip = 0; ip < sp->num_orbitals; ip++)
    {
        /*Apply the phase factor */
        for (int idx = 0; idx < pbasis; idx++) gbptr[idx] = fptr[idx] * std::conj(fftw_phase[idx]);

        /*Do the backwards transform */
        coarse_pwaves->FftInverse(gbptr, beptr);

        std::complex<double> *orbit_C = (std::complex<double> *)orbit;
        double *orbit_R = (double *)orbit;

        if(ct.is_gamma)
        {
            for (int idx = 0; idx < pbasis; idx++) orbit_R[idx] = std::real(beptr[idx]);
        }
        else
        {
            for (int idx = 0; idx < pbasis; idx++) orbit_C[idx] = beptr[idx];
        }

        /*Advance the temp pointers */
        fptr += pbasis;
        orbit += pbasis;

    } 


    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);


} // end get_ionic_orbitals


// Used to generate a set of atomic orbitals that span the full grid. To use these orbitals
// as projectors you need to create a projector object with a stride of ct.max_orbitals
// and use the provided orbitals array. This routine should not be used to create a set
// of atomic orbitals for LCAO starts since the stride padding is neither needed or desired there.
template <class KpointType> void Kpoint<KpointType>::get_orbitals(KpointType *orbitals)
{

    std::complex<double> I_t(0.0, 1.0);

    int pbasis = this->pbasis;
    int nlxdim = get_NX_GRID();
    int nlydim = get_NY_GRID();
    int nlzdim = get_NZ_GRID();

    /* Pointer to the result of forward transform on the coarse grid */
    std::complex<double> *fptr;
    std::complex<double> *beptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);
    std::complex<double> *gbptr = (std::complex<double> *)fftw_malloc(sizeof(std::complex<double>) * pbasis);

    if ((beptr == NULL) || (gbptr == NULL))
        rmg_error_handler (__FILE__, __LINE__, "can't allocate memory\n");

    std::complex<double> *fftw_phase = new std::complex<double>[pbasis];

    /* Loop over ions */
    for (size_t ion = 0, i_end = Atoms.size(); ion < i_end; ++ion)
    {
        KpointType *orbit = &orbitals[ion * ct.max_orbitals * pbasis];

        /* Generate ion pointer */
        ION *iptr = &Atoms[ion];

        /* Get species type */
        SPECIES *sp = &Species[iptr->species];


        /*Calculate the phase factor */
        FindPhaseKpoint (this->kp.kvec, nlxdim, nlydim, nlzdim, this->BetaProjector->nlcrds[ion].data(), fftw_phase, false);

        /*Temporary pointer to the already calculated forward transform */
        fptr = (std::complex<double> *)&sp->forward_orbital[this->kidx * sp->num_orbitals * pbasis];

        /* Loop over atomic orbitals */
        for (int ip = 0; ip < sp->num_orbitals; ip++)
        {

            /*Apply the phase factor */
            for (int idx = 0; idx < pbasis; idx++) gbptr[idx] = fptr[idx] * std::conj(fftw_phase[idx]);

            /*Do the backwards transform */
            coarse_pwaves->FftInverse(gbptr, beptr);

            std::complex<double> *orbit_C = (std::complex<double> *)orbit;
            double *orbit_R = (double *)orbit;

            if(ct.is_gamma)
            {
                for (int idx = 0; idx < pbasis; idx++) orbit_R[idx] = std::real(beptr[idx]);
            }
            else
            {
                for (int idx = 0; idx < pbasis; idx++) orbit_C[idx] = beptr[idx];
            }

            /*Advance the temp pointers */
            fptr += pbasis;
            orbit += pbasis;

        } 

    }                           /* end for */

    delete [] fftw_phase;
    fftw_free (gbptr);
    fftw_free (beptr);


} // end get_orbitals


// Sets up weight and sint arrays for the beta functions
template <class KpointType> void Kpoint<KpointType>::get_nlop(int projector_type)
{

    if(this->BetaProjector) delete this->BetaProjector;
    reset_beta_arrays ();

    this->BetaProjector = new Projector<KpointType>(projector_type, ct.max_nl, BETA_PROJECTOR);
    int num_nonloc_ions = this->BetaProjector->get_num_nonloc_ions();

    if(ct.nvme_weights)
    {
        if(nvme_weight_fd != -1) close(nvme_weight_fd);

        nvme_weight_path = ct.nvme_weights_path + std::string("rmg_weight") + std::to_string(pct.spinpe) + "_" +
            std::to_string(pct.kstart + this->kidx) + "_" + std::to_string(pct.gridpe);
        nvme_weight_fd = FileOpenAndCreate(nvme_weight_path, O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);
    }

    this->nl_weight_size = (size_t)this->BetaProjector->get_num_tot_proj() * (size_t)this->pbasis + 128;
    ct.beta_alloc[0] = this->nl_weight_size * sizeof(KpointType);
    MPI_Allreduce(&ct.beta_alloc[0], &ct.beta_alloc[1], 1, MPI_LONG, MPI_MIN, grid_comm);
    MPI_Allreduce(&ct.beta_alloc[0], &ct.beta_alloc[2], 1, MPI_LONG, MPI_MAX, grid_comm);
    MPI_Allreduce(MPI_IN_PLACE, &ct.beta_alloc[0], 1, MPI_LONG, MPI_SUM, grid_comm);
    ct.beta_alloc[0] *= (size_t)(ct.num_kpts_pe * pct.pe_kpoint * ct.nspin);
    ct.beta_alloc[1] *= (size_t)(ct.num_kpts_pe * pct.pe_kpoint * ct.nspin);
    ct.beta_alloc[2] *= (size_t)(ct.num_kpts_pe * pct.pe_kpoint * ct.nspin);

    // when we calculate stress, we need beta and x*beta, y*beta, z*beta, so that allocate memory for 4 timrs of nl_weight_size
    int stress_factor = 1;
    if(ct.stress) stress_factor = 4;

#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    gpuError_t custat;
    // Managed memory is faster when gpu memory is not constrained but 
    // pinned memory works better when it is constrained.
    if(ct.pin_nonlocal_weights)
    {
        custat = gpuMallocHost((void **)&this->nl_weight, stress_factor * this->nl_weight_size * sizeof(KpointType));
        RmgGpuError(__FILE__, __LINE__, custat, "Error: gpuMallocHost failed.\n");
    }
    else
    {
        this->nl_weight = (KpointType *)RmgMallocHost(stress_factor * this->nl_weight_size * sizeof(KpointType));
#if CUDA_ENABLED
        int device = -1;
        gpuGetDevice(&device);
        cudaMemAdvise ( this->nl_weight, stress_factor * this->nl_weight_size * sizeof(KpointType), cudaMemAdviseSetReadMostly, device);
        custat = cudaMalloc((void **)&this->nl_weight_gpu, stress_factor * this->nl_weight_size * sizeof(KpointType));
        RmgGpuError(__FILE__, __LINE__, custat, "Error: gpuMalloc failed.\n");
#elif HIP_ENABLED
        int device = -1;
        gpuGetDevice(&device);
        hipMemAdvise ( this->nl_weight, stress_factor * this->nl_weight_size * sizeof(KpointType), hipMemAdviseSetReadMostly, device);
        custat = hipMalloc((void **)&this->nl_weight_gpu, stress_factor * this->nl_weight_size * sizeof(KpointType));
        RmgGpuError(__FILE__, __LINE__, custat, "Error: gpuMalloc failed.\n");
#elif SYCL_ENABLED
        custat = gpuMalloc((void **)&this->nl_weight_gpu, stress_factor * this->nl_weight_size * sizeof(KpointType));

#endif
    }
    for(size_t idx = 0;idx < stress_factor * this->nl_weight_size;idx++) this->nl_weight[idx] = 0.0;
#else
    if(ct.nvme_weights)
    {
        this->nl_weight = (KpointType *)CreateMmapArray(nvme_weight_fd, this->nl_weight_size*sizeof(KpointType));
        if(!this->nl_weight) rmg_error_handler(__FILE__,__LINE__,"Error: CreateMmapArray failed for weights. \n");
        madvise(this->nl_weight, stress_factor * this->nl_weight_size*sizeof(KpointType), MADV_NORMAL);
    }
    else
    {
        this->nl_weight = new KpointType[stress_factor * this->nl_weight_size]();
    }
#endif

    int factor = 2;
    if(ct.is_gamma) factor = 1; 
    size_t sint_alloc = (size_t)(factor * num_nonloc_ions * this->BetaProjector->get_pstride() * ct.noncoll_factor);
    sint_alloc *= (size_t)ct.max_states;
    sint_alloc += 16;    // In case of lots of vacuum make sure something is allocated otherwise allocation routine may fail
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    this->newsint_local = (KpointType *)RmgMallocHost(sint_alloc * sizeof(KpointType));
#else
    this->newsint_local = new KpointType[sint_alloc]();
#endif

    MPI_Barrier(grid_comm);

} 


template <class KpointType> void Kpoint<KpointType>::reset_beta_arrays(void)
{

    if (this->nl_weight != NULL) {
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
        if(ct.pin_nonlocal_weights)
        {
            gpuFreeHost(this->nl_weight);
        }
        else
        {
            gpuFreeHost(this->nl_weight);
            gpuFree(this->nl_weight_gpu);
        }
        if (this->newsint_local)
            RmgFreeHost(this->newsint_local);
#else
        int stress_factor = 1;
        if(ct.stress) stress_factor = 4;
        if(ct.nvme_weights)
        {
            munmap(this->nl_weight, stress_factor * this->nl_weight_size*sizeof(double));
        }
        else
        {
            delete [] this->nl_weight;
        }
        if (this->newsint_local)
            delete [] this->newsint_local;
#endif
        this->nl_weight = NULL;
        this->newsint_local = NULL;
    }

}


template <class KpointType> void Kpoint<KpointType>::reset_orbital_arrays(void)
{

    if (this->orbital_weight != NULL) {
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
        RmgFreeHost(this->orbital_weight);
        RmgFreeHost(this->orbitalsint_local);
#else
        if(ct.nvme_weights)
        {
            munmap(this->orbital_weight, this->orbital_weight_size*sizeof(double));
        }
        else
        {
            delete [] this->orbital_weight;
            delete [] this->orbitalsint_local;
        }
#endif
        this->orbital_weight = NULL;
        this->orbitalsint_local = NULL;
    }
}


// Sets up weight and sint arrays for the beta functions
template <class KpointType> void Kpoint<KpointType>::get_ldaUop(int projector_type)
{

    if(this->OrbitalProjector) delete this->OrbitalProjector;
    reset_orbital_arrays ();

    if(this->ldaU) delete this->ldaU;

    this->ldaU = new LdaU<KpointType>(*this);

    //  Can make this more efficient at some point by restricting to ct.num_ldaU_ions but that does not yet work
    this->OrbitalProjector = new Projector<KpointType>(projector_type, ct.max_ldaU_orbitals, ORBITAL_PROJECTOR);
    int num_nonloc_ions = this->OrbitalProjector->get_num_nonloc_ions();

    if(ct.nvme_weights)
    {
        if(nvme_ldaU_fd != -1) close(nvme_ldaU_fd);

        nvme_ldaU_path = ct.nvme_weights_path + std::string("rmg_orbital_weight") + std::to_string(pct.spinpe) + "_" +
            std::to_string(pct.kstart + this->kidx) + "_" + std::to_string(pct.gridpe);
        nvme_ldaU_fd = FileOpenAndCreate(nvme_ldaU_path, O_RDWR|O_CREAT|O_TRUNC, (mode_t)0600);

    }

    this->orbital_weight_size = (size_t)this->OrbitalProjector->get_num_tot_proj() * (size_t)this->pbasis * (size_t)ct.noncoll_factor;

#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    this->orbital_weight = (KpointType *)RmgMallocHost(this->orbital_weight_size * sizeof(KpointType));
#if CUDA_ENABLED
    int device = -1;
    gpuGetDevice(&device);
    cudaMemAdvise ( this->orbital_weight, this->orbital_weight_size * sizeof(KpointType), cudaMemAdviseSetReadMostly, device);
#elif HIP_ENABLED
    int device = -1;
    gpuGetDevice(&device);
    hipMemAdvise ( this->orbital_weight, this->orbital_weight_size * sizeof(KpointType), hipMemAdviseSetReadMostly, device);
#endif

    for(size_t idx = 0;idx < this->orbital_weight_size;idx++) this->orbital_weight[idx] = 0.0;

#else
    if(ct.nvme_weights)
    {
        this->orbital_weight = (KpointType *)CreateMmapArray(nvme_weight_fd, this->orbital_weight_size*sizeof(KpointType));
        if(!this->orbital_weight) rmg_error_handler(__FILE__,__LINE__,"Error: CreateMmapArray failed for weights. \n");
        madvise(this->orbital_weight, this->orbital_weight_size*sizeof(KpointType), MADV_NORMAL);

    }
    else
    {
        this->orbital_weight = new KpointType[this->orbital_weight_size]();
    }
#endif


    int factor = 2;
    if(ct.is_gamma) factor = 1; 
    size_t sint_alloc = (size_t)(factor * num_nonloc_ions * this->OrbitalProjector->get_pstride());
    sint_alloc *= (size_t)ct.max_states;
    sint_alloc += 16;    // In case of lots of vacuum make sure something is allocated otherwise allocation routine may fail
#if CUDA_ENABLED || HIP_ENABLED || SYCL_ENABLED
    this->orbitalsint_local = (KpointType *)RmgMallocHost(sint_alloc * sizeof(KpointType));
#else
    this->orbitalsint_local = new KpointType[sint_alloc]();
#endif

    MPI_Barrier(grid_comm);

} 

// Cleans up nvme arrays if they have been used
template <class KpointType> void Kpoint<KpointType>::DeleteNvmeArrays(void)
{
    reset_beta_arrays();
    reset_orbital_arrays();
    if(nvme_weight_fd > 0)
    {
        std::string weight_path = ct.nvme_weights_path + std::string("rmg_weight") + std::to_string(pct.spinpe) + "_" +
            std::to_string(pct.kstart + this->kidx) + "_" + std::to_string(pct.gridpe);

        if(boost::filesystem::exists(weight_path.c_str()))
        {
            unlink(weight_path.c_str());
        }
        nvme_weight_fd = -1;
    }

    if(nvme_ldaU_fd > 0)
    {
        std::string orbital_path = ct.nvme_weights_path + std::string("rmg_orbital_weight") + std::to_string(pct.spinpe) + "_" +
            std::to_string(pct.kstart + this->kidx) + "_" + std::to_string(pct.gridpe);

        if(boost::filesystem::exists(orbital_path.c_str()))
        {
            unlink(orbital_path.c_str());
        }
        nvme_ldaU_fd = -1;
    }

}

template <class KpointType> void Kpoint<KpointType>::ClearPotentialAcceleration(void)
{
    for(size_t idx=0;idx < this->dvh_size;idx++)
    {
        this->dvh[idx] = 0.0;
    }
}

// Instantiate the allocator
template<typename T>
std::vector<boost::pool<rmg_user_allocator> *> Kpoint<T>::kalloc;


