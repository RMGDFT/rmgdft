
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "Kpoint.h"
#include "blas.h"
#include <complex>
#include <omp.h>
#include "../Headers/prototypes.h"

extern "C" void zaxpy(int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);

template Kpoint<double>::Kpoint(double*, double, int, int, MPI_Comm, BaseGrid *, TradeImages *, Lattice *);
template Kpoint<std::complex <double> >::Kpoint(double*, double, int, int, MPI_Comm, BaseGrid *, TradeImages *, Lattice *);
template void Kpoint<double>::sort_orbitals(void);
template void Kpoint<std::complex <double> >::sort_orbitals(void);
template void Kpoint<double>::set_pool(double *pool);
template void Kpoint<std::complex <double> >::set_pool(std::complex<double> *pool);
template int Kpoint<double>::get_nstates(void);
template int Kpoint<std::complex <double> >::get_nstates(void);
template void Kpoint<double>::random_init(void);
template void Kpoint<std::complex <double> >::random_init(void);
template void Kpoint<double>::orthogonalize(double *tpsi);
template void Kpoint<std::complex <double> >::orthogonalize(std::complex <double> *tpsi);


template <class KpointType> Kpoint<KpointType>::Kpoint(double *kkpt, double kkweight, int knstates, int kindex, MPI_Comm newcomm, BaseGrid *newG, TradeImages *newT, Lattice *newL)
{

    this->kpt[0] = kkpt[0];
    this->kpt[1] = kkpt[1];
    this->kpt[2] = kkpt[2];
    this->comm = newcomm;
    this->kidx = kindex;
    this->kweight = kkweight;
    this->nstates = knstates;
    this->Kstates = new State<KpointType>[this->nstates];
    this->G = newG;
    this->T = newT;
    this->L = newL;
    this->pbasis = this->G->get_P0_BASIS(1);
}

template <class KpointType> void Kpoint<KpointType>::set_pool(KpointType *pool)
{
    KpointType *tptr;
    int state;

    this->orbital_storage = pool;

    tptr = pool;

    for(state = 0;state < nstates;state++) {
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
    KpointType *tmp_orbital;

    tmp_orbital = new KpointType[this->pbasis];

    for(state = 0;state < this->nstates - 1;state++)
    {
        sp = this->Kstates + state;
        sp1 = sp++;
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

    delete [] tmp_orbital;
}


// Generates random initial wavefunctions for this Kpoint
template <class KpointType> void Kpoint<KpointType>::random_init(void)
{

    int PX0_GRID = Rmg_G->get_PX0_GRID(1);
    int PY0_GRID = Rmg_G->get_PY0_GRID(1);
    int PZ0_GRID = Rmg_G->get_PZ0_GRID(1);

    int pbasis = PX0_GRID * PY0_GRID * PZ0_GRID;
    double *tmp_psiR = new double[pbasis];
    double *tmp_psiI = new double[pbasis];

    double *xrand = new double[2 * Rmg_G->get_NX_GRID(1)];
    double *yrand = new double[2 * Rmg_G->get_NY_GRID(1)];
    double *zrand = new double[2 * Rmg_G->get_NZ_GRID(1)];

    // Set state 0 to a constant 
    for(int idx = 0;idx < pbasis;idx++) {

        this->Kstates[0].psi[idx] = 1.0;

    }

    /* If random start and Fermi occupation, start with
       each state equally occupied  */
    if (ct.occ_flag && (ct.runflag != 1))
    {
        /* Set occupation for the first state */
        for (int idx = 0; idx < (ct.spin_flag+1); idx++) {
            this->Kstates[0].occupation[idx] = ct.nel / ((ct.spin_flag+1) * ct.num_states);
            this->kstates[0].occupation[idx] = ct.nel / ((ct.spin_flag+1) * ct.num_states);
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
        for (int idx = 0; idx < 2*Rmg_G->get_NX_GRID(1); idx++)
            xrand[idx] = rand0 (&idum) - 0.5;
        for (int idx = 0; idx < 2*Rmg_G->get_NY_GRID(1); idx++)
            yrand[idx] = rand0 (&idum) - 0.5;
        for (int idx = 0; idx < 2*Rmg_G->get_NZ_GRID(1); idx++)
            zrand[idx] = rand0 (&idum) - 0.5;


        /* If random start and Fermi occupation, start with
           each state equally occupied  */

        if (ct.occ_flag && (ct.runflag != 1))
        {
            for (int idx = 0; idx < (ct.spin_flag+1); idx++) {
                this->Kstates[state].occupation[idx] = ct.nel / ((ct.spin_flag+1) * ct.num_states);
                this->kstates[state].occupation[idx] = ct.nel / ((ct.spin_flag+1) * ct.num_states);
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

                    tmp_psiI[idx] = xrand[Rmg_G->get_NX_GRID(1) + xoff + ix] * 
                                    yrand[Rmg_G->get_NY_GRID(1) + yoff + iy] * 
                                    zrand[Rmg_G->get_NZ_GRID(1) + zoff + iz];
                    tmp_psiI[idx] = tmp_psiI[idx] * tmp_psiI[idx];

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
                a[1] = tmp_psiI[idx];
            }

        }

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
       const char *transt = "t";
       const char *uplo = "l";

       double *tarr = new double[this->nstates];
       double *global_matrix = new double[this->nstates * this->nstates];

       ssyrk( uplo, transt, &this->nstates, &this->pbasis, &one, this->orbital_storage, &this->pbasis,
                   &zero, global_matrix, &this->nstates);

       /* get the global part */
       length = this->nstates * this->nstates;
       MPI_Allreduce(MPI_IN_PLACE, global_matrix, length, MPI_DOUBLE, MPI_SUM, this->comm);


       /* compute the cholesky factor of the overlap matrix */
       cholesky(global_matrix, this->nstates);


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


    }
    
}

template <class KpointType> void Kpoint<KpointType>::orthogonalize(std::complex<double> *tpsi)
{

   RmgTimer RT("Orthogonalization");

   if(ct.norm_conserving_pp) {

       int ione = 1;
       std::complex<double> *dr = new std::complex<double>[this->nstates];

       // compute the lower-triangular part of the overlap matrix
       for(int st = 0;st < this->nstates;st++) {

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

           for(int st1 = st + 1;st1 < this->nstates;st1++) {
               zaxpy(&this->pbasis, &dr[st1], &this->orbital_storage[st*this->pbasis], &ione, &this->orbital_storage[st1*this->pbasis], &ione);
           }            
           
       }

       delete [] dr;
   }
}
