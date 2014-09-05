/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/

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
template void Kpoint<double>::mix_betaxpsi(int mix);
template void Kpoint<std::complex <double> >::mix_betaxpsi(int mix);
template void Kpoint<double>::mix_betaxpsi1(int istate);
template void Kpoint<std::complex <double> >::mix_betaxpsi1(int istate);


template <class KpointType> Kpoint<KpointType>::Kpoint(double *kkpt, double kkweight, int knstates, int kindex, MPI_Comm newcomm, BaseGrid *newG, TradeImages *newT, Lattice *newL)
{

    this->kpt[0] = kkpt[0];
    this->kpt[1] = kkpt[1];
    this->kpt[2] = kkpt[2];
    this->comm = newcomm;
    this->kidx = kindex;
    this->kweight = kkweight;
    this->nl_weight = NULL;
    this->nl_Bweight = NULL;
    this->nstates = knstates;
    this->Kstates = new State<KpointType>[this->nstates];
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
    STATE *spl, *spl1;

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

    // Legacy
    for(state = 0;state < this->nstates - 1;state++)
    {
        sp = &this->Kstates[state];
        sp1 = &this->Kstates[state+1];
        spl = &this->kstates[state];
        spl1 = &this->kstates[state+1];
        if (spl->eig[0] > spl1->eig[0])
        {

            if (((spl->occupation[0] > 0.1) && (spl1->occupation[0] > 0.1))
                || ((spl->occupation[0] < 0.1) && (spl1->occupation[0] < 0.1)))
            {
                KpointType tmp;
                for(int idx = 0;idx < this->pbasis;idx++) {
                    tmp = sp->psi[idx];
                    sp->psi[idx] = sp1->psi[idx];
                    sp1->psi[idx] = tmp;
                }

                t1 = spl->eig[0];
                spl->eig[0] = spl1->eig[0];
                spl1->eig[0] = t1;

                t1 = spl->oldeig[0];
                spl->oldeig[0] = spl1->oldeig[0];
                spl1->oldeig[0] = t1;

                t1 = spl->occupation[0];
                spl->occupation[0] = spl1->occupation[0];
                spl1->occupation[0] = t1;


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
    if (ct.occ_flag && (ct.runflag != 1))
    {
        /* Set occupation for the first state */
        for (int idx = 0; idx < (ct.spin_flag+1); idx++) {
            this->Kstates[0].occupation[idx] = ct.nel / ((ct.spin_flag+1) * this->nstates);
            this->kstates[0].occupation[idx] = ct.nel / ((ct.spin_flag+1) * this->nstates);
        }
    }

    long int idum = 3356;
    int xoff = Rmg_G->get_PX_OFFSET(1);
    int yoff = Rmg_G->get_PY_OFFSET(1);
    int zoff = Rmg_G->get_PZ_OFFSET(1);

    /* Initialize the random number generator */
    rand0 (&idum);

    for (int state = 0; state < this->nstates; state++)
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
                this->Kstates[state].occupation[idx] = ct.nel / ((ct.spin_flag+1) * this->nstates);
                this->kstates[state].occupation[idx] = ct.nel / ((ct.spin_flag+1) * this->nstates);
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
                    //tmp_psiR[idx] = tmp_psiR[idx] * tmp_psiR[idx];

                    tmp_psiI[idx] = xrand[Rmg_G->get_NX_GRID(1) + xoff + ix] * 
                                    yrand[Rmg_G->get_NY_GRID(1) + yoff + iy] * 
                                    zrand[Rmg_G->get_NZ_GRID(1) + zoff + iz];
                    //tmp_psiI[idx] = tmp_psiI[idx] * tmp_psiI[idx];

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
                //a[1] = 0.0;

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

#if 0
        // RRRRR Check orthogonalization
        double *varr = new double[this->nstates*this->nstates];
        for(int idx=0;idx<this->nstates*this->nstates;idx++) varr[idx] = 0.0;
        for(int st = 0;st < this->nstates;st++) {
            for(int st1 = 0;st1 < this->nstates;st1++) {
                for(int idx = 0;idx < this->pbasis;idx++) {
                    varr[st*this->nstates + st1] = varr[st*this->nstates + st1] + vel*this->orbital_storage[st*this->pbasis + idx] * this->orbital_storage[st1*this->pbasis + idx];
                }

            }

        }

        int ll = this->nstates*this->nstates;
        global_sums (varr, &ll, pct.grid_comm);

        for(int st = 0;st < this->nstates;st++) {
            for(int st1 = 0;st1 < this->nstates;st1++) {
                rmg_printf("\nOOOOOOO  %d  %d  (%14.6f)\n", st, st1, varr[st*this->nstates + st1]);
            }
        }
#endif

    }
    else {

        int incx=1;
        double *cR = new double[this->nstates];
        STATE *st;

        st = this->kstates;
        for(int ist1 = 0;ist1 < this->nstates;ist1++) {


            // Normalize this orbital
            this->Kstates[ist1].normalize(this->Kstates[ist1].psi, ist1);

            /* This will calculate cR coefficients */
            for (int ist2 = ist1 + 1; ist2 < this->nstates; ist2++) {

                int sidx1 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist1 * ct.max_nl;
                int sidx2 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist2 * ct.max_nl;
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
                    double *sint1R = &pct.newsintR_local[sidx1 + nidx * this->nstates * ct.max_nl];
                    double *sint2R = &pct.newsintR_local[sidx2 + nidx * this->nstates * ct.max_nl];


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

                    int lsidx1 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist1 * ct.max_nl;
                    int lsidx2 = this->kidx * pct.num_nonloc_ions * this->nstates * ct.max_nl + ist2 * ct.max_nl;

                    double *ptr1R = &pct.newsintR_local[lsidx1 + ion * this->nstates * ct.max_nl];
                    double *ptr2R = &pct.newsintR_local[lsidx2 + ion * this->nstates * ct.max_nl];

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

#if 0
       // RRRRR Check orthogonalization
       std::complex<double> *tarr = new std::complex<double>[this->nstates*this->nstates];
       for(int idx=0;idx<this->nstates*this->nstates;idx++) tarr[idx] = std::complex<double>(0.0,0.0);
       for(int st = 0;st < this->nstates;st++) {
           for(int st1 = 0;st1 < this->nstates;st1++) {
               for(int idx = 0;idx < this->pbasis;idx++) {
                   tarr[st*this->nstates + st1] = tarr[st*this->nstates + st1] + vel*std::conj(this->orbital_storage[st*this->pbasis + idx]) * this->orbital_storage[st1*this->pbasis + idx];
               }

           }

       }

       int ll = 2*this->nstates*this->nstates;
       global_sums ((double *)tarr, &ll, pct.grid_comm);

       for(int st = 0;st < this->nstates;st++) {
           for(int st1 = 0;st1 < this->nstates;st1++) {
               rmg_printf("\nOOOOOOO  %d  %d  (%14.6f,%14.6f)\n", st, st1, std::real(tarr[st*this->nstates + st1]), std::imag(tarr[st*this->nstates + st1]));
           }
       }
       delete [] tarr;
#endif
   }
   else {

      int incx=1;
      double *cR = new double[this->nstates];
      double *cI = new double[this->nstates];
      STATE *st;

      st = this->kstates;
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
                  KpointType *sint1 = &this->newsint_local[sidx1 + nidx * this->nstates * ct.max_nl];
                  //double *sint1I = &pct.newsintI_local[sidx1 + nidx * this->nstates * ct.max_nl];
                  KpointType *sint2 = &this->newsint_local[sidx2 + nidx * this->nstates * ct.max_nl];
                  //double *sint2I = &pct.newsintI_local[sidx2 + nidx * this->nstates * ct.max_nl];


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

                  KpointType *ptr1 = &this->newsint_local[lsidx1 + ion * this->nstates * ct.max_nl];
                  KpointType *ptr2 = &this->newsint_local[lsidx2 + ion * this->nstates * ct.max_nl];

                  double *ptr1R = &pct.newsintR_local[lsidx1 + ion * this->nstates * ct.max_nl];
                  //double *ptr1I = &pct.newsintI_local[lsidx1 + ion * this->nstates * ct.max_nl];
                  double *ptr2R = &pct.newsintR_local[lsidx2 + ion * this->nstates * ct.max_nl];
                  //double *ptr2I = &pct.newsintI_local[lsidx2 + ion * this->nstates * ct.max_nl];
#if 0
                  for(int inh=0;inh < ct.max_nl;inh++) {
                      ptr2[inh] = ptr2[inh] + 
                  }
                  QMD_daxpy (ct.max_nl, -cR[ist2], ptr1R, incx, ptr2R, incx);
                  QMD_daxpy (ct.max_nl, cI[ist2], ptr1I, incx, ptr2R, incx);
  
                  QMD_daxpy (ct.max_nl, -cR[ist2], ptr1I, incx, ptr2I, incx);
                  QMD_daxpy (ct.max_nl, -cI[ist2], ptr1R, incx, ptr2I, incx);
#endif

              }

          }

      }
      delete [] cI;
      delete [] cR;
   }

}

template <class KpointType> void Kpoint<KpointType>::mix_betaxpsi(int mix)
{
    if(mix) {

        double scale = 1.0 - ct.prjmix;
        for(int idx = 0;idx < this->sint_size;idx++) {
            this->oldsint_local[idx] = scale * this->oldsint_local[idx];
            this->oldsint_local[idx] = this->oldsint_local[idx] + ct.prjmix * this->newsint_local[idx];
        }

    }
    else {

        for(int idx=0;idx < this->sint_size;idx++) {
            this->oldsint_local[idx] = this->newsint_local[idx];
        }

    }
}

template <class KpointType> void Kpoint<KpointType>::mix_betaxpsi1(int istate)
{

    double scale = pow(1.0 - ct.prjmix, (double)istate);
    if(istate == 0) scale = 1.0 - ct.prjmix;

    for (int ion = 0; ion < pct.num_nonloc_ions; ion++) {

        /* For localized <beta|psi>, there is offset due to the ion*/
        int loffset = ion * this->nstates * ct.max_nl;
        for(int idx=0;idx < ct.max_nl;idx++) 
            this->oldsint_local[loffset + istate * ct.max_nl + idx] *= scale;

        for(int idx=0;idx < ct.max_nl;idx++) {
            this->oldsint_local[loffset + istate * ct.max_nl + idx] += 
                (1.0 - scale) * this->newsint_local[loffset + istate * ct.max_nl + idx];

        }

    }

}
