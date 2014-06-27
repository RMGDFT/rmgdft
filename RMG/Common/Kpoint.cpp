
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
#include "../Headers/prototypes.h"


template class Kpoint<float>;
template Kpoint<double>::Kpoint(double*, double, int, int, int);
template Kpoint<std::complex <double> >::Kpoint(double*, double, int, int, int);
template void Kpoint<double>::sort_orbitals(void);
template void Kpoint<std::complex <double> >::sort_orbitals(void);
template void Kpoint<double>::set_pool(double *pool);
template void Kpoint<std::complex <double> >::set_pool(std::complex<double> *pool);
template int Kpoint<double>::get_nstates(void);
template int Kpoint<std::complex <double> >::get_nstates(void);


template <class KpointType> Kpoint<KpointType>::Kpoint(double *kkpt, double kkweight, int knstates, int kpbasis, int kindex)
{

    this->kpt[0] = kkpt[0];
    this->kpt[1] = kkpt[1];
    this->kpt[2] = kkpt[2];
    this->kidx = kindex;
    this->kweight = kkweight;
    this->nstates = knstates;
    this->pbasis = kpbasis;
    this->Kstates = new State<KpointType>[this->nstates];
}

template <class KpointType> void Kpoint<KpointType>::set_pool(KpointType *pool)
{
    KpointType *tptr;
    int state;

    this->orbital_storage = pool;

    tptr = pool;

    for(state = 0;state < nstates;state++) {
        Kstates[state].set_storage(tptr); 
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

