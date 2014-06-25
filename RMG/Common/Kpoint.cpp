
#include <Kpoint.h>
#include <complex>

using namespace std;

template class Kpoint<float>;
template Kpoint<double>::Kpoint(double*, double, int, int, int);
template Kpoint<complex <double> >::Kpoint(double*, double, int, int, int);

template <class KpointType> Kpoint<KpointType>::Kpoint(double *kkpt, double kkweight, int knstates, int kpbasis, int kindex)
{

    this->kpt[0] = kkpt[0];
    this->kpt[1] = kkpt[1];
    this->kpt[2] = kkpt[2];
    this->kidx = kindex;
    this->kweight = kkweight;
    this->nstates = knstates;
    this->pbasis = kpbasis;

}

template <class KpointType> void Kpoint<KpointType>::set_pool(KpointType *pool)
{
    KpointType *tptr;
    int state;

    this->orbital_storage = pool;

    tptr = pool;

    Kstates = new State<KpointType> [nstates];

    for(state = 0;state < nstates;state++) {
        Kstates[state].set_storage(tptr); 
        tptr += pbasis;
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

    int state, pbasis;
    double t1;
    State<KpointType> *sp, *sp1;
    KpointType *tmp_orbital;

    tmp_orbital = new KpointType[pbasis];

    for(state = 0;state < nstates;state++)
    {
        sp = Kstates + state;
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

