
#include <Kpoint.h>
#include <complex>

template <class KpointType> Kpoint<KpointType>::Kpoint(double *kpt, double kweight, KpointType *pool, int nstates, int pbasis)
{

    int state;
    KpointType *tptr;

    State<KpointType> *sp;
  
    orbital_storage = pool;
    tptr = pool;

    Kstates = new State<KpointType> [nstates];

    sp = Kstates;

    for(state = 0;state < nstates;state++) {
        Kstates[state].set_storage(tptr); 
        tptr += pbasis;
    }

}

template <class KpointType> void Kpoint<KpointType>::sort_orbitals(void)
{

    int state, pbasis;
    double t1;
    State<KpointType> *sp, *sp1;
    KpointType *tmp_orbital;

    tmp_orbital = new KpointType[pbasis];

    for(state = 0;state < this.num_states;state++)
    {
        sp = this.Kstate + state;
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

void TestKpoint(void)
{

  float t1[20000];
  double kpt[3];
  Kpoint<float> K1(kpt, 0.5, t1, 100, 10000);

}
