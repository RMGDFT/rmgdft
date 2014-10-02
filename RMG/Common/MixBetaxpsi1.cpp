#include <stdlib.h>
#include <complex>
#include "transition.h"
#include "const.h"
#include "RmgTimer.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "../Headers/prototypes.h"


template void MixBetaxpsi1(State<double> *);
template void MixBetaxpsi1(State<std::complex<double > > *);

// Mixes projectors for a single orbital
// using a power function.
template <typename OrbitalType>
void MixBetaxpsi1 (State<OrbitalType> *sp)
{

    int size, koffset, loffset, idx1, kpt, istate;
    double scale;

    kpt = sp->kidx;
    istate = sp->istate;

    scale = pow(1.0 - ct.prjmix, (double)istate);
    if(istate == 0) scale = 1.0 - ct.prjmix;
    size = ct.max_nl;
    koffset = kpt * pct.num_nonloc_ions * ct.num_states * ct.max_nl;

    for (idx1 = 0; idx1 < pct.num_nonloc_ions; idx1++) {

        /* For localized <beta|psi>, there is offset due to both k-point and ion*/
        loffset = koffset + idx1 * ct.num_states * ct.max_nl;
        for(int ix = 0;ix < size;ix++) pct.oldsintR_local[loffset + istate * ct.max_nl + ix] *= scale;
//        my_scal( scale, &pct.oldsintR_local[loffset + istate * ct.max_nl], size);

        scale = 1.0 - scale;
        for(int ix = 0;ix < size;ix++) 
            pct.oldsintR_local[loffset + istate * ct.max_nl + ix] += scale * pct.newsintR_local[loffset + istate * ct.max_nl + ix];

//        my_axpy(1.0 - scale, &pct.newsintR_local[loffset + istate * ct.max_nl],
//                &pct.oldsintR_local[loffset + istate * ct.max_nl], size);


    }
}

