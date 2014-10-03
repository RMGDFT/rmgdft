#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <complex>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "State.h"
#include "Kpoint.h"
#include "transition.h"


template void OutputBandPlot(Kpoint<double> **);
template void OutputBandPlot(Kpoint<std::complex<double> > **);

template <typename KpointType>
void OutputBandPlot(Kpoint<KpointType> ** Kptr)
{
    FILE *bs_f = fopen ("band.dat", "w");
    if(!bs_f) {
        rmg_printf("Unable to write band plot data.\n");
        return;
    }

    for (int is = 0; is < ct.num_states; is++)
    {
        for(int ik = 0; ik < ct.num_kpts; ik++)
        {

            fprintf (bs_f, "\n %4d  %16.8f ", ik, Kptr[ik]->Kstates[is].eig[0] * Ha_eV);

        }
        
        fprintf (bs_f, "\n &&");
    }

    fclose (bs_f);
}
