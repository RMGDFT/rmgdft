
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <unordered_map>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "rmg_error.h"
#include "transition.h"
#include "LaplacianCoeff.h"




void SetLaplacian()
{
    FiniteDiff FD(&Rmg_L);
    double a[3][3];
    int Ngrid[3], dim[3];
    int Lorder;
    Lorder = ct.kohn_sham_fd_order;
    for(int i = 0; i < 3; i++)
    {
        a[0][i] = Rmg_L.a0[i];
        a[1][i] = Rmg_L.a1[i];
        a[2][i] = Rmg_L.a2[i];
    }

    Ngrid[0] = Rmg_G->get_NX_GRID(1);
    Ngrid[1] = Rmg_G->get_NY_GRID(1);
    Ngrid[2] = Rmg_G->get_NZ_GRID(1);
    dim[0] = Rmg_G->get_PX0_GRID(1);
    dim[1] = Rmg_G->get_PY0_GRID(1);
    dim[2] = Rmg_G->get_PZ0_GRID(1);

    LC = new LaplacianCoeff(a, Ngrid, Lorder, dim);
    LC->SetBrav(Rmg_L.get_ibrav_type());
    LC->SetOffdiag(ct.laplacian_offdiag);
    LC->CalculateCoeff();
    FD.FdCoeffs.insert({dim[0]*dim[1]*dim[2]+8, LC});

    LC_6 = new LaplacianCoeff(a, Ngrid, 6, dim);
    LC_6->SetBrav(Rmg_L.get_ibrav_type());
    LC_6->SetOffdiag(ct.laplacian_offdiag);
    LC_6->CalculateCoeff();
    FD.FdCoeffs.insert({dim[0]*dim[1]*dim[2]+6, LC_6});

    LC_4 = new LaplacianCoeff(a, Ngrid, 4, dim);
    LC_4->SetBrav(Rmg_L.get_ibrav_type());
    LC_4->SetOffdiag(ct.laplacian_offdiag);
    LC_4->CalculateCoeff();
    FD.FdCoeffs.insert({dim[0]*dim[1]*dim[2]+4, LC_4});

    Ngrid[0] = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    Ngrid[1] = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    Ngrid[2] = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);
    dim[0] = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    dim[1] = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    dim[2] = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);

    HLC = new LaplacianCoeff(a, Ngrid, Lorder, dim);
    HLC->SetBrav(Rmg_L.get_ibrav_type());
    HLC->SetOffdiag(ct.laplacian_offdiag);
    HLC->CalculateCoeff();
    FD.FdCoeffs.insert({dim[0]*dim[1]*dim[2]+8, LC_4});

}

