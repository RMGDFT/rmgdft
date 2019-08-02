#include "TradeImages.h"
#include "RmgSumAll.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "BlasWrappers.h"
#include "vhartree.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_error.h"
#include "transition.h"
#include "packfuncs.h"




// C wrappers
void mgrid_solv (double * v_mat, double * f_mat, double * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundary_flag)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.mgrid_solv<double>( v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz,
                   level, max_levels, pre_cyc, post_cyc, mu_cyc, step, 0.0, k, NULL,
                   gxsize, gysize, gzsize,
                   gxoffset, gyoffset, gzoffset,
                   pxdim, pydim, pzdim, boundary_flag);

}

void mg_restrict (double * full, double * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.mg_restrict<double>(full, half, dimx, dimy, dimz, dx2, dy2, dz2, xoffset, yoffset, zoffset);
}

void mg_prolong (double * full, double * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.mg_prolong<double>(full, half, dimx, dimy, dimz, dx2, dy2, dz2, xoffset, yoffset, zoffset);
}

void eval_residual (double * mat, double * f_mat, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, double * res)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.eval_residual<double>(mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, res, NULL);
}

int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    return MG.MG_SIZE(curdim, curlevel, global_dim, global_offset, global_pdim, roffset, bctype);
}

void solv_pois (double * vmat, double * fmat, double * work,
                int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double k)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.solv_pois<double>(vmat, fmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, 0.0, k, NULL);
}
