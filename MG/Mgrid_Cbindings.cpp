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
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "transition.h"
#include "packfuncs.h"




// C wrappers
extern "C" void mgrid_solv (double * v_mat, double * f_mat, double * work,
                 int dimx, int dimy, int dimz,
                 double gridhx, double gridhy, double gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, double step, double k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundary_flag)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.mgrid_solv<double>( v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz,
                   level, nb_ids, max_levels, pre_cyc, post_cyc, mu_cyc, step, 0.0, k, NULL,
                   gxsize, gysize, gzsize,
                   gxoffset, gyoffset, gzoffset,
                   pxdim, pydim, pzdim, boundary_flag);

}

extern "C" void mg_restrict (double * full, double * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.mg_restrict<double>(full, half, dimx, dimy, dimz, dx2, dy2, dz2, xoffset, yoffset, zoffset);
}

extern "C" void mg_prolong (double * full, double * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.mg_prolong<double>(full, half, dimx, dimy, dimz, dx2, dy2, dz2, xoffset, yoffset, zoffset);
}

extern "C" void eval_residual (double * mat, double * f_mat, int dimx, int dimy, int dimz,
                    double gridhx, double gridhy, double gridhz, double * res)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.eval_residual<double>(mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, res, NULL);
}

extern "C" int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    return MG.MG_SIZE(curdim, curlevel, global_dim, global_offset, global_pdim, roffset, bctype);
}

extern "C" void solv_pois (double * vmat, double * fmat, double * work,
                int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double k)
{
    Mgrid MG(&Rmg_L, Rmg_T);
    MG.solv_pois<double>(vmat, fmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, 0.0, k, NULL);
}

extern "C" double get_vh (double * rho, double * rhoc, double * vh_eig, int min_sweeps, int max_sweeps, int maxlevel, double rms_target, int boundaryflag)
{
    int dimx = Rmg_G->get_PX0_GRID(Rmg_G->get_default_FG_RATIO()), dimy = Rmg_G->get_PY0_GRID(Rmg_G->get_default_FG_RATIO()), dimz = Rmg_G->get_PZ0_GRID(Rmg_G->get_default_FG_RATIO());
    int pbasis = dimx * dimy * dimz;
    double *rho_neutral = new double[pbasis];

    /* Subtract off compensating charges from rho */
    for (int i = 0; i < pbasis; i++)
        rho_neutral[i] = rho[i] - rhoc[i];

    // Make sure it is completely neutral
    double sum = 0.0;
    for(int i=0;i < pbasis;i++) sum += rho_neutral[i];
    sum = RmgSumAll(sum, pct.grid_comm) / (double)Rmg_G->get_GLOBAL_BASIS(Rmg_G->get_default_FG_RATIO());
    for(int i=0;i < pbasis;i++) rho_neutral[i] -= sum;

    double residual = CPP_get_vh (Rmg_G, &Rmg_L, Rmg_T, rho_neutral, ct.vh_ext, min_sweeps, max_sweeps, maxlevel, ct.poi_parm.gl_pre, 
                ct.poi_parm.gl_pst, ct.poi_parm.mucycles, rms_target, 
                ct.poi_parm.gl_step, ct.poi_parm.sb_step, boundaryflag, Rmg_G->get_default_FG_RATIO(), true);

    /* Pack the portion of the hartree potential used by the wavefunctions
     * back into the wavefunction hartree array. */
    CPP_pack_dtos (Rmg_G, vh_eig, ct.vh_ext, dimx, dimy, dimz, boundaryflag);

    delete [] rho_neutral;
    return residual;
}

