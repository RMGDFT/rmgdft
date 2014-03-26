#include "Mgrid.h"
#include "common_prototypes.h"

using namespace std;


// C wrappers
extern "C" void mgrid_solv (rmg_double_t * v_mat, rmg_double_t * f_mat, rmg_double_t * work,
                 int dimx, int dimy, int dimz,
                 rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz,
                 int level, int *nb_ids, int max_levels, int *pre_cyc,
                 int *post_cyc, int mu_cyc, rmg_double_t step, rmg_double_t k,
                 int gxsize, int gysize, int gzsize,
                 int gxoffset, int gyoffset, int gzoffset,
                 int pxdim, int pydim, int pzdim, int boundary_flag)
{
    Lattice L;
    Mgrid MG(&L);
    MG.mgrid_solv<double>( v_mat, f_mat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz,
                   level, nb_ids, max_levels, pre_cyc, post_cyc, mu_cyc, step, k,
                   gxsize, gysize, gzsize,
                   gxoffset, gyoffset, gzoffset,
                   pxdim, pydim, pzdim, boundary_flag);

}

extern "C" void mg_restrict (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{
    Lattice L;
    Mgrid MG(&L);
    MG.mg_restrict<double>(full, half, dimx, dimy, dimz, dx2, dy2, dz2, xoffset, yoffset, zoffset);
}

extern "C" void mg_prolong (rmg_double_t * full, rmg_double_t * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset)
{
    Lattice L;
    Mgrid MG(&L);
    MG.mg_prolong<double>(full, half, dimx, dimy, dimz, dx2, dy2, dz2, xoffset, yoffset, zoffset);
}

extern "C" void eval_residual (rmg_double_t * mat, rmg_double_t * f_mat, int dimx, int dimy, int dimz,
                    rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t * res)
{
    Lattice L;
    Mgrid MG(&L);
    MG.eval_residual<double>(mat, f_mat, dimx, dimy, dimz, gridhx, gridhy, gridhz, res);
}

extern "C" int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype)
{
    Lattice L;
    Mgrid MG(&L);
    return MG.MG_SIZE(curdim, curlevel, global_dim, global_offset, global_pdim, roffset, bctype);
}

extern "C" void solv_pois (rmg_double_t * vmat, rmg_double_t * fmat, rmg_double_t * work,
                int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t step, rmg_double_t k)
{
    Lattice L;
    Mgrid MG(&L);
    MG.solv_pois<double>(vmat, fmat, work, dimx, dimy, dimz, gridhx, gridhy, gridhz, step, k);
}
