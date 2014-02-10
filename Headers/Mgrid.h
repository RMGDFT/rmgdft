#ifndef RMG_Mgrid_H
#define RMG_Mgrid_H 1

#include "Grid.h"

class Mgrid {

public:
    template <typename RmgType> void mg_restrict (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

    template <typename RmgType> void mg_prolong (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

template <typename RmgType> void eval_residual (RmgType * mat, RmgType * f_mat, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, RmgType * res);

template <typename RmgType> void solv_pois (RmgType * vmat, RmgType * fmat, RmgType * work,
                int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz, rmg_double_t step, rmg_double_t k);



};

#endif

