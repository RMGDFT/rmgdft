#ifndef RMG_FiniteDiff_H
#define RMG_FiniteDiff_H 1

template <typename RmgType>
rmg_double_t FD_app_cil_sixth_standard (RmgType *rptr, RmgType *b, int dimx, int dimy, int dimz, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
template <typename RmgType>
rmg_double_t FD_app_cil_sixth_global (RmgType * rptr, RmgType * b, rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);
template <typename RmgType>
void FD_app_cir_sixth_standard (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz);
template <typename RmgType>
void FD_app_cir_sixth_global (RmgType * rptr, RmgType * b);
template <typename RmgType>
rmg_double_t FD_app_del2c (RmgType * a, RmgType * b, int dimx, int dimy, int dimz,
                rmg_double_t gridhx, rmg_double_t gridhy, rmg_double_t gridhz);



#endif
