#ifndef RMG_FiniteDiff_H
#define RMG_FiniteDiff_H 1

/* Order of finite differencing for driver routines */
#define APP_CI_FOURTH 4
#define APP_CI_SIXTH 6

#ifdef __cplusplus
#include "Lattice.h"
#include "TradeImages.h"

template <typename RmgType>
void CPP_app_cir_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, int order);
template <typename RmgType>
double CPP_app_cil_driver (Lattice *L, TradeImages *T, RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, int order);


#include "rmg_error.h"

class FiniteDiff {

private:
    Lattice *L;

public:
    FiniteDiff(Lattice *lptr);

    bool check_anisotropy(double hx, double hy, double hz, double limit);

    template <typename RmgType>
    double app_cil_sixth (RmgType *rptr, RmgType *b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    void app_cir_sixth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz);

    template <typename RmgType>
    double app_del2c (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app6_del2 (RmgType * a, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    double app_cil_fourth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz);

    template <typename RmgType>
    void app_cir_fourth (RmgType * rptr, RmgType * b, int dimx, int dimy, int dimz);

    template <typename RmgType>
    void app_cir_fcc (RmgType * a, RmgType * b, int dimx, int dimy, int dimz);

    template <typename RmgType>
    void app_cir_bcc (RmgType * a, RmgType * b, int dimx, int dimy, int dimz);

    template <typename RmgType>
    void app_cir_hex (RmgType * a, RmgType * b, int dimx, int dimy, int dimz);


};
#endif

#endif
