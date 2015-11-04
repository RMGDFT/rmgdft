#ifndef RMG_vdW_H
#define RMG_vdW_H 1

#include "BaseGrid.h"

void xc_vdW_DF (BaseGrid &G, Lattice &L, TradeImages &T, int type, double *rho_valence, double *rho_core, double &etxc, double &vtxc, double *v);


#endif
