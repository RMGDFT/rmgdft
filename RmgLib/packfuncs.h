#ifndef RMG_packfuncs_H
#define RMG_packfuncs_H 1

#include "BaseGrid.h"

template <typename RmgType>
void CPP_pack_stop_axpy (RmgType * sg, RmgType * pg, double alpha, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_pack_stop (RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_pack_ptos(RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_app_smooth1 (RmgType * f, RmgType * work, int dimx, int dimy, int dimz);

void CPP_pack_dtos (BaseGrid *G, double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag);

void CPP_pack_stod (BaseGrid *G, double * s, double * d, int dimx, int dimy, int dimz, int boundaryflag);

#endif
