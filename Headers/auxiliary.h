#ifndef RMG_Auxiliary_H
#define RMG_Auxiliary_H 1

template <typename RmgType>
void CPP_genvpsi (RmgType * psi, RmgType * sg_twovpsi, rmg_double_t * vtot, rmg_double_t * vnl, rmg_double_t * kd,
              rmg_double_t kmag, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_pack_stop_axpy (RmgType * sg, RmgType * pg, rmg_double_t alpha, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_pack_stop (RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_pack_ptos(RmgType * sg, RmgType * pg, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_app_smooth (RmgType * f, RmgType * work, int dimx, int dimy, int dimz);

template <typename RmgType>
void CPP_app_smooth1 (RmgType * f, RmgType * work, int dimx, int dimy, int dimz);

void CPP_pack_dtos (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz, int boundaryflag);

void CPP_pack_stod (rmg_double_t * s, rmg_double_t * d, int dimx, int dimy, int dimz, int boundaryflag);

#endif
