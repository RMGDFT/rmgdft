#include "transition.h"
#include "common_prototypes.h"

/// C interface function
extern "C" void set_ibrav_type(int newtype)
{
  Rmg_L.set_ibrav_type(newtype);
}
/// C interface function
extern "C" int get_ibrav_type(void)
{
  return Rmg_L.get_ibrav_type();
}
extern "C" double get_xside(void)
{
    return Rmg_L.get_xside();
}
extern "C" double get_yside(void)
{
    return Rmg_L.get_yside();
}
extern "C" double get_zside(void)
{
    return Rmg_L.get_zside();
}
extern "C" double get_vel(void)
{
    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    return Rmg_L.get_omega() / t1;
}
extern "C" double get_omega(void)
{
    return Rmg_L.get_omega();
}
extern "C" double get_vel_f(void)
{
    double t1 = (double) (Rmg_G->get_NX_GRID(Rmg_G->get_default_FG_RATIO()) * Rmg_G->get_NY_GRID(Rmg_G->get_default_FG_RATIO()) * Rmg_G->get_NZ_GRID(Rmg_G->get_default_FG_RATIO()));
    return Rmg_L.get_omega() / t1;
}
extern "C" double get_celldm(int which)
{
    return Rmg_L.get_celldm(which);
}
extern "C" double get_a0(int which)
{
    return Rmg_L.get_a0(which);
}
extern "C" double get_a1(int which)
{
    return Rmg_L.get_a1(which);
}
extern "C" double get_a2(int which)
{
    return Rmg_L.get_a2(which);
}
extern "C" double get_b0(int which)
{
    return Rmg_L.get_b0(which);
}
extern "C" double get_b1(int which)
{
    return Rmg_L.get_b1(which);
}
extern "C" double get_b2(int which)
{
    return Rmg_L.get_b2(which);
}
extern "C" void to_crystal (double *crystal, double *cartesian)
{
    Rmg_L.to_crystal(crystal, cartesian);
}
extern "C" void to_cartesian (double *crystal, double *cartesian)
{
    Rmg_L.to_cartesian(crystal, cartesian);
}
extern "C" void cross_product (double * a, double * b, double * c)
{
    Rmg_L.cross_product(a, b, c);
}
extern "C" double metric (double * crystal)
{
    return Rmg_L.metric(crystal);
}
extern "C" void recips(void)
{
    Rmg_L.recips();
}
extern "C" void latgen (double *celldm, double *a0, double *a1, double *a2, double *OMEGAI, int *flag)
{
    Rmg_L.latgen(celldm, OMEGAI, a0, a1, a2, flag);
}
extern "C" void latgen_f_ (int *ibrav, double *celldm, double *a0, double *a1, double *a2, double *OMEGAI, int *flag)
{
    //*ibrav = Rmg_L.get_ibrav_type();
    Rmg_L.latgen(celldm, OMEGAI, a0, a1, a2, flag);
}
extern "C" double get_anisotropy(void)
{
    double hmaxgrid = get_xside() * get_hxgrid();
    if (get_yside() * get_hygrid() > hmaxgrid)
        hmaxgrid = get_yside() * get_hygrid();
    if (get_zside() * get_hzgrid() > hmaxgrid)
        hmaxgrid = get_zside() * get_hzgrid();

    double hmingrid = get_xside() * get_hxgrid();
    if (get_yside() * get_hygrid() < hmingrid)
        hmingrid = get_yside() * get_hygrid();
    if (get_zside() * get_hzgrid() < hmingrid)
        hmingrid = get_zside() * get_hzgrid();

    return hmaxgrid / hmingrid;
}
