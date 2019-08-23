#include "transition.h"
#include "common_prototypes.h"

/// C interface function
void set_ibrav_type(int newtype)
{
  Rmg_L.set_ibrav_type(newtype);
}
/// C interface function
int get_ibrav_type(void)
{
  return Rmg_L.get_ibrav_type();
}
double get_xside(void)
{
    return Rmg_L.get_xside();
}
double get_yside(void)
{
    return Rmg_L.get_yside();
}
double get_zside(void)
{
    return Rmg_L.get_zside();
}
double get_vel(void)
{
    double t1 = (double) (Rmg_G->get_NX_GRID(1) * Rmg_G->get_NY_GRID(1) * Rmg_G->get_NZ_GRID(1));
    return Rmg_L.get_omega() / t1;
}
double get_omega(void)
{
    return Rmg_L.get_omega();
}
double get_vel_f(void)
{
    double t1 = (double) (Rmg_G->get_NX_GRID(Rmg_G->get_default_FG_RATIO()) * Rmg_G->get_NY_GRID(Rmg_G->get_default_FG_RATIO()) * Rmg_G->get_NZ_GRID(Rmg_G->get_default_FG_RATIO()));
    return Rmg_L.get_omega() / t1;
}
double get_celldm(int which)
{
    return Rmg_L.get_celldm(which);
}
double get_a0(int which)
{
    return Rmg_L.get_a0(which);
}
double get_a1(int which)
{
    return Rmg_L.get_a1(which);
}
double get_a2(int which)
{
    return Rmg_L.get_a2(which);
}
double get_b0(int which)
{
    return Rmg_L.get_b0(which);
}
double get_b1(int which)
{
    return Rmg_L.get_b1(which);
}
double get_b2(int which)
{
    return Rmg_L.get_b2(which);
}
void to_crystal (double *crystal, double *cartesian)
{
    Rmg_L.to_crystal(crystal, cartesian);
}
void to_cartesian (double *crystal, double *cartesian)
{
    Rmg_L.to_cartesian(crystal, cartesian);
}
void cross_product (double * a, double * b, double * c)
{
    Rmg_L.cross_product(a, b, c);
}
double metric (double * crystal)
{
    return Rmg_L.metric(crystal);
}
void recips(void)
{
    Rmg_L.recips();
}
double get_anisotropy(void)
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
