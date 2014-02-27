#ifndef RMG_Lattice_H
#define RMG_Lattice_H 1

#include "const.h"
#include "rmgtypes.h"
#include "rmg_error.h"

class Lattice {

public:
    // lengths of the sides of the supercell
    static rmg_double_t xside;
    static rmg_double_t yside;
    static rmg_double_t zside;

    // lattice vectors
    static rmg_double_t a0[3];
    static rmg_double_t a1[3];
    static rmg_double_t a2[3];

    // reciprocal lattice vectors
    static rmg_double_t b0[3];
    static rmg_double_t b1[3];
    static rmg_double_t b2[3];

    // cell dimensions
    static rmg_double_t celldm[6];

    // Total cell volume
    static rmg_double_t omega;

    // Volume elements on coarse and fine grids
    static rmg_double_t vel;
    static rmg_double_t vel_f;

    // Global uniform grid spacing in x
    static rmg_double_t hxgrid;

    // Global uniform grid spacing in y
    static rmg_double_t hygrid;

    // Global uniform grid spacing in z
    static rmg_double_t hzgrid;

    // The fine uniform grid spacing in x
    static rmg_double_t hxxgrid;

    // The fine uniform grid spacing in y
    static rmg_double_t hyygrid;

    // The fine uniform grid spacing in z
    static rmg_double_t hzzgrid;

    void latgen (rmg_double_t * celldm, rmg_double_t * OMEGAI, int *flag);

    void cross_product (rmg_double_t * a, rmg_double_t * b, rmg_double_t * c);
    void to_crystal (rmg_double_t *crystal, rmg_double_t *cartesian);
    void to_cartesian (rmg_double_t *crystal, rmg_double_t *cartesian);
    void recips (void);
    rmg_double_t metric (rmg_double_t * crystal);

};

#endif
