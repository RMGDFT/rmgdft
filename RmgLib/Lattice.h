#ifndef RMG_Lattice_H
#define RMG_Lattice_H 1

//#include "rmgtypes.h"
#include "rmg_error.h"


/* Crystal lattice types */
#define CUBIC_PRIMITIVE                 1
#define CUBIC_FC                        2
#define CUBIC_BC                        3
#define HEXAGONAL                       4
#define TRIGONAL_PRIMITIVE              5
#define TETRAGONAL_PRIMITIVE            6
#define TETRAGONAL_BC                   7
#define ORTHORHOMBIC_PRIMITIVE          8
#define ORTHORHOMBIC_BASE_CENTRED       9
#define ORTHORHOMBIC_BC                 10
#define ORTHORHOMBIC_FC                 11
#define MONOCLINIC_PRIMITIVE            12
#define MONOCLINIC_BASE_CENTRED         13
#define TRICLINIC_PRIMITIVE             14

class Lattice {

private:

    // Grid bravais lattice type 
    int ibrav;

    // lengths of the sides of the supercell
    double xside;
    double yside;
    double zside;

public:

    // lattice vectors
    double a0[3];
    double a1[3];
    double a2[3];

    // reciprocal lattice vectors
    double b0[3];
    double b1[3];
    double b2[3];

    // cell dimensions
    double celldm[6];

    // Total cell volume
    double omega;

    void latgen (double * celldm, double * OMEGAI, double *a0, double *a1, double *a2, int *flag);

    void cross_product (double * a, double * b, double * c);
    void to_crystal (double *crystal, double *cartesian);
    void to_cartesian (double *crystal, double *cartesian);
    void recips (void);
    int get_ibrav_type(void);
    void set_ibrav_type(int newtype);
    double metric (double * crystal);
    double get_omega(void);
    double get_celldm(int which);
    double get_a0(int which);
    double get_a1(int which);
    double get_a2(int which);
    double get_b0(int which);
    double get_b1(int which);
    double get_b2(int which);

    double get_xside(void);
    double get_yside(void);
    double get_zside(void);

};

#endif
