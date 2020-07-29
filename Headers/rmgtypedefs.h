
#ifndef RMG_TYPEDEFS_H
#define RMG_TYPEDEFS_H 1

#include "params.h"

#ifdef __cplusplus
    #include <complex>
    #include <boost/multi_array.hpp>
    typedef std::complex<double> DoubleC;
    typedef boost::multi_array<int, 1> int_1d_array;
    typedef boost::multi_array<int, 2> int_2d_array;
    typedef boost::multi_array<int, 3> int_3d_array;
    typedef boost::multi_array<int, 4> int_4d_array;
    typedef boost::multi_array<float, 1> float_1d_array;
    typedef boost::multi_array<float, 2> float_2d_array;
    typedef boost::multi_array<float, 3> float_3d_array;
    typedef boost::multi_array<float, 4> float_4d_array;
    typedef boost::multi_array<double, 1> double_1d_array;
    typedef boost::multi_array<double, 2> double_2d_array;
    typedef boost::multi_array<double, 3> double_3d_array;
    typedef boost::multi_array<double, 4> double_4d_array;
    typedef boost::multi_array<std::complex<double>, 1> doubleC_1d_array;
    typedef boost::multi_array<std::complex<double>, 2> doubleC_2d_array;
    typedef boost::multi_array<std::complex<double>, 3> doubleC_3d_array;
    typedef boost::multi_array<std::complex<double>, 4> doubleC_4d_array;
#else
    #include <complex.h>
    typedef complex double   DoubleC;
#endif

/*Structure for storing PDB information
 * Each ion should have it*/
typedef struct
{

/* 1 -  6  Record name*/
char record_name[7];

/* 7 - 11 Atom serial number*/
int serial_num;

/*13 - 16  Atom name*/
char name[5];

/* 17 Alternate location indicator.*/
char altLoc[2];

/* 18 - 20 Residue name*/
char resName[4];

/* 22 Chain identifier*/
char chainID[2];

/* 23 - 26 Residue sequence number*/
int resSeq;

/* 27 Code for insertion of residues*/
char iCode[2];

/* 55 - 60 Occupancy*/
double occupancy;

/* 61 - 66 Temperature factor*/
double tempFactor;

/* 77 - 78  Element symbol, right-justified. */
char element[3];

/*79 - 80  Charge on the atom.*/
char charge[3];

} PDB_INFO;



/* Ion structure */
#include "ION.h"

/* structure for TF ions (simple solvent model) */
typedef struct
{

    /* Actual Physical coordinates at current time step */
    double crds[3];


    /* Actual crystal coordinates at current time step */
    double xtal[3];

    /* q and alpha parameters for gaussian representing charge density and ions of TF molecules
     * Both are gaussians: q * pow( (alpha/PI), 1.5) * exp (-1.0*alpha*r*r) 
     * Gaussians representing ions should be sharper than for charge density
     * For water good starting point is: 
     * O q=7.05  alpha=0.75 for electrons and q0=6.0 alpha0= 1.5 for ions
     * H q=0.475 alpha=0.80 for electrons and q0=1.0 alpha0= 1.6 for ions */
    double q;
    
    /* alpha parameter for gaussian representing charge density*/
    /* q * pow( (alpha/PI), 1.5) * exp (-1.0*alpha*r*r) */
    double alpha;
    
    /* q parameter for gaussian representing ions*/
    /* q * pow( (alpha/PI), 1.5) * exp (-1.0*alpha*r*r) */
    double q0;
    
    /* alpha parameter for gaussian representing charge density*/
    /* q * pow( (alpha/PI), 1.5) * exp (-1.0*alpha*r*r) */
    double alpha0;
} TF_ION;

typedef struct {int species; int ip; int l; int m; int proj_index;} PROJ_INFO;

#endif

