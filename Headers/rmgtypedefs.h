
#include "params.h"

#ifndef RMG_TYPEDEFS_H
#define RMG_TYPEDEFS_H 1

#ifdef __cplusplus
    #include <complex>
    typedef std::complex<double> DoubleC;
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
typedef struct
{

    /* Initial physical coordinates at start of run */
    double icrds[3];

    /* Actual Physical coordinates at current time step */
    double crds[3];

    /* Positions at the previous time step */
    double ocrds1[3];

    /* Positions at  2dt back */
    double ocrds2[3];

    /* Positions at  3dt back */
    double ocrds3[3];

    /* Initial crystal coordinates at start of run */
    double ixtal[3];

    /* Actual crystal coordinates at current time step */
    double xtal[3];

    /* Crystal coordinates  at the previous time step */
    double oxtal[3];

    /*Position of ion relative to the middle of non-local box around the ion 
     *          * determined in get_nlop, AIget_cindex sets this up*/
    double nlcrds[3];


    /* Coordinates of the corner of the grid that the local */
    /* difference potential is nonzero on.                  */
    double lxcstart;
    double lycstart;
    double lzcstart;


    /* Coordinates of the corner of the grid that the non-local */
    /* potential is nonzero on.                                 */
    double nlxcstart;
    double nlycstart;
    double nlzcstart;
 
    int nl_global_grid_xstart;
    int nl_global_grid_ystart;
    int nl_global_grid_zstart;


    /* Coordinates of the corner of the grid that the Qfunction */
    /* potential is nonzero on.                                 */
    double Qxcstart;
    double Qycstart;
    double Qzcstart;


    /* Integer species type when using a raw pseudopotential */
    int species;

    /* Forces on the ion */
    double force[4][3];

    /* Current velocity of the ion */
    double velocity[3];

    /* Kleinman-Bylander normalization coefficients */
    double pd[(MAX_L + 1) * (MAX_L + 1)];

    /* Milliken normalization coefficients */
    double mnorm[(MAX_L + 1) * (MAX_L + 1)];

    /* Total number of projectors */
    int prjcount;

    /* Movable flag */
    int movable;

        /* Force modifier parameters */
        struct {
                double setA_weight;
                double setA_coord[3];
                double setB_weight;
                double setB_coord[3];
        double forcemask[3];
        } constraint;

    /*  number of local orbitals on the ion */
    int n_loc_states;


    int ixstart;
    int iystart;
    int izstart;
    int ixend;
    int iyend;
    int izend;

    int frozen;
 
       /* Localization mask */
    char *lmask[4];

    

    int first_state_index;

    /*Stores PDB information*/
    PDB_INFO pdb;



    int ixstart_loc;
    int iystart_loc;
    int izstart_loc;
    int ixend_loc;
    int iyend_loc;
    int izend_loc;

    double xcstart_loc;
    double ycstart_loc;
    double zcstart_loc;

    int orbital_start_index;
    int num_orbitals;

    double partial_charge;


} ION;


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

#endif

