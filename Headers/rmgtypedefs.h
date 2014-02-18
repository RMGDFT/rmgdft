
#include "rmgtypes.h"
#include "params.h"

#ifndef RMG_TYPEDEFS_H
#define RMG_TYPEDEFS_H 1


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
rmg_double_t occupancy;

/* 61 - 66 Temperature factor*/
rmg_double_t tempFactor;

/* 77 - 78  Element symbol, right-justified. */
char element[3];

/*79 - 80  Charge on the atom.*/
char charge[3];

} PDB_INFO;



/* Ion structure */
typedef struct
{

    /* Initial physical coordinates at start of run */
    rmg_double_t icrds[3];

    /* Actual Physical coordinates at current time step */
    rmg_double_t crds[3];

    /* Positions at the previous time step */
    rmg_double_t ocrds1[3];

    /* Positions at  2dt back */
    rmg_double_t ocrds2[3];

    /* Positions at  3dt back */
    rmg_double_t ocrds3[3];

    /* Initial crystal coordinates at start of run */
    rmg_double_t ixtal[3];

    /* Actual crystal coordinates at current time step */
    rmg_double_t xtal[3];

    /* Crystal coordinates  at the previous time step */
    rmg_double_t oxtal[3];

    /*Position of ion relative to the middle of non-local box around the ion 
     *          * determined in get_nlop, AIget_cindex sets this up*/
    rmg_double_t nlcrds[3];


    /* Coordinates of the corner of the grid that the local */
    /* difference potential is nonzero on.                  */
    rmg_double_t lxcstart;
    rmg_double_t lycstart;
    rmg_double_t lzcstart;


    /* Coordinates of the corner of the grid that the non-local */
    /* potential is nonzero on.                                 */
    rmg_double_t nlxcstart;
    rmg_double_t nlycstart;
    rmg_double_t nlzcstart;


    /* Coordinates of the corner of the grid that the Qfunction */
    /* potential is nonzero on.                                 */
    rmg_double_t Qxcstart;
    rmg_double_t Qycstart;
    rmg_double_t Qzcstart;


    /* Integer species type when using a raw pseudopotential */
    int species;

    /* Forces on the ion */
    rmg_double_t force[4][3];

    /* Current velocity of the ion */
    rmg_double_t velocity[3];

    /* Kleinman-Bylander normalization coefficients */
    rmg_double_t pd[(MAX_L + 1) * (MAX_L + 1)];

    /* Milliken normalization coefficients */
    rmg_double_t mnorm[(MAX_L + 1) * (MAX_L + 1)];

    /* Total number of projectors */
    int prjcount;

    /* Movable flag */
    int movable;

        /* Force modifier parameters */
        struct {
                rmg_double_t setA_weight;
                rmg_double_t setA_coord[3];
                rmg_double_t setB_weight;
                rmg_double_t setB_coord[3];
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
    double *fftw_phase_sin;
    double *fftw_phase_cos;

    /*Stores PDB information*/
    PDB_INFO pdb;

} ION;



#endif

