#ifndef RMG_ION_H
#define RMG_ION_H 1

#include "species.h"
#include "const.h"

/* Ion structure */
class ION
{

public:

    void RotateCoordinates(void)
    {
        ocrds3[0] = ocrds2[0];
        ocrds3[1] = ocrds2[1];
        ocrds3[2] = ocrds2[2];

        ocrds2[0] = ocrds1[0];
        ocrds2[1] = ocrds1[1];
        ocrds2[2] = ocrds1[2];

        ocrds1[0] = crds[0];
        ocrds1[1] = crds[1];
        ocrds1[2] = crds[2];
    }

    void RotateForces(void)
    {
        for(int i = 3;i > 0;i--)
        {
            force[i][0] = force[i-1][0];
            force[i][1] = force[i-1][1];
            force[i][2] = force[i-1][2];
        }
    }

    void ZeroForces(void)
    {
        for (int ic = 0; ic < 4; ic++)
        {
            force[ic][0] = 0.0;
            force[ic][1] = 0.0;
            force[ic][2] = 0.0;
        }
    }

    void ZeroVelocity(void)
    {
        velocity[0] = 0.0;
        velocity[1] = 0.0;
        velocity[2] = 0.0;
    }

    double GetKineticEnergy(void)
    {
        /* Get ionic mass */
        double mass = this->Type->atomic_mass * mu_me;
        double ke = 0.0;
        if (movable[0])
            ke += velocity[0] * velocity[0] * mass;
        if (movable[1])
            ke += velocity[1] * velocity[1] * mass;
        if (movable[2])
            ke += velocity[2] * velocity[2] * mass;
        return 0.5 * ke;
    }

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

    /* Pointer to species class. Would prefer a reference to a pointer but
     * Species is setup after Atoms so can't use a reference right now. Fix later.
     */
    SPECIES *Type;

    /* Forces on the ion */
    double force[4][3];

    /* Current velocity of the ion */
    double velocity[3];

    /* Milliken normalization coefficients */
    double mnorm[(MAX_L + 1) * (MAX_L + 1)];

    /* Total number of projectors */
    int prjcount;

    /* Movable flag */
    int movable[3];

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
    double init_spin_rho;
    double init_spin_x;
    double init_spin_y;
    double init_spin_z;

    // Augmentation charges associated with this ion
    std::vector<double> augfunc;
    std::vector<double> augfunc_xyz[3];

    // An index array which maps the q-functions onto the 3-d grid associated with each processor.
    std::vector<int> Qindex;

};


#endif

