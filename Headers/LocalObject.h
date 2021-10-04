/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#ifndef LocalObject_H
#define LocalObject_H 1

//  localized objects, such orbitals, projectors, or atomic wavefunctions, can distributed
//  in two methods. One way is to project the localized object onto the whole real space 
//  with processor domain decomposition. Another way is that localized objects are distributed
// over proc.

#include <vector>
#include <iterator>
#include "BaseGrid.h"

template <typename KpointType> class LocalObject {

public:
//    LocalObject(int num_objects, double *center_crd, double *radius, 
//                BaseGrid *Rmg_G, int coarse_fine, MPI_Comm comm);
    LocalObject(int num_objects, int *ixmin, int *iymin, int *izmin, 
                int *dimx, int *dimy, int *dimz, bool delocalized,
                BaseGrid &Rmg_G, int density, MPI_Comm comm);
    LocalObject(const LocalObject &);
    ~LocalObject(void);

    size_t storage_size;
    KpointType *storage_cpu;
    KpointType *storage_ptr;
    KpointType *storage_gpu;
    char *mask;

    int num_thispe;
    int num_tot;
    int density;
    int *index_proj_to_global;
    int *index_global_to_proj;
    bool delocalized;
    MPI_Comm comm;
    int pbasis;

    //OrbitalsOwnedProce[st][pe_list][pe, offset_x, offset_y, offset_z]
    //which processor has the data for a orbital, had its offset in x,y, z index.]
    std::vector<int> *OrbitalOwnedProcs_pelist;
    std::vector<int> *OrbitalOwnedProcs_xoffset;
    std::vector<int> *OrbitalOwnedProcs_yoffset;
    std::vector<int> *OrbitalOwnedProcs_zoffset;

    // Type LOCALIZED or DELOCALIZED
    int type;
    void ReadOrbitalsFromSingleFiles(std::string filename, BaseGrid &Rmg_G);
    void WriteOrbitals(std::string filename, BaseGrid &Rmg_G);
    void ReadProjectedOrbitals(std::string filename, BaseGrid &Rmg_G);
    void ReadProjectors(int num_ions, int max_nlpoint, int *num_proj_perion, BaseGrid &Rmg_G);
    void GetAtomicOrbitals(int num_ions, BaseGrid &Rmg_G);
    void SetZeroBoundary(BaseGrid &Rmg_G, int multi_grid_level, int fd_order);
    void ReAssign(BaseGrid &BG);
    void AssignOrbital(int st, KpointType *psi);
    void Normalize();
    void SetOrbitalProjToSingle(BaseGrid &);
    void WriteOrbitalsToSingleFiles(std::string filename, BaseGrid &Rmg_G);



private:
    int *ixmin, *iymin, *izmin, *dimx, *dimy, *dimz;
    int max_dimx, max_dimy, max_dimz;

};

extern LocalObject<double> *LocalOrbital;
extern LocalObject<double> *H_LocalOrbital;
extern LocalObject<double> *LocalProj;
extern LocalObject<double> *LocalAtomicOrbital;

#endif
