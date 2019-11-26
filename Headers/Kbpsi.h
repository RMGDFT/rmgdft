/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#ifndef Kbpsi_H
#define Kbpsi_H 1

#include <vector>
#include <iterator>
//  @name Kbpsi 
//  struct for each ion's nonlocal  <beta|psi>
typedef struct
{
    int send_size, recv_size;
    int send_to_pe, recv_from_pe;
    std::vector<int> send_ions;
    std::vector<int> recv_ions;
} KBPSI_COMM_INFO;

typedef struct
{
    

    int kbpsi_comm_loop;
    

    int max_send_size, max_recv_size;

    // vector has a length of kbpsi_comm_loop
    KBPSI_COMM_INFO *comm_info;

    int *num_orbital_thision;

// pair.first:  orbital index
// pair.second: list of kbpsi values for diffenent projectors. length = num_projectors.
   std::vector<double> *kbpsi_ion;
   std::vector<double> *kbpsi_res_ion;
   std::vector<int> *orbital_index;

} KBPSI;

extern KBPSI Kbpsi_str;
extern double *Kbpsi_mat;
extern double *Upsi_mat;
extern double *Upsi_mat_dist;
extern double *Upsi_mat_local;

// used for NEGF now, potentially can be used in ON when the Localized orbitals can be blocked. 
extern std::vector<double *> Kbpsi_mat_block;
#endif
