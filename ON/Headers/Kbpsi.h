/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#ifndef Kbpsi_H
#define Kbpsi_H 1

#include <vector>
#include <list>
//  @name Kbpsi 
//  struct for each ion's nonlocal  <beta|psi>
typedef struct
{
    long int send_size, recv_size;
    int send_to_pe, recv_from_pe;
    std::vector<int> send_ions;
    std::vector<int> recv_ions;
} KBPSI_COMM_INFO;

typedef struct
{
    int orbital_index;
    std::vector<double> kbpsi;
} KBPSI_ION;
    
typedef struct
{
    

    int kbpsi_comm_loop;

    // vector has a length of kbpsi_comm_loop
    KBPSI_COMM_INFO *comm_info;

// pair.first:  orbital index
// pair.second: list of kbpsi values for diffenent projectors. length = num_projectors.
   std::vector<KBPSI_ION> *kbpsi_ion;

} KBPSI;

extern KBPSI Kbpsi_str;
#endif
