/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 

#ifndef Kbpsi_H
#define Kbpsi_H 1

//  @name Kbpsi 
//  struct for each ion's nonlocal  <beta|psi>
typedef struct
{
    
//  original ion index
    int ion_index;  

// num of projectors (beta) for this ion
    int num_projectors;

// num of orbitals overlaped with this ion;
    int num_overlap_orbitals;  

// num of orbitals on this pe overlaped with this ion;
    int num_overlap_orbitals_thispe;  

    int *orbital_index;     
// which orbital overlaps with this ion. 
//dimenstion of num_overlap_orbitals

    int *orbital_index_thispe;     
// which orbital in this pe overlaps with this ion. 
//dimenstion of num_overlap_orbitals_thispe

//kbpsi values <beta|phi>
// dimension of num_projectors * num_overlap_orbitals.
    double *kbpsi_value; 

// dimension of num_projectors * num_overlap_orbitals_thispe.
    double *kbpsi_value_thispe; 

} KBPSI;

#endif
