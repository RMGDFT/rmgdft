/************************** SVN Revision Information **************************
 **    $Id: InitPe4image.cpp 3142 2015-08-07 16:16:21Z luw $    **
 ******************************************************************************/
/*   init  processor topology for multiple images, spin/no spin, kpoint
 *   , and grid
 *   MPI_COMM_WORLD  :   first level
 *   pct.img_comm :      second level: number of proc: NPES=pct.npes_images[pct.thisimag]
 *                               it split into the following
 *                               3dimensional cartesian group
 *                       (ct.spin_flag +1) * pct.pe_kpoint *
 *                       numpeforgrid = npes of pct.img_comm
 *                 
 *   pct.spin_comm:      third level:  2 proc for spin-polarized, 1 for  nonspin
 *                           this one contains NPES/(ct.spin_flag+1)
 *                           communicators
 *                                  this will split into two car groups
 *  pct.kpsub_comm:      third level:   num_proc is pct.pe_kpoint, this comm should be used for sum over kpoint    
 *                                      it contains NPES/pct.pe_kpoint communicators 
 *  pct.grid_comm:       third level:   NPES/(ct.spinflag+1)/pct.pe_kpoint     
 *                                       must be divisible 
 *                                  this is the grid communicator, size
 *                                  should be pex * pey * pez from
 *                                  regular input file
 *     Wenchang Lu, NCSU, June, 2014
 */
    

#include "grid.h"
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "transition.h"


void InitPe4kpspin()
{

    // Initialize coalesce factors. Will adjust later
    for(int level = 0;level < MAX_MG_LEVELS;level++)
    {
        pct.coalesce_factors[level][0] = 1;
        pct.coalesce_factors[level][1] = 1;
        pct.coalesce_factors[level][2] = 1;
    }

 // Get the prime factors of kpt
    std::vector<int> kpt_factors = {1};
    GetPrimeFactors(kpt_factors, ct.num_kpts, ct.num_kpts);
    int nspin = 1+ct.spin_flag;
    if(pct.image_npes[pct.thisimg] % nspin !=0)
    {
        printf("\n npes %d for image needs to be even for spin-polarized case spin = %d", pct.image_npes[pct.thisimg], ct.spin_flag);
        exit(0);
    }

    int npes;
    npes = pct.image_npes[pct.thisimg]/nspin;

    if(pct.pe_kpoint <1 || pct.pe_kpoint > npes )
    {
        pct.pe_kpoint = 1;

        while(kpt_factors.size()) {
            if(npes % kpt_factors.back() == 0) {
                pct.pe_kpoint *= kpt_factors.back();
                npes /= kpt_factors.back();
            }
            kpt_factors.pop_back();
        }
    }

    npes = pct.image_npes[pct.thisimg]/nspin;

    if(npes % pct.pe_kpoint != 0) 
    {
        printf("\n npes %d for k-parallelization needs to be divisible by pe_kpoint %d", npes, pct.pe_kpoint);
        exit(0);

    }
    //pct.pe_kpoint = std::min(2, npes);
    //pct.pe_kpoint = 2;
    // set up communicator for spin
    int ndims = 3;

    int dims[] = { nspin, pct.pe_kpoint, pct.image_npes[pct.thisimg]/nspin/pct.pe_kpoint};

    int periods[] = { 0, 0, 0 };
    int reorder = 1;
    MPI_Cart_create (pct.img_comm, ndims, dims, periods, reorder, &pct.grid_topo_comm);

    int remains[] = { 1, 0, 0 };
    MPI_Cart_sub (pct.grid_topo_comm, remains, &pct.spin_comm);

    remains[0] = 0;
    remains[1] = 1;
    remains[2] = 0;
    MPI_Cart_sub (pct.grid_topo_comm, remains, &pct.kpsub_comm);
    /* set spinpe rank value to local spin rank value */
    MPI_Comm_rank (pct.spin_comm, &pct.spinpe);

    remains[0] = 0;
    remains[1] = 0;
    remains[2] = 1;
    MPI_Cart_sub (pct.grid_topo_comm, remains, &pct.grid_comm);

    MPI_Barrier(MPI_COMM_WORLD);
    /* set gridpe rank value to local grid rank value */
    int kpsub_rank;
    MPI_Comm_rank (pct.grid_comm, &pct.gridpe);
    MPI_Comm_rank (pct.kpsub_comm, &kpsub_rank);
    MPI_Comm_size (pct.grid_comm, &NPES);
    MPI_Comm_size (pct.grid_comm, &pct.grid_npes);

    ct.num_kpts_pe = ct.num_kpts / pct.pe_kpoint;
    pct.kstart = ct.num_kpts_pe * kpsub_rank;

    int kpt_mode = ct.num_kpts % pct.pe_kpoint;
    if( kpt_mode > 0) 
    {
        if(kpsub_rank < kpt_mode) 
        {
            ct.num_kpts_pe++;
            pct.kstart = ct.num_kpts_pe * kpsub_rank;
        }
        else
        {
            pct.kstart = ct.num_kpts_pe * kpsub_rank + kpt_mode;
        }
    }

    // Set up grids and neighbors
    //set_rank(pct.gridpe);

}                               /* end init_pe */


