/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/*  pmo.h: define structure and parameters for parallel matrix operation  */


struct parallel_matrix_operation
{

    int npe_energy;
    int ncol;
    int nrow;
    int mblock;
    int *ictxt;

    int *mxllda_cond;
    int *mxlocc_cond;
    int *mxllda_lead;
    int *mxlocc_lead;
    int *desc_lead;
    int *desc_cond;
    int *desc_cond_lead;
    int *desc_lead_cond;
 
    int myblacs;
        
    int *orb_index;
    int *diag_begin; /* start address of the diagonal blocks in Htri ...*/
    int *offdiag_begin; /* start address of the diagonal blocks in Htri ...*/
    int *lowoffdiag_begin;
    int ntot; /* size of distributed matrix Htri, ... */
    int ntot_low;
};
typedef struct  parallel_matrix_operation parallel_matrix_operation;

extern parallel_matrix_operation pmo;

void zero_lead_image(double *);
