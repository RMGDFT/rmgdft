/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include "main.h"

void init_img_topo (int ndims)
{
    int count, dims[ndims+1], periods[ndims+1];
    char buf[MAX_CHAR];

    /* only group masters partake of inter-image communications */
    if ( pct.thispe == 0 ) 
    {
        if ( ndims > 0 ) {
            for ( count = 0; get_data( "image_communicator_topology", &buf, LINE, NULL); count++ )
                sscanf ( buf, "%d %d", &dims[count], &periods[count] );
        }
        else
        {
            /* default to a 1D non-periodic topology */
            ndims = 1;
            dims[0] = pct.images;
            periods[0] = 0;
        }

    MPI_Cart_create ( pct.rmg_comm, ndims, dims, periods, 1, &pct.img_topo_comm );
    }
}
