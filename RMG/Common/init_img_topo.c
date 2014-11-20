/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/



#include "main.h"

void init_img_topo (int ndims)
{
    int count, dims[ndims+1], periods[ndims+1];
    char buf[MAX_CHAR];

    /* only group masters partake of inter-image communications */
    if ( pct.imgpe == 0 ) 
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
