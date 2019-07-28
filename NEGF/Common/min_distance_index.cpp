#include "negf_prototypes.h"
/*
 **    $Id$    **
******************************************************************************/
 

/*
 *  find the which distance is the minimum and return the index
 */




int min_distance_index(double *distance, int n )
{

    int index, i;
    index = 0;

    double min;
    min = distance[0];
    for (i = 1; i< n; i++)
    {
        if(min > distance[i]) 
        {
            index = i;
            min = distance[i];
        }
    }
    return index;

}

