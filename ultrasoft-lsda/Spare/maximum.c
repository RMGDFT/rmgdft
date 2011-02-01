#include "main.h"


int maximum (REAL *array)
{
	int idx, i;
        REAL max;

	idx=0; 
	max= *array;
	for (i=1; i< FP0_BASIS; i++)
		if (*(array+i) > max)
		{
			max = *(array+i);
			idx = i;
		}
	return idx;
	
}
