#include "main.h"


int minimum (REAL *array)
{
	int idx, i;
        REAL min;

	idx=0; 
	min= *array;
	for (i=1; i< FP0_BASIS; i++)
		if (*(array+i)< min)
		{
			min=*(array+i);
			idx=i;
		}
	return idx;
	
}
