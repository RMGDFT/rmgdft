#include "negf_prototypes.h"
/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* splint.c   ***
 * Interpolation routine ... spline.c is also required.
 * Input: xa[i], ya[i], y2a[i], n, x
 * Output: y (interpolated data for a given x)
 */


#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "main.h"
#include "init_var.h"
#include "LCR.h"



void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{

    int klo,khi,k;
    double h,b,a;

    klo=0;   
    khi=n-1;   
    while (khi-klo > 1) 
    {   
        k=(khi+klo) >> 1;  
        if (xa[k] > x) khi=k; 
        else klo=k;        
    }                      
    h=xa[khi]-xa[klo];
    if (h == 0.0) 
    {
        printf (" Bad xa input to routine splint \n");
        exit(0);
    }
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
   
}
