/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/
 
/****f* spline.c   ***
 * This is a Cubic Spline interpolation routine (calculates derivative).
 * splint.c needs to call for completing the interpolation.
 * Input: x[i], y[i], n, yp1, ypn
 * Output: y2[i] 
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
#include "init_var_negf.h"
#include "LCR.h"



void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{


    int i,k;
    double p,qn,sig,un,*u;
/*
    for (i = 0; i < n; i++)
    {
    if(pct.gridpe ==0) printf (" hello  %d    %f    %f  \n", i, y[i], x[i] );
    }
*/

    my_malloc_init(u, n-1, double);   
    if (yp1 > 0.99e30)   
    {
        y2[0]=0.0; 
        u[0]=0.0;
    } 
    else 
    {                
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for (i=1;i<=n-2;i++) 
    { 
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0; 
        y2[i]=(sig-1.0)/p;
        u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))
                  /(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30) 
    {    
        qn=0.0;        
        un=0.0; 
    }       
    else 
    {               
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k=n-2;k>=0;k--) 
    { 
        y2[k]=y2[k]*y2[k+1]+u[k];     
    }
    my_free(u);
}
