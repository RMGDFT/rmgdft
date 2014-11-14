
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <float.h>
#include <math.h>

#include <iostream>
#include <cstring>
#include <cmath>
#include "make_conf.h"

#if PLPLOT_LIBS
#include <plstream.h>

#include "Plots.h"

// Version that takes std::vector args.
void LinePlotLog10y(const char * filename, const char * xlegend, const char * ylegend, const char * title, std::vector<double>& x, std::vector<double>& y)
{
    int nsize = y.size();
    std::vector<double> log10y(nsize);

    for(int i = 0;i < nsize;i++) log10y[i] = log10(y[i]);

    LinePlot(filename, xlegend, ylegend, title, x, log10y);

}

// Version that takes double array args. If x is NULL then the index of each element will be used as the x-value.
void LinePlotLog10y(const char * filename, const char * xlegend, const char * ylegend, const char * title, double *x, double *y, int nsize)
{
    std::vector<double> yn(nsize);
    std::vector<double> xn;

    for(int i = 0;i < nsize;i++) yn[i] = y[i];
    if(x) {
        for(int i = 0;i < nsize;i++) xn.push_back((double)i);
    }

    LinePlotLog10y(filename, xlegend, ylegend, title, xn, yn);

}

// If x has no entries then the index of each element will be used as the x-value.
void LinePlot(const char * filename, const char * xlegend, const char * ylegend, const char * title, std::vector<double>& x, std::vector<double>& y)
{

    int nsize = y.size();
    if(nsize < 2) return;   // Not much point in plotting only 1 point

    PLFLT xmin = DBL_MAX, xmax = -DBL_MAX, ymin = DBL_MAX, ymax = -DBL_MAX;

    PLFLT *xp = new PLFLT[nsize]; 
    PLFLT *yp = new PLFLT[nsize]; 

    // Output file name
    plsetopt("o", filename);
    plsetopt("dev", "pngcairo");


    // Background color is white
    plscolbg(255, 255, 255);

    // set color 15 to black
    plscol0(15, 0, 0, 0);

    // Line width
    //plwidth( 2 );


    // Initialize plot library
    plinit();

    if(!x.size()) {
        for(int i = 0;i < nsize;i++) xp[i] = (PLFLT)i;
        xmin = 0.0;
        xmax = (double)(5*(nsize/5) + 5);
    }
    else {
        for(int i = 0;i < nsize;i++) xp[i] = (PLFLT)x[i];
    }
    for(int i = 0;i < nsize;i++) {
        yp[i] = (PLFLT)y[i];
        if(yp[i] < ymin) ymin = yp[i];
        if(yp[i] > ymax) ymax = yp[i];
    }
    

    // Create a labelled box to hold the plot.
    plcol0(15);
    plenv( xmin, xmax, ymin, ymax, 0, 0 );
    pllab( xlegend, ylegend, title);


    // Plot the data that was prepared above in red.
    plcol0(1);
    plline( nsize, xp, yp );

    // Close PLplot library
    plend();

    delete [] yp;
    delete [] xp; 
}

// If x has no entries then the index of each element will be used as the x-value.
void MultiLinePlot(const char * filename, const char * xlegend, const char * ylegend, const char * title, 
double *x, double *y, int nsize, int nlines)
{

    if(nsize < 2) return;   // Not much point in plotting only 1 point

    PLFLT xmin = DBL_MAX, xmax = -DBL_MAX, ymin = DBL_MAX, ymax = -DBL_MAX;

    PLFLT *xp = new PLFLT[nsize]; 
    PLFLT *yp = new PLFLT[nsize]; 

    // Output file name
    plsetopt("o", filename);
    plsetopt("dev", "pngcairo");


    // Background color is white
    plscolbg(255, 255, 255);

    // set color 15 to black
    plscol0(15, 0, 0, 0);

    // Line width
    //plwidth( 2 );

    // Initialize plot library
    plinit();

    for(int i = 0;i < nsize;i++) {
        xp[i] = (PLFLT)x[i];
        if(xp[i] < xmin) xmin = xp[i];
        if(xp[i] > xmax) xmax = xp[i];
    }
    for(int i = 0;i < nsize * nlines;i++) {
        if(y[i] < ymin) ymin = y[i];
        if(y[i] > ymax) ymax = y[i];
    }
    

    // Create a labelled box to hold the plot.
    plcol0(15);
    plenv( xmin, xmax, ymin, ymax, 0, 0 );
    pllab( xlegend, ylegend, title);


    // Plot the data that was prepared above in red.
    for(int i = 0; i < nlines; i++)
    {
        for(int j = 0; j < nsize; j++)
            yp[j] = (PLFLT)y[i * nsize + j];

        plcol0(1);
        plline( nsize, xp, yp );
    }

    // Close PLplot library
    plend();

    delete [] yp;
    delete [] xp; 
}


#endif
