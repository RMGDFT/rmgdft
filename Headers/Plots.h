#ifndef RMG_Plot_H
#define RMG_Plot_H 1

void LinePlot(const char * filename, const char * xlegend, const char * ylegend, const char * title, std::vector<double>& x, std::vector<double>& y);
void LinePlot(const char * filename, const char * xlegend, const char * ylegend, const char * title, double *x, double *y, int nsize);
void LinePlotLog10y(const char * filename, const char * xlegend, const char * ylegend, const char * title, std::vector<double>& x, std::vector<double>& y);
void LinePlotLog10y(const char * filename, const char * xlegend, const char * ylegend, const char * title, double *x, double *y, int nsize);

void MultiLinePlot(const char * filename, const char * xlegend, const char * ylegend, const char * title,
double *x, double *y, int nsize, int nlines);


#endif
