#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cfloat>
#include <climits>
#include <set>
#include <list>
#include <vector>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


#include <boost/algorithm/string.hpp>

/* electron charge in Coulomb */
#define     e_C         1.60217646e-19

// Boltzmann constant in unit of eV/K
#define  Kb  8.617332478e-5; 

/* Planck constant in SI */
#define     h_SI        6.626068e-34

int main (int argc, char **argv)
{
    std::string line;

    std::string whitespace_delims = " \n\t\"";
    std::vector<double> ener, cond;
    std::vector<std::string> line_values;
    double bias = 0.10, temperature = 300.0;
   
    double KT ;
    KT = temperature * Kb;
    int E_POINTS = 0;
    
    int amode = S_IREAD | S_IWRITE;
    int item;


    std::ifstream file_cond (argv[1]);

    int line_end = 1;
    while(1){
        getline(file_cond, line); //read number
        boost::trim(line);
        boost::split( line_values, line, boost::is_any_of(whitespace_delims), boost::token_compress_on );
        if(line_values.size() != 2)
        {
            break;
        }
        E_POINTS ++;
    //    std::cout << E_POINTS <<line_values.size()<< line +'\n';
        ener.push_back(std::atof(line_values[0].c_str()) );
        cond.push_back(std::atof(line_values[1].c_str()) );
    }

    //printf("\n %d ", E_POINTS);
    //for (int i = 0; i < E_POINTS; i++)
     //   printf("\n %f %f ", ener[i], cond[i]);

    printf ("&& gate    current(microA)\n");
    double f1, f2;
    for(int igate = 0; igate < E_POINTS; igate++)
    {
        double gate = ener[igate];
        double EF1 = gate - 0.5 * bias;
        double EF2 = gate + 0.5 * bias;

        double current = 0.0;
        f1 = 1.0 / (1.0 + exp ((ener[0]-EF1) / KT));
        f2 = 1.0 / (1.0 + exp ((ener[0]-EF2) / KT));
        if( f2 - f1 > 1.0e-2) continue;

        f1 = 1.0 / (1.0 + exp ((ener[E_POINTS-1]-EF1) / KT));
        f2 = 1.0 / (1.0 + exp ((ener[E_POINTS-1]-EF2) / KT));
        if( f2 - f1 > 1.0e-2) continue;

        for (int iene = 1; iene < E_POINTS; iene++)
        {
            f1 = 1.0 / (1.0 + exp ((ener[iene]-EF1) / KT));
            f2 = 1.0 / (1.0 + exp ((ener[iene]-EF2) / KT));
            current += (f2 - f1) * cond[iene] * (ener[iene] - ener[iene - 1]);
        }


        current *= 2.0 * e_C * e_C / h_SI * 1e6;

       printf("%e      %e\n",  gate, current);

    }



}
