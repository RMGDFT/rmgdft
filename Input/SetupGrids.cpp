#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cfloat>
#include <climits>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include "BaseGrid.h"
#include "transition.h"
#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "CheckValue.h"
#include "RmgException.h"
#include "RmgInputFile.h"
#include "InputOpts.h"
#include "grid.h"



// Sets up both processor grid and wavefunction grid at the same time
void SetupGrids(int npes, int& NX_GRID, int& NY_GRID, int &NZ_GRID, double *celldm, double h, PE_CONTROL& pelc)
{

 
    int P_GRID[3] = {1, 1, 1};

    double AX = celldm[0];
    double AY = celldm[1];
    double AZ = celldm[2];

    // Get the prime factors of npes
    std::vector<int> npe_factors = {1};
    GetPrimeFactors(npe_factors, npes, npes);

    int iptr = 0;  // used to rotate over the three coordinate dimensions
    while(npe_factors.size()) {

        // Get the estimated values of NX_GRID, NY_GRID and NZ_GRID 
        // may need some adjustments for hexagonal grids
        std::fesetround(FE_TONEAREST);
        NX_GRID = std::rint(AX / h);
        NY_GRID = std::rint(AY / h);
        NZ_GRID = std::rint(AZ / h);

        // Adjust to make sure they are divisible by 4
        int ix = NX_GRID % 4;
        int iy = NY_GRID % 4;
        int iz = NZ_GRID % 4;
        NX_GRID = (ix <= 2) ? NX_GRID - ix: NX_GRID + 4 - ix;
        NY_GRID = (iy <= 2) ? NY_GRID - iy: NY_GRID + 4 - iy;
        NZ_GRID = (iz <= 2) ? NZ_GRID - iz: NZ_GRID + 4 - iz;

        // Find all of the prime factors of each
        std::vector<int> n_factors[3] = {{1},{1},{1}};
        GetPrimeFactors(n_factors[0], NX_GRID, NX_GRID);
        GetPrimeFactors(n_factors[1], NY_GRID, NY_GRID);
        GetPrimeFactors(n_factors[2], NZ_GRID, NZ_GRID);

        // If the largest prime factor of the dim references by iptr is also a prime
        // factor of npes then make that the processor x-grid dimension
        int osize = npe_factors.size();
        for(int i = 0;i < 3;i++) {
            if(n_factors[iptr].back() == npe_factors.back()) {
                P_GRID[iptr] *= npe_factors.back();
                npe_factors.pop_back();
                n_factors[iptr].pop_back();
                iptr++;
                if(iptr == 3) iptr = 0;
            }
        }

        // If nothing matched then absorb the prime factor into whichever dim
        // iptr references
        if(osize == npe_factors.size()) {
            P_GRID[iptr] *= npe_factors.back();
            npe_factors.pop_back();
            n_factors[iptr].pop_back();
            iptr++;
            if(iptr == 3) iptr = 0;
        }

        pelc.pe_x = P_GRID[0];
        pelc.pe_y = P_GRID[1];
        pelc.pe_z = P_GRID[2];
 
    }

}


// Sets up processor grid if wavefunction grid has already been specified
void SetupProcessorGrid(int npes, int NX_GRID, int NY_GRID, int NZ_GRID, PE_CONTROL& pelc)
{
    // We take this path if wavefunction grid is specified but processor grid is not
    std::vector<int> npe_factors = {1};
    std::vector<int> nx_factors = {1};
    std::vector<int> ny_factors = {1};
    std::vector<int> nz_factors = {1};
    GetPrimeFactors(npe_factors, npes, npes);
    GetPrimeFactors(nx_factors, NX_GRID, NX_GRID);
    GetPrimeFactors(ny_factors, NY_GRID, NY_GRID);
    GetPrimeFactors(nz_factors, NZ_GRID, NZ_GRID);
    pelc.pe_x = 1;
    pelc.pe_y = 1;
    pelc.pe_z = 1;
    int token = 0;
    npe_factors.erase(npe_factors.begin());
    std::reverse(npe_factors.begin(),npe_factors.end()); 
    for(auto it = npe_factors.begin();it != npe_factors.end(); ++it) {
        std::vector<int>::iterator nx_it, ny_it, nz_it;
        nx_it = std::find(nx_factors.begin(), nx_factors.end(), *it);
        ny_it = std::find(ny_factors.begin(), ny_factors.end(), *it);
        nz_it = std::find(nz_factors.begin(), nz_factors.end(), *it);
        if(token == 0) {
            if(nx_it != nx_factors.end()) {
                pelc.pe_x *= *it;
                nx_factors.erase(nx_it);
            }
            else if(ny_it != ny_factors.end()) {
                pelc.pe_y *= *it;
                ny_factors.erase(ny_it);
            }
            else if(nz_it != nz_factors.end()) {
                pelc.pe_z *= *it;
                nz_factors.erase(nz_it);
            }
            else {
                pelc.pe_x *= *it;
            }
        }
        else if(token == 1) {
            if(ny_it != ny_factors.end()) {
                pelc.pe_y *= *it;
                ny_factors.erase(ny_it);
            }
            else if(nz_it != nz_factors.end()) {
                pelc.pe_z *= *it;
                nz_factors.erase(nz_it);
            }
            else if(nx_it != nx_factors.end()) {
                pelc.pe_x *= *it;
                nx_factors.erase(nx_it);
            }
            else {
                pelc.pe_y *= *it;
            }
        }
        else if(token == 2) {
            if(nz_it != nz_factors.end()) {
                pelc.pe_z *= *it;
                nz_factors.erase(nz_it);
            }
            else if(nx_it != nx_factors.end()) {
                pelc.pe_x *= *it;
                nx_factors.erase(nx_it);
            }
            else if(ny_it != ny_factors.end()) {
                pelc.pe_y *= *it;
                ny_factors.erase(ny_it);
            }
            else {
                pelc.pe_z *= *it;
            }
        }
        token++;
        if(token == 3) token = 0;
    }
    if(pct.imgpe == 0) std::cout << "Auto processor grid settings: npes=" <<  npes << " PE_X=" << pelc.pe_x << " PE_Y=" << pelc.pe_y << " PE_Z=" << pelc.pe_z << std::endl;

}


// Sets up the wavefunction grid assuming the processor grid was specified
void SetupWavefunctionGrid(int npes, int& NX_GRID, int& NY_GRID, int &NZ_GRID, double *celldm, double h)
{

    // Get the estimated values of NX_GRID, NY_GRID and NZ_GRID 
    // may need some adjustments for hexagonal grids
    std::fesetround(FE_TONEAREST);
    NX_GRID = std::rint(celldm[0] / h);
    NY_GRID = std::rint(celldm[1] / h);
    NZ_GRID = std::rint(celldm[2] / h);

    // Adjust to make sure they are divisible by 4
    int ix = NX_GRID % 4;
    int iy = NY_GRID % 4;
    int iz = NZ_GRID % 4;
    NX_GRID = (ix <= 2) ? NX_GRID - ix: NX_GRID + 4 - ix;
    NY_GRID = (iy <= 2) ? NY_GRID - iy: NY_GRID + 4 - iy;
    NZ_GRID = (iz <= 2) ? NZ_GRID - iz: NZ_GRID + 4 - iz;

}

