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



#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;
#include <iostream> 
#include <fstream>
#include <sstream>
#include <iterator>
#include <string> 
#include <cfloat> 
#include <climits> 
#include <unordered_map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include "portability.h"
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


/**********************************************************************


**********************************************************************/


namespace Ri = RmgInput;

void AutoSet(CONTROL& lc, PE_CONTROL& pelc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    static Ri::ReadVector<int> ProcessorGrid;
    static Ri::ReadVector<int> WavefunctionGrid;


    double *celldm = Rmg_L.celldm;
    static double grid_spacing;

    InputKey *ik;
    ik = InputMap["grid_spacing"];
    grid_spacing = ik->Readdoubleval[0];

    ik->Readdoubleval = &grid_spacing;
    ik = InputMap["wavefunction_grid"];
    WavefunctionGrid = ik->Vint;

    int NX_GRID = WavefunctionGrid.vals.at(0);
    int NY_GRID = WavefunctionGrid.vals.at(1);
    int NZ_GRID = WavefunctionGrid.vals.at(2);

    int ReqNPES = pelc.pe_x * pelc.pe_y * pelc.pe_z; 
    bool autoset_processor_grid = ((NPES != ReqNPES));

    // if this is equal to 1 then user did not try to manually set the
    // wavefunction grid so we have to set it ourselves which this test
    // will determine.
    bool autoset_wavefunction_grid = (1 == (NX_GRID * NY_GRID * NZ_GRID));

    // We take this path if wavefunction grid is specified but processor grid is not
    if(autoset_processor_grid && (!autoset_wavefunction_grid)) {

        SetupProcessorGrid(NPES, NX_GRID, NY_GRID, NZ_GRID, pelc);

    }
    else if(autoset_processor_grid && autoset_wavefunction_grid) {

        // Neither the wavefunction grid or the processor grid was set so we do them both here.
        SetupGrids(NPES, NX_GRID, NY_GRID, NZ_GRID, celldm, grid_spacing, pelc);

    }
    else if (autoset_wavefunction_grid) {

        // Processor grid specified but wavefunction grid was not
        SetupWavefunctionGrid(NPES, NX_GRID, NY_GRID, NZ_GRID, celldm, grid_spacing);

    }

    // Save results in Input map. Will clean up in the future when all code branches
    // are converted to C++
    ik = InputMap["processor_grid"];
    ProcessorGrid = ik->Vint;
    ProcessorGrid.vals[0] = pelc.pe_x;
    ProcessorGrid.vals[1] = pelc.pe_y;
    ProcessorGrid.vals[2] = pelc.pe_z;
    ik->Vint = ProcessorGrid;

    WavefunctionGrid.vals[0] = NX_GRID;
    WavefunctionGrid.vals[1] = NY_GRID;
    WavefunctionGrid.vals[2] = NZ_GRID;
    ik = InputMap["wavefunction_grid"];
    ik->Vint = WavefunctionGrid;
    ik = InputMap["grid_spacing"];
    ik->Readdoubleval = &grid_spacing;


    // Sanity check
    CheckAndTerminate(pelc.pe_x, 1, INT_MAX, "The value given for the processor grid X dimension is " + boost::lexical_cast<std::string>(pelc.pe_x) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_y, 1, INT_MAX, "The value given for the processor grid Y dimension is " + boost::lexical_cast<std::string>(pelc.pe_y) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_z, 1, INT_MAX, "The value given for the processor grid Z dimension is " + boost::lexical_cast<std::string>(pelc.pe_z) + " and only postive values are allowed.");


    // Set grid object up
    Rmg_G = new BaseGrid(NX_GRID, NY_GRID, NZ_GRID, pelc.pe_x, pelc.pe_y, pelc.pe_z, 0, lc.FG_RATIO);

    int FNX_GRID = NX_GRID * lc.FG_RATIO;
    int FNY_GRID = NY_GRID * lc.FG_RATIO;
    int FNZ_GRID = NZ_GRID * lc.FG_RATIO;



    // If the user has not specifically set the number of poisson multigrid levels use the max
    if(lc.poi_parm.levels == -1) {
        for(lc.poi_parm.levels = 6;lc.poi_parm.levels >= 0;lc.poi_parm.levels--) {
            bool poi_level_err = false;
            if ((FNX_GRID / (1 << (lc.poi_parm.levels+1))) < 2) poi_level_err = true;
            if ((FNY_GRID / (1 << (lc.poi_parm.levels+1))) < 2) poi_level_err = true;
            if ((FNZ_GRID / (1 << (lc.poi_parm.levels+1))) < 2) poi_level_err = true;
            if ((FNX_GRID % (1 << (lc.poi_parm.levels))) != 0) poi_level_err = true;
            if ((FNY_GRID % (1 << (lc.poi_parm.levels))) != 0) poi_level_err = true;
            if ((FNZ_GRID % (1 << (lc.poi_parm.levels))) != 0) poi_level_err = true;
            if (!poi_level_err) break;
        }
    }

    // Cutoff parameter -- do we want separate settings for each of these?
    lc.rhocparm = lc.cparm;
    ct.betacparm = ct.cparm * (ct.nxfgrid - 1);


    // If the user has not specifically set the number of kohn-sham multigrid levels use 2
    if(lc.eig_parm.levels == -1) lc.eig_parm.levels = 3;

    int checklevel;

    for(checklevel = 1;checklevel < lc.eig_parm.levels;checklevel++) {
        bool eig_level_err = false;
        if ((NX_GRID / (1 << checklevel)) < 3) eig_level_err = true;
        if ((NY_GRID / (1 << checklevel)) < 3) eig_level_err = true;
        if ((NZ_GRID / (1 << checklevel)) < 3) eig_level_err = true;
        if ((NX_GRID % (1 << checklevel)) != 0) eig_level_err = true;
        if ((NY_GRID % (1 << checklevel)) != 0) eig_level_err = true;
        if ((NZ_GRID % (1 << checklevel)) != 0) eig_level_err = true;
        if (eig_level_err) {
            lc.eig_parm.levels = checklevel - 1;
            if(lc.eig_parm.levels > 2) lc.eig_parm.levels = 2;
            if(pct.imgpe == 0) std::cout << "Too many eigenvalue multigrid levels specified. Resetting to " << lc.eig_parm.levels << std::endl;
            break;
        }
    }

}
