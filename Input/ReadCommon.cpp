

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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
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


/**********************************************************************

   The RMG input file consists of a set of key-value pairs of the form.

   name=scalar    where scalar can be an integer, double or boolean value
   period_of_diagonalization=1
   charge_density_mixing = 0.5
   initial_diagonalization = true
   
   
   There are also strings and arrays which are delineated by double quotes so an integer array
   with three elements would be
   processor_grid = "2 2 2"

   while a string example would be
   description = "64 atom diamond cell test run at the gamma point"

   strings can span multiple lines so
   description = "
     64 atom diamond cell test run at the gamma point
     Uses Vanderbilt ultrasoft pseudopotential"

   would also be valid.

   Strong error checking is enforced on scalar types. This is done by
   registering the type before the input file is loaded and parsed. A class
   of type RmgInputFile If(cfile) is declared and the keys registered using
   the RegisterInputKey method. The method is overloaded and takes different
   arguments depending on what type of data is being read. For example scalar
   values are registered in the following way.

   RmgInputFile If(inputfile);
   If.RegisterInputKey("potential_grid_refinement", &FG_RATIO, 0, 3, 2,
                     CHECK_AND_FIX, OPTIONAL,
                     "Ratio of the fine grid to the coarse grid.",
                     "potential_grid_refinement must be in the range (1 <= ratio <= 2). Resetting to the default value of 2.\n");

   No matter what type of key is being registered the first argument of RegisterInputKey
   always consists of the key name as used in the input file. The second argument is
   always the destination where the data will be written. The last two arguments are a help
   message and an error message respectively.

   For integer and double scalar values the 3rd, 4th, and 5th values are the min, max and default
   values respectively. The sixth argument specifies the action to take if the value obtained
   from the input file does not fall within the acceptable range. CHECK_AND_FIX means set the 
   value to the default and continue execution while CHECK_AND_TERMINATE ends the program. The
   7th argument specifies whether the key is OPTIONAL or REQUIRED. If the key is REQUIRED the
   default value is ignored.

**********************************************************************/

namespace Ri = RmgInput;

void ReadCommon(int argc, char *argv[], char *cfile, CONTROL& lc)
{


//    CONTROL lc;
    PE_CONTROL pelc;
    RmgInputFile If(cfile);
 
    Ri::ReadVector<int> ProcessorGrid;
    Ri::ReadVector<int> CoarseGrid;
    int FG_RATIO;
    double celldm[6];
 
    If.RegisterInputKey("processor_grid", &ProcessorGrid, 3, REQUIRED, 
                     "Three-D (x,y,z) layout of the MPI processes.\n", 
                     "You must specify a triplet of (X,Y,Z) dimensions for the processor grid.\n");

    If.RegisterInputKey("coarse_grid", &CoarseGrid, 3, REQUIRED, 
                     "Three-D (x,y,z) layout of the MPI processes.\n", 
                     "You must specify a triplet of (X,Y,Z) dimensions for the coarse grid.\n");

    // Deault of zero is OK because this means to try to set it automatically later on.
    // The value of 64 covers any possible hardware scenario I can imagine currently but might
    // need to be adjusted at some point in the future.
    If.RegisterInputKey("threads_per_node", &lc.THREADS_PER_NODE, 0, 64, 0, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "Number of threads each MPI process will use.\n", 
                     "threads_per_node cannnot be a negative number and must be less than 64.\n");

    If.RegisterInputKey("potential_grid_refinement", &FG_RATIO, 0, 3, 2, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "Ratio of the fine grid to the coarse grid.", 
                     "potential_grid_refinement must be in the range (1 <= ratio <= 2). Resetting to the default value of 2.\n");

    If.RegisterInputKey("potential_acceleration_constant_step", &lc.potential_acceleration_constant_step, 0.0, 2.0, 1.111, 
                      CHECK_AND_FIX, OPTIONAL, 
                     "Time step used for constant potential acceleration.\n",
                     "potential_acceleration_constant_step must lie in the range (0.0, 2.0). Resetting to the default value of 0.0.\n");

    If.RegisterInputKey("potential_acceleration_poisson_step", &lc.potential_acceleration_poisson_step, 0.0, 3.0, 0.0, 
                      CHECK_AND_FIX, OPTIONAL, 
                     "Time step used for poisson potential acceleration.\n",
                     "potential_acceleration_poisson_step must lie in the range (0.0, 3.0). Resetting to the default value of 0.0.\n");

    If.RegisterInputKey("ionic_time_step", &lc.iondt, 0.0, DBL_MAX, 50.0, 
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "Ionic time step for use in molecular dynamics and structure optimizations.\n",
                     "ionic_time_step must be greater than 0.0.\n");

    If.RegisterInputKey("unoccupied_states_per_kpoint", &lc.num_unocc_states, 0, INT_MAX, 10, 
                     CHECK_AND_TERMINATE, OPTIONAL, 
                     "The number of unoccupied orbitals.\n", 
                     "unoccupied_states_per_kpoint must be greater than 0. Terminating.\n");

    If.RegisterInputKey("period_of_diagonalization", &lc.diag, 0, INT_MAX, 1, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "Diagonalization period (per scf step).\n", 
                     "Diagonalization period must be greater than 0. Resetting to the default value of 1.\n");

    If.RegisterInputKey("end_diagonalization_step", &lc.end_diag, 0, INT_MAX, 1000000,
                     CHECK_AND_FIX, OPTIONAL, 
                     "Stop diagonalizing after end_diag steps.\n", 
                     "end_diag must be greater than 0. Resetting to the default value of 1000000\n");

    If.RegisterInputKey("max_scf_steps", &lc.max_scf_steps, 0, INT_MAX, 500,
                     CHECK_AND_FIX, OPTIONAL, 
                     "Maximum number of self consistent steps to perform.\n", 
                     "max_scf_steps must be greater than 0. Resetting to the default value of 500\n");

    If.RegisterInputKey("charge_pulay_order", &lc.charge_pulay_order, 5, 5, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_pulay_refresh", &lc.charge_pulay_refresh, 0, 0, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("write_data_period", &lc.checkpoint, 5, 5, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("write_eigvals_period", &lc.write_eigvals_period, 1, 10, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("max_md_steps", &lc.max_md_steps, 100, 100, 100,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("hartree_max_sweeps", &lc.hartree_max_sweeps, 5, 100, 30,
                     CHECK_AND_FIX, OPTIONAL,
                     "Maximum number of hartree iterations to perform per scf step.\n",
                     "hartree_max_sweeps must lie in the range (5,100). Resetting to the default value of 30.\n");

    If.RegisterInputKey("hartree_min_sweeps", &lc.hartree_min_sweeps, 0, 5, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "Minimum number of hartree iterations to perform per scf step.\n",
                     "hartree_min_sweeps must lie in the range (0.5). Resetting to the default value of 5.\n");

    If.RegisterInputKey("kohn_sham_pre_smoothing", &lc.eig_parm.gl_pre, 1, 5, 2,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of global grid pre-smoothing steps to perform before a multigrid preconditioner iteration.\n",
                     "kohn_sham_pre_smoothing must lie in the range (1,5). Resetting to the default value of 2.\n");

    If.RegisterInputKey("kohn_sham_post_smoothing", &lc.eig_parm.gl_pst, 1, 5, 2,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of global grid post-smoothing steps to perform after a multigrid preconditioner iteration.\n",
                     "kohn_sham_post_smoothing must lie in the range (1,5). Resetting to the default value of 2.\n");

    If.RegisterInputKey("kohn_sham_mucycles", &lc.eig_parm.mucycles, 1, 3, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of mu (also known as W) cycles to use in the kohn-sham multigrid preconditioner.\n",
                     "kohn_sham_mucycles must lie in the range (1,3). Resetting to the default value of 1.\n");

    If.RegisterInputKey("kohn_sham_fd_order", &lc.kohn_sham_fd_order, 4, 8, 6,
                     CHECK_AND_FIX, OPTIONAL,
                     "Order of the global grid finite difference operators to be used in the kohn-sham multigrid preconditioner.\n ",
                     "kohn_sham_fd_order must lie in the range (4,8). Resetting to the default value of 6.\n");

    If.RegisterInputKey("poisson_pre_smoothing", &lc.poi_parm.gl_pre, 1, 6, 3,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of global hartree grid pre-smoothing steps to perform before a multigrid iteration.\n",
                     "poisson_pre_smoothing must lie in the range (1,6). Resetting to the default value of 3.\n");

    If.RegisterInputKey("poisson_post_smoothing", &lc.poi_parm.gl_pst, 1, 6, 3,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of global hartree grid post-smoothing steps to perform after a multigrid iteration.\n",
                     "");

    If.RegisterInputKey("poisson_mucycles", &lc.poi_parm.mucycles, 1, 3, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of mu (also known as W) cycles to use in the hartree multigrid solver.\n",
                     "poisson_mucycles must lie in the range (1,3). Resetting to the default value of 1.\n");

    If.RegisterInputKey("poisson_coarsest_steps", &lc.poi_parm.coarsest_steps, 10, 100, 25,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of smoothing steps to use on the coarsest level in the hartree multigrid solver.\n",
                     "poisson_coarsest_steps must lie in the range (10,100). Resetting to the default value of 25.\n");

    If.RegisterInputKey("kohn_sham_mg_levels", &lc.eig_parm.levels, 0, 2, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of multigrid levels to use in the kohn-sham multigrid preconditioner.\n",
                     "kohn_sham_mg_levels must lie in the range (0,2). Resetting to the default value of 1.\n");

    If.RegisterInputKey("poisson_mg_levels", &lc.poi_parm.levels, -1, 6, -1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of multigrid levels to use in the hartree multigrid solver.\n",
                     "poisson_mg_levels must lie in the range (-1,6) where -1=automatic. Resetting to the default value of automatic (-1).\n");

    If.RegisterInputKey("fine_grid_non_local_pp", &lc.nxfgrid, 1, 4, 4,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("scalapack_block_factor", &lc.scalapack_block_factor, 32, 512,64,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("E_POINTS", &lc.E_POINTS, 201, 201, 201,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("max_rmg_steps", &lc.max_rmg_steps, 1, 1, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("md_number_of_nose_thermostats", &lc.nose.m, 5, 5, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("dynamic_time_delay", &lc.relax_steps_delay, 5, 5, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("dynamic_time_counter", &lc.relax_steps_counter, 0, 0 , 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("scf_steps_offset", &lc.scf_steps, 0, 0, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("total_scf_steps_offset", &lc.total_scf_steps, 0, 0, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("md_steps_offset", &lc.md_steps, 0, 0, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("b_spline_order", &lc.interp_order, 5, 5, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("b_spline_trade_order", &lc.interp_trade, 3, 3, 3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_density_mixing", &lc.mix, 0.0, 1.0, 0.5,
                     CHECK_AND_FIX, OPTIONAL,
                     "Proportion of the current charge density to replace with the new density after each scf step.\n",
                     "charge_density_mixing must lie in the range (0.0, 1.0)\n");


    // Booleans next. Booleans are never required.
    If.RegisterInputKey("charge_pulay_special_metrics", &lc.charge_pulay_special_metrics, false,
                        "Flag to test whether or not the modified metrics should be used in Pulay mixing.");

//    If.RegisterInputKey("write_pseudopotential_plots", NULL, false,
//                        "");

    If.RegisterInputKey("equal_initial_density", &lc.init_equal_density_flag, false,
                        "Specifies whether to set initial up and down density to be equal.");

    If.RegisterInputKey("write_pdos", &lc.pdos_flag, false,
                        "Flag to write partial density of states.");

    If.RegisterInputKey("mask_function_filtering", &lc.mask_function, false,
                        "");

    If.RegisterInputKey("sort_wavefunctions", &lc.sortflag, false, 
                        "Sort wavefunctions by eigenvalue. Not needed if using subspace diagonalization.");

    If.RegisterInputKey("initial_diagonalization", &lc.initdiag, false, 
                        "Perform initial subspace diagonalization.");

    If.RegisterInputKey("folded_spectrum", &lc.use_folded_spectrum, false, 
                         "Use folded spectrum.");

//    If.RegisterInputKey("relax_dynamic_timestep", NULL, false,
//                        "");

    If.RegisterInputKey("md_randomize_velocity", &lc.nose.randomvel, true,
                        "");

    If.RegisterInputKey("scalapack_global_sums", &lc.scalapack_global_sums, true,
                        "");

#if 0
    If.RegisterInputKey("charge_pulay_scale", &lc.charge_pulay_scale, min, max, 0.50,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_pulay_special_metrics_weight", &lc.charge_pulay_special_metrics_weight, min, max, 100.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("rms_convergence_criterion", &lc.thr_rms, min, max, 1.0E-7,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("relax_max_force", &lc.thr_frc, min, max, 2.5E-3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("occupation_electron_temperature_eV", &lc.occ_width, min, max, 0.04,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("occupation_number_mixing", &lc.occ_mix, min, max, 0.3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("electric_field_magnitude", &lc.e_field, min, max, 0.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("hartree_rms_ratio", &lc.hartree_rms_ratio, min, max, 1000.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("kohn_sham_time_step", &lc.eig_parm.gl_step, min, max, 0.3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_finest_time_step", &lc.poi_parm.gl_step, min, max, 0.6,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_coarse_time_step", &lc.poi_parm.sb_step, min, max, 0.6,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("energy_cutoff_parameter", &lc.cparm, min, max, 1.75,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("projector_mixing", &lc.prjmix, min, max, 0.5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("Emin", &lc.Emin, min, max, -6.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("Emax", &lc.Emax, min, max, 0.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("folded_spectrum_width", &lc.folded_spectrum_width, min, max, 0.3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("neb_spring_constant", &lc.neb_spring_constant, min, max, 0.5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("md_temperature", &lc.nose.temp, min, max, 300,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("md_nose_oscillation_frequency_THz", &lc.nose.fNose, min, max, 15.59,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("max_ionic_time_step", &lc.iondt_max, min, max, 150,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("ionic_time_step_increase", &lc.iondt_inc, min, max, 1.1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("system_charge", &lc.background_charge, min, max, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

#endif

    If.LoadInputKeys();


    // Check items that require custom handling
    try {
        pelc.pe_x = ProcessorGrid.vals.at(0);
        pelc.pe_y = ProcessorGrid.vals.at(1);
        pelc.pe_z = ProcessorGrid.vals.at(2);
    }
    catch (const std::out_of_range& oor) {
        throw RmgFatalException() << "You must specify a triplet of (X,Y,Z) dimensions for the processor grid.\n";
    }
    CheckAndTerminate(pelc.pe_x, 1, INT_MAX, "The value given for the processor grid X dimension is " + std::to_string(pelc.pe_x) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_y, 1, INT_MAX, "The value given for the processor grid Y dimension is " + std::to_string(pelc.pe_y) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_z, 1, INT_MAX, "The value given for the processor grid Z dimension is " + std::to_string(pelc.pe_z) + " and only postive values are allowed.");

    int NX_GRID=1, NY_GRID=1, NZ_GRID=1;
    try {
        NX_GRID = CoarseGrid.vals.at(0);
        NY_GRID = CoarseGrid.vals.at(1);
        NZ_GRID = CoarseGrid.vals.at(2);
    }
    catch (const std::out_of_range& oor) {
        throw RmgFatalException() << "You must specify a triplet of (X,Y,Z) dimensions for the coarse grid.\n";
    }
    CheckAndTerminate(pelc.pe_x, 1, INT_MAX, "The value given for the global wavefunction grid X dimension is " + std::to_string(NX_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_y, 1, INT_MAX, "The value given for the global wavefunction grid Y dimension is " + std::to_string(NY_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_z, 1, INT_MAX, "The value given for the global wavefunction grid Z dimension is " + std::to_string(NZ_GRID) + " and only postive values are allowed.");

    int FNX_GRID = NX_GRID * FG_RATIO;
    int FNY_GRID = NY_GRID * FG_RATIO;
    int FNZ_GRID = NZ_GRID * FG_RATIO;

    // If the user has not specifically set the number of poisson multigrid levels use the max
    if(lc.poi_parm.levels == -1) {
        for(ct.poi_parm.levels = 6;ct.poi_parm.levels >= 0;ct.poi_parm.levels--) {
            bool poi_level_err = false;
            if ((FNX_GRID / (1 << ct.poi_parm.levels)) < 3) poi_level_err = true;
            if ((FNY_GRID / (1 << ct.poi_parm.levels)) < 3) poi_level_err = true;
            if ((FNZ_GRID / (1 << ct.poi_parm.levels)) < 3) poi_level_err = true;
            if ((FNX_GRID % (1 << ct.poi_parm.levels)) != 0) poi_level_err = true;
            if ((FNY_GRID % (1 << ct.poi_parm.levels)) != 0) poi_level_err = true;
            if ((FNZ_GRID % (1 << ct.poi_parm.levels)) != 0) poi_level_err = true;
            if (!poi_level_err) break;
        }
    }


exit(0);

#if 0

    // Parse the input file
    auto parsedOptions = po::parse_config_file(ss, control, true);
    po::store(parsedOptions, vm);
    po::notify(vm);
    auto unregistered = po::collect_unrecognized(parsedOptions.options, po::include_positional);

    // Process results

    control.add_options()
       //("description", po::value<std::string>(&ii)->default_value("RMG RUN"), "description of run")

        ("a_length", po::value<double>(&celldm[0]), "")
        ("b_length", po::value<double>(&celldm[1]), "")
        ("c_length", po::value<double>(&celldm[2]), "")
        ("alpha", po::value<double>(&celldm[3]), "")
        ("beta", po::value<double>(&celldm[4]), "")
        ("gamma", po::value<double>(&celldm[5]), "")
//        ("ionic_time_step_decrease", po::value<double>(&lc.iondt_dec) DBL, s_ttt);

    ;


    // Set grid info up
//    Rmg_G = new BaseGrid(NX_GRID, NY_GRID, NZ_GRID, pct.pe_x, pct.pe_y, pct.pe_z, 0, FG_RATIO);

    if (vm.count("description")) {
        //std::cout << "AAdescription=" << ii << std::endl;
    }

#endif


}
