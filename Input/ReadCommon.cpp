

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

CONTROL *ReadCommon(int argc, char *argv[], char *cfile, CONTROL *pp)
{


    CONTROL lc;
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

    If.RegisterInputKey("potential_grid_refinement", &FG_RATIO, 0, 3, 2, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "Ratio of the fine grid to the coarse grid.", 
                     "potential_grid_refinement must be in the range (1 <= ratio <= 2). Resetting to the default value of 2.\n");

    If.RegisterInputKey("potential_acceleration_constant_step", &lc.potential_acceleration_constant_step, 0.0 - DBL_MIN, 2.0 + DBL_MIN, 1.111, 
                      CHECK_AND_FIX, OPTIONAL, 
                     "Time step used for constant potential acceleration.\n",
                     "potential_acceleration_constant_step must lie in the range (0.0, 2.0). Resetting to the default value of 0.0.\n");

    If.RegisterInputKey("potential_acceleration_poisson_step", &lc.potential_acceleration_poisson_step, 0.0 - DBL_MIN, 3.0 + DBL_MIN, 0.0, 
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
                     "",
                     "");

    If.RegisterInputKey("hartree_min_sweeps", &lc.hartree_min_sweeps, 5, 5 , 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("kohn_sham_pre_smoothing", &lc.eig_parm.gl_pre, 1, 5, 2,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("kohn_sham_post_smoothing", &lc.eig_parm.gl_pst, 1, 5, 2,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("kohn_sham_mucycles", &lc.eig_parm.mucycles, 1, 3, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("kohn_sham_fd_order", &lc.kohn_sham_fd_order, 4, 8, 6,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_pre_smoothing", &lc.poi_parm.gl_pre, 1, 6, 3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_post_smoothing", &lc.poi_parm.gl_pst, 1, 6, 3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_mucycles", &lc.poi_parm.mucycles, 1, 3, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_coarsest_steps", &lc.poi_parm.coarsest_steps, 10, 100, 25,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("kohn_sham_mg_levels", &lc.eig_parm.levels, 0, 2, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("poisson_mg_levels", &lc.poi_parm.levels, 0, 10, -1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

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

exit(0);

#if 0
    // Load them into the program options control structure
    po::options_description control("RMG control options");
    LoadInputKeys(&InputMap, &control);


    po::variables_map vm;
    std::stringstream ss;
    ss << sfile;

    // Parse the input file
    auto parsedOptions = po::parse_config_file(ss, control, true);
    po::store(parsedOptions, vm);
    po::notify(vm);
    auto unregistered = po::collect_unrecognized(parsedOptions.options, po::include_positional);

    // Process results

    control.add_options()
       //("description", po::value<std::string>(&ii)->default_value("RMG RUN"), "description of run")

        ("processor_grid", po::value<Ri::ReadVector<int> >(&ProcessorGrid), "Three-D (x,y,z) layout of the MPI processes.")
        ("coarse_grid", po::value<Ri::ReadVector<int> >(&CoarseGrid), "Three-D (x,y,z) dimensions of the global wavefunction grid.")
        ("threads_per_node", po::value<int>(&lc.THREADS_PER_NODE)->default_value(-1),"")

        ("charge_pulay_scale", po::value<double>(&lc.charge_pulay_scale)->default_value(0.50), "")
        ("charge_pulay_special_metrics_weight", po::value<double>(&lc.charge_pulay_special_metrics_weight)->default_value(100.0), "")
        ("charge_density_mixing", po::value<double>(&lc.mix)->default_value(0.5), "")
        ("rms_convergence_criterion", po::value<double>(&lc.thr_rms)->default_value(1.0E-7), "")
        ("relax_max_force", po::value<double>(&lc.thr_frc)->default_value(2.5E-3), "")
        ("occupation_electron_temperature_eV", po::value<double>(&lc.occ_width)->default_value(0.04), "")
        ("occupation_number_mixing", po::value<double>(&lc.occ_mix)->default_value(0.3), "")
        ("electric_field_magnitude", po::value<double>(&lc.e_field)->default_value(0.0), "")
        ("ionic_time_step", po::value<double>(&lc.iondt)->default_value(50), "")
        ("hartree_rms_ratio", po::value<double>(&lc.hartree_rms_ratio)->default_value(1000.0), "")
        ("a_length", po::value<double>(&celldm[0]), "")
        ("b_length", po::value<double>(&celldm[1]), "")
        ("c_length", po::value<double>(&celldm[2]), "")
        ("alpha", po::value<double>(&celldm[3]), "")
        ("beta", po::value<double>(&celldm[4]), "")
        ("gamma", po::value<double>(&celldm[5]), "")
        ("kohn_sham_time_step", po::value<double>(&lc.eig_parm.gl_step)->default_value(0.3), "")
        ("poisson_finest_time_step", po::value<double>(&lc.poi_parm.gl_step)->default_value(0.6), "")
        ("poisson_coarse_time_step", po::value<double>(&lc.poi_parm.sb_step)->default_value(0.6), "")
        ("energy_cutoff_parameter", po::value<double>(&lc.cparm)->default_value(1.75), "")
        ("potential_acceleration_constant_step", po::value<double>(&lc.potential_acceleration_constant_step)->default_value(0.0), "")
        ("potential_acceleration_poisson_step", po::value<double>(&lc.potential_acceleration_poisson_step)->default_value(0.0), "")
        ("projector_mixing", po::value<double>(&lc.prjmix)->default_value(0.5), "")
        ("Emin", po::value<double>(&lc.Emin)->default_value(-6.0), "")
        ("Emax", po::value<double>(&lc.Emax)->default_value(0.0), "")
        ("folded_spectrum_width", po::value<double>(&lc.folded_spectrum_width)->default_value(0.3), "")
        ("neb_spring_constant", po::value<double>(&lc.neb_spring_constant)->default_value(0.5), "")
        ("md_temperature", po::value<double>(&lc.nose.temp)->default_value(300), "")
        ("md_nose_oscillation_frequency_THz", po::value<double>(&lc.nose.fNose)->default_value(15.59), "")
        ("ionic_time_step", po::value<double>(&lc.iondt)->default_value(50), "")
        ("max_ionic_time_step", po::value<double>(&lc.iondt_max)->default_value(150), "")
        ("ionic_time_step_increase", po::value<double>(&lc.iondt_inc)->default_value(1.1), "")
//        ("ionic_time_step_decrease", po::value<double>(&lc.iondt_dec) DBL, s_ttt);
        ("system_charge", po::value<double>(&lc.background_charge)->default_value(0), "")
        ("energy_cutoff_parameter", po::value<double>(&lc.cparm)->default_value(1.75), "")

#if GPU_ENABLED
        ("gpu_direct_collectives", po::value<bool>(&lc.gpu_direct_collectives)->default_value(false), "")
#endif
        ("charge_pulay_special_metrics", po::value<bool>(&lc.charge_pulay_special_metrics)->default_value(false), "Flag to test whether or not the modified metrics should be used in Pulay mixing.")
//        ("write_pseudopotential_plots", po::value<bool>(NULL)->default_value(false), "")
        ("equal_initial_density", po::value<bool>(&lc.init_equal_density_flag)->default_value(false), "Specifies whether to set initial up and down density to be equal.")
        ("write_pdos", po::value<bool>(&lc.pdos_flag)->default_value(false), "Flag to write partial density of states.")
        ("mask_function_filtering", po::value<bool>(&lc.mask_function)->default_value(false), "")
        ("write_memory_report", po::value<bool>(&lc.write_memory_report)->default_value(false), "")
        ("sort_wavefunctions", po::value<bool>(&lc.sortflag)->default_value(false), "Sort wavefunctions by eigenvalue. Not needed if using subspace diagonalization.")
        ("initial_diagonalization", po::value<bool>(&lc.initdiag)->default_value(false), "Perform initial subspace diagonalization.")
        ("folded_spectrum", po::value<bool>(&lc.use_folded_spectrum)->default_value(false), "Use folded spectrum.")
//        ("relax_dynamic_timestep", po::value<bool>(NULL)->default_value(false), "")
        ("write_pdos", po::value<bool>(&lc.pdos_flag)->default_value(false), "")
// Depends on start flag
//        ("md_randomize_velocity", po::value<bool>(&lc.nose.randomvel)->default_value(false), "")
//        ("md_randomize_velocity", po::value<bool>(&lc.nose.randomvel)->default_value(true), "")
        ("scalapack_global_sums", po::value<bool>(&lc.scalapack_global_sums)->default_value(true), "")
    ;


    // Set grid info up
//    Rmg_G = new BaseGrid(NX_GRID, NY_GRID, NZ_GRID, pct.pe_x, pct.pe_y, pct.pe_z, 0, FG_RATIO);

    if (vm.count("description")) {
        //std::cout << "AAdescription=" << ii << std::endl;
    }

#endif


}
