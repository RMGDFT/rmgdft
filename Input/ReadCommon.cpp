

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




CONTROL *ReadCommon(int argc, char *argv[], char *cfile, CONTROL *pp)
{


    CONTROL lc;
    PE_CONTROL pelc;
    RmgInputFile If(cfile);
 
//    Ri::ReadVector<int> ProcessorGrid;
//    Ri::ReadVector<int> CoarseGrid;
    int FG_RATIO;
    double celldm[6];
 
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

    If.LoadInputKeys();

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
        ("ionic_time_step", po::value<double>(&lc.iondt)->default_value(50.0), "Ionic time step for use in molecular dynamics and structure optimizations.")
        ("unoccupied_states_per_kpoint", po::value<int>(&lc.num_unocc_states)->default_value(10), "The number of unoccupied orbitals.")
        ("period_of_diagonalization", po::value<int>(&lc.diag)->default_value(1), "Diagonalization period (per scf step).")
        ("end_diagonalization_step", po::value<int>(&lc.end_diag)->default_value(9999999), "Stop diagonalizing after end_diag steps.")
        ("threads_per_node", po::value<int>(&lc.THREADS_PER_NODE)->default_value(-1),"")
        ("charge_pulay_order", po::value<int>(&lc.charge_pulay_order)->default_value(5),"")
        ("charge_pulay_refresh", po::value<int>(&lc.charge_pulay_refresh)->default_value(0),"")
        ("max_scf_steps", po::value<int>(&lc.max_scf_steps)->default_value(500),"")
        ("write_data_period", po::value<int>(&lc.checkpoint)->default_value(5),"")
        ("write_eigvals_period", po::value<int>(&lc.write_eigvals_period)->default_value(5),"")
        ("max_md_steps", po::value<int>(&lc.max_md_steps)->default_value(100),"")
        ("hartree_max_sweeps", po::value<int>(&lc.hartree_max_sweeps)->default_value(100),"")
        ("hartree_min_sweeps", po::value<int>(&lc.hartree_min_sweeps)->default_value(5),"")
        ("kohn_sham_pre_smoothing", po::value<int>(&lc.eig_parm.gl_pre)->default_value(2),"")
        ("kohn_sham_post_smoothing", po::value<int>(&lc.eig_parm.gl_pst)->default_value(1),"")
        ("kohn_sham_mucycles", po::value<int>(&lc.eig_parm.mucycles)->default_value(1),"")
        ("kohn_sham_fd_order", po::value<int>(&lc.kohn_sham_fd_order)->default_value(6),"")
        ("poisson_pre_smoothing", po::value<int>(&lc.poi_parm.gl_pre)->default_value(3),"")
        ("poisson_post_smoothing", po::value<int>(&lc.poi_parm.gl_pst)->default_value(3),"")
        ("poisson_mucycles", po::value<int>(&lc.poi_parm.mucycles)->default_value(1),"")
        ("poisson_coarsest_steps", po::value<int>(&lc.poi_parm.coarsest_steps)->default_value(80),"")
        ("kohn_sham_mg_levels", po::value<int>(&lc.eig_parm.levels)->default_value(1),"")
        ("poisson_mg_levels", po::value<int>(&lc.poi_parm.levels)->default_value(-1),"")
        ("fine_grid_non_local_pp", po::value<int>(&lc.nxfgrid)->default_value(4),"")
        ("scalapack_block_factor", po::value<int>(&lc.scalapack_block_factor)->default_value(32),"")
        ("E_POINTS", po::value<int>(&lc.E_POINTS)->default_value(201),"")
        ("max_rmg_steps", po::value<int>(&lc.max_rmg_steps)->default_value(1),"")
        ("md_number_of_nose_thermostats", po::value<int>(&lc.nose.m)->default_value(5),"")
        ("dynamic_time_delay", po::value<int>(&lc.relax_steps_delay)->default_value(5),"")
        ("dynamic_time_counter", po::value<int>(&lc.relax_steps_counter)->default_value(0),"")
        ("scf_steps_offset", po::value<int>(&lc.scf_steps)->default_value(0),"")
        ("total_scf_steps_offset", po::value<int>(&lc.total_scf_steps)->default_value(0),"")
        ("md_steps_offset", po::value<int>(&lc.md_steps)->default_value(0),"")
        ("b_spline_order", po::value<int>(&lc.interp_order)->default_value(5),"")
        ("b_spline_trade_order", po::value<int>(&lc.interp_trade)->default_value(3),"")

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


    
    try {
        pelc.pe_x = ProcessorGrid.vals.at(0);
        pelc.pe_y = ProcessorGrid.vals.at(1);
        pelc.pe_z = ProcessorGrid.vals.at(2);
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "You must specify a triplet of (X,Y,Z) dimensions for the processor grid." << std::endl;
        rmg_error_handler(__FILE__, __LINE__, "Terminating.\n");
    }
    CheckValueAndTerminate(pelc.pe_x, 1, INT_MAX, "The value given for the processor grid X dimension is " + std::to_string(pelc.pe_x) + " and only postive values are allowed.");
    CheckValueAndTerminate(pelc.pe_y, 1, INT_MAX, "The value given for the processor grid Y dimension is " + std::to_string(pelc.pe_y) + " and only postive values are allowed.");
    CheckValueAndTerminate(pelc.pe_z, 1, INT_MAX, "The value given for the processor grid Z dimension is " + std::to_string(pelc.pe_z) + " and only postive values are allowed.");


    int NX_GRID=1, NY_GRID=1, NZ_GRID=1;
    try {
        NX_GRID = CoarseGrid.vals.at(0);
        NY_GRID = CoarseGrid.vals.at(1);
        NZ_GRID = CoarseGrid.vals.at(2);
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "You must specify a triplet of (X,Y,Z) dimensions for the coarse grid." << std::endl;
        rmg_error_handler(__FILE__, __LINE__, "Terminating.\n");
    }
    CheckValueAndTerminate(pelc.pe_x, 1, INT_MAX, "The value given for the global wavefunction grid X dimension is " + std::to_string(NX_GRID) + " and only postive values are allowed.");
    CheckValueAndTerminate(pelc.pe_y, 1, INT_MAX, "The value given for the global wavefunction grid Y dimension is " + std::to_string(NY_GRID) + " and only postive values are allowed.");
    CheckValueAndTerminate(pelc.pe_z, 1, INT_MAX, "The value given for the global wavefunction grid Z dimension is " + std::to_string(NZ_GRID) + " and only postive values are allowed.");

    // Set grid info up
//    Rmg_G = new BaseGrid(NX_GRID, NY_GRID, NZ_GRID, pct.pe_x, pct.pe_y, pct.pe_z, 0, FG_RATIO);

    if (vm.count("potential_grid_refinement"))
        FG_RATIO = CheckValueAndFix(FG_RATIO, 0, 3, 2, "potential_grid_refinement must be in the range (1 <= ratio <= 2). Resetting to the default value of 2.\n");

    if (vm.count("potential_acceleration_constant_step"))
        lc.potential_acceleration_constant_step = 
            CheckValueAndFix(lc.potential_acceleration_constant_step, 0.0 - DBL_MIN, 2.0 + DBL_MIN, 0.0, 
            "potential_acceleration_constant_step must lie in the range (0.0, 2.0). Resetting to the default value of 0.0\n");

    if (vm.count("potential_acceleration_poisson_step"))
        lc.potential_acceleration_constant_step = 
            CheckValueAndFix(lc.potential_acceleration_poisson_step, 0.0 - DBL_MIN, 3.0 + DBL_MIN, 0.0, 
            "potential_acceleration_poisson_step must lie in the range (0.0, 3.0). Resetting to the default value of 0.0\n");

    if (vm.count("ionic_time_step")) 
        ct.iondt = CheckValueAndFix(ct.iondt, 0.0 - DBL_MIN, DBL_MAX, 50.0, "ionic_time_step must be >= 0.0. Resetting to the default value of 50.0.\n");

    if (vm.count("unoccupied_states_per_kpoint"))
        ct.num_unocc_states = CheckValueAndFix(ct.num_unocc_states, 0, INT_MAX, 10, "unoccupied_states_per_kpoint must be > 0. Resetting to the default value of 10\n");

    if (vm.count("description")) {
        //std::cout << "AAdescription=" << ii << std::endl;
    }

#endif


}
