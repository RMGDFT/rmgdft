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
   If.RegisterInputKey("potential_grid_refinement", &lc.FG_RATIO, 0, 3, 2,
                     CHECK_AND_FIX, OPTIONAL,
                     "Ratio of the fine grid to the wavefunction grid.",
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

void ReadCommon(int argc, char *argv[], char *cfile, CONTROL& lc, PE_CONTROL& pelc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    RmgInputFile If(cfile, InputMap, pct.img_comm);
    std::string CalculationMode;
    std::string DiscretizationType;
    std::string EnergyOutputType;
    std::string Description;
    std::string Infile;
    std::string Outfile;
    std::string Infile_tddft;
    std::string Outfile_tddft;
    std::string wfng_file;
    std::string rhog_file;
    std::string vxc_file;
 
    static Ri::ReadVector<int> ProcessorGrid;
    Ri::ReadVector<int> DefProcessorGrid({{1,1,1}});
    static Ri::ReadVector<int> WavefunctionGrid;
    Ri::ReadVector<int> DefWavefunctionGrid({{1,1,1}});
    Ri::ReadVector<int> kpoint_mesh;
    Ri::ReadVector<int> def_kpoint_mesh({{1,1,1}});
    Ri::ReadVector<int> kpoint_is_shift;
    Ri::ReadVector<int> def_kpoint_is_shift({{0,0,0}});

    static Ri::ReadVector<int> DipoleCorrection;
    Ri::ReadVector<int> DefDipoleCorrection({{0,0,0}});

    double celldm[6];
    static double grid_spacing;
    double a0[3], a1[3], a2[3], omega;

 
    If.RegisterInputKey("description", &Description, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Description of the run.\n", 
                     "");

    If.RegisterInputKey("pseudopotential", NULL , "",
                     CHECK_AND_FIX, OPTIONAL,
                     "External pseudopotentials.\n", 
                     "");

    If.RegisterInputKey("input_wave_function_file", &Infile, "Waves/wave.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Input file/path to  read wavefunctions and other binary data from on a restart.\n", 
                     "");
    If.RegisterInputKey("input_tddft_file", &Infile_tddft, "Waves/wave_tddft.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Input file/path to  read wavefunctions and other binary data from on a restart.\n", 
                     "");

    If.RegisterInputKey("output_wave_function_file", &Outfile, "Waves/wave.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Output file/path to store wavefunctions and other binary data.\n", 
                     "");
    If.RegisterInputKey("output_tddft_file", &Outfile_tddft, "Waves/wave_tddft.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Output file/path to store wavefunctions and other binary data.\n", 
                     "");
    If.RegisterInputKey("restart_tddft", &lc.restart_tddft, false, 
                        "restart TDDFT");

    If.RegisterInputKey("processor_grid", &ProcessorGrid, &DefProcessorGrid, 3, OPTIONAL, 
                     "Three-D (x,y,z) layout of the MPI processes.\n", 
                     "You must specify a triplet of (X,Y,Z) dimensions for the processor grid.\n");

    If.RegisterInputKey("wavefunction_grid", &WavefunctionGrid, &DefWavefunctionGrid, 3, OPTIONAL, 
                     "Three-D (x,y,z) dimensions of the grid the wavefunctions are defined on.\n", 
                     "You must specify a triplet of (X,Y,Z) dimensions for the wavefunction grid.\n");

    If.RegisterInputKey("dipole_correction", &DipoleCorrection, &DefDipoleCorrection, 3, OPTIONAL, 
                     "(1,1,1) for molecule, dipole correction in all directions. \n", 
                     "(0,0,0) means no correction by default, (1,0,0) or others have not programed\n");

    If.RegisterInputKey("kpoint_mesh", &kpoint_mesh, &def_kpoint_mesh, 3, OPTIONAL, 
                     "Three-D layout of the kpoint mesh.\n", 
                     "You must specify a triplet of coordinate dimensions for the kpoint_mesh.\n");

    If.RegisterInputKey("kpoint_is_shift", &kpoint_is_shift, &def_kpoint_is_shift, 3, OPTIONAL, 
                     "Three-D layout of the kpoint shift.\n", 
                     "You must specify a triplet of coordinate dimensions for kpoint_is_shift.\n");

    int ibrav;
    If.RegisterInputKey("bravais_lattice_type", NULL, &ibrav, "Orthorhombic Primitive",
                     CHECK_AND_TERMINATE, OPTIONAL, bravais_lattice_type,
                     "Bravais Lattice Type.\n", 
                     "bravais_lattice_type not found.\n");

    If.RegisterInputKey("start_mode", NULL, &lc.runflag, "LCAO Start",
                     CHECK_AND_TERMINATE, OPTIONAL, start_mode,
                     "Type of run. Choices are \"Random Start\", \"Restart From File\", or \"LCAO Start\".\n", 
                     "start_mode must be one of  \"Random Start\", \"Restart From File\", or \"LCAO Start\". Terminating.\n");

    If.RegisterInputKey("subdiag_driver", NULL, &lc.subdiag_driver, "auto",
                     CHECK_AND_FIX, OPTIONAL, subdiag_driver,
                     "Driver type used for subspace diagonalization of the eigenvectors.\n", 
                     "subdiag_driver must be lapack, scalapack, magma or auto. Resetting to auto.\n");

    If.RegisterInputKey("kohn_sham_solver", NULL, &lc.kohn_sham_solver, "multigrid",
                     CHECK_AND_FIX, OPTIONAL, kohn_sham_solver,
                     "Kohn-Sham solver.\n", 
                     "kohn_sham_solver must be multigrid or davidson. Resetting to multigrid.\n");

    If.RegisterInputKey("poisson_solver", NULL, &lc.poisson_solver, "pfft",
                     CHECK_AND_FIX, OPTIONAL, poisson_solver,
                     "poisson solver.\n", 
                     "poisson_solver must be multigrid or pfft. Resetting to pfft.\n");


    If.RegisterInputKey("crds_units", NULL, NULL, "Bohr",
                     CHECK_AND_FIX, OPTIONAL, crds_units,
                     "Units for the atomic coordinates.\n", 
                     "Coordinates must be specified in either Bohr or Angstrom.\n");

    If.RegisterInputKey("charge_mixing_type", NULL, &lc.charge_mixing_type, "Pulay",
                     CHECK_AND_TERMINATE, OPTIONAL, charge_mixing_type,
                     "Type of charge density mixing to use. Linear and Pulay are the available options.\n", 
                     "charge_mixing_type must be either \"Linear\" or \"Pulay\". Terminating.\n");
    
    If.RegisterInputKey("charge_analysis", NULL, &lc.charge_analysis_type, "Voronoi",
                     CHECK_AND_TERMINATE, OPTIONAL, charge_analysis,
                     "Type of charge nalysis to use. Only Voronoi deformation density is currently available.\n", 
                     "charge_analysis must be either \"Voronoi\" or \"None\". Terminating.\n");
    
    If.RegisterInputKey("charge_analysis_period", &lc.charge_analysis_period, 0, 500, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "How often to  perform and write out charge analysis.",
                     "charge_analysis_write_period must lie in the range (1,500). Resetting to the default value of 0.\n");
    
    If.RegisterInputKey("dipole_moment", &lc.dipole_moment, false, 
                        "Turns on calculation of dipole moment for the entire cell.");

    If.RegisterInputKey("vdwdf_grid_type", NULL, NULL, "Coarse",
                     CHECK_AND_TERMINATE, OPTIONAL, vdwdf_grid_type,
                     "Type of grid to use when computing vdw-df correlation.\n", 
                     "vdwdf_grid_type be either \"Coarse\" or \"Fine\". Terminating.\n");

    If.RegisterInputKey("relax_mass", NULL, &lc.relax_mass, "Atomic",
                     CHECK_AND_TERMINATE, OPTIONAL, relax_mass,
                     "Mass to use for structural relaxation, either atomic masses, or use the mass of carbon for all atoms.\n", 
                     "relax_mass must be either \"Atomic\" or \"Equal\". Terminating.\n");

    If.RegisterInputKey("md_integration_order", NULL, &lc.mdorder, "5th Beeman-Velocity Verlet",
                     CHECK_AND_TERMINATE, OPTIONAL, md_integration_order,
                     "Integration order for molecular dynamics.\n", 
                     "md_integration_order must be either \"2nd Velocity Verlet\", \"3rd Beeman-Velocity Verlet\" or \"5th Beeman-Velocity Verlet\". Terminating.\n");

    If.RegisterInputKey("z_average_output_mode", NULL, &lc.zaverage, "None",
                     CHECK_AND_TERMINATE, OPTIONAL, z_average_output_mode,
                     "z_average_output_mode.\n", 
                     "z_average_output_mode not supported. Terminating.\n");

    If.RegisterInputKey("atomic_coordinate_type", NULL, &lc.crd_flag, "Absolute",
                     CHECK_AND_TERMINATE, OPTIONAL, atomic_coordinate_type,
                     "Flag indicated whether or not atomic coordinates are absolute or cell relative.\n", 
                     "atomic_coordinate_type must be either \"Absolute\" or \"Cell Relative\". Terminating.\n");

    If.RegisterInputKey("calculation_mode", NULL, &lc.forceflag, "Quench Electrons",
                     CHECK_AND_TERMINATE, OPTIONAL, calculation_mode,
                     "Type of calculation to perform.\n", 
                     "calculation_mode not available.\n");

    If.RegisterInputKey("relax_method", NULL, &lc.relax_method, "Fast Relax",
                     CHECK_AND_TERMINATE, OPTIONAL, relax_method,
                     "Type of relaxation method to use for structural optimizations.\n", 
                     "relax_method not supported.\n");

    If.RegisterInputKey("md_temperature_control", NULL, &lc.tcontrol, "Nose Hoover Chains",
                     CHECK_AND_TERMINATE, OPTIONAL, md_temperature_control,
                     "Type of temperature control method to use in molecular dynamics.\n", 
                     "md_temperature_control type not supported.\n");

    If.RegisterInputKey("md_temperature", &lc.nose.temp, 0.0, DBL_MAX, 300.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "Target MD Temperature.\n",
                     "md_temperature must be a positive number.\n");

    If.RegisterInputKey("md_nose_oscillation_frequency_THz", &lc.nose.fNose, 0.0, DBL_MAX, 15.59,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "md_nose_oscillation_frequency_THz must be a positive real number.");

    If.RegisterInputKey("discretization_type", &DiscretizationType, &lc.discretization, "Central",
                     CHECK_AND_FIX, OPTIONAL, discretization_type,
                     "Type of discretization to use for the Kohn-Sham equations. Mehrstellen or Central types are implemented.\n", 
                     "discretization_type must be either \"Mehrstellen\" or \"Central\". Setting to \"Central\".\n");

    If.RegisterInputKey("energy_output_units", &EnergyOutputType, &lc.energy_output_units, "Hartrees",
                     CHECK_AND_FIX, OPTIONAL, energy_output_units,
                     "Units to be used when writing energy values to the output file. Hartrees or Rydbergs are available.\n", 
                     "energy_output_units must be either \"Hartrees\" or \"Rydbergs\". Setting to \"Hartrees\".\n");

    If.RegisterInputKey("boundary_condition_type", NULL, &lc.boundaryflag, "Periodic",
                     CHECK_AND_TERMINATE, OPTIONAL, boundary_condition_type,
                     "Boundary condition type Only periodic is currently implemented.\n", 
                     "discretization_type must be Periodic.\n");

    If.RegisterInputKey("exchange_correlation_type", NULL, &lc.xctype, "AUTO_XC",
                     CHECK_AND_TERMINATE, OPTIONAL, exchange_correlation_type,
                     "Type of functional for exchange-correlation.\n", 
                     "exchange_correlation_type not supported. Terminating.\n");

    If.RegisterInputKey("occupations_type", NULL, &lc.occ_flag, "Fermi Dirac",
                     CHECK_AND_TERMINATE, OPTIONAL, occupations_type,
                     "Method used to set the occupations of the electronic orbitals.\n", 
                     "occupations_type not supported. Terminating.\n");

    If.RegisterInputKey("interpolation_type", NULL, &lc.interp_flag, "FFT",
                     CHECK_AND_TERMINATE, OPTIONAL, interpolation_type,
                     "Interpolation method for transferring data between the potential grid and the wavefunction grid.\n", 
                     "interpolation_type not supported. Terminating.\n");

    If.RegisterInputKey("a_length", &celldm[0], 0.0, DBL_MAX, 0.0, 
                     CHECK_AND_TERMINATE, REQUIRED, 
                     "First lattice constant.\n", 
                     "a_length must be a positive number. Terminating.\n");

    If.RegisterInputKey("b_length", &celldm[1], 0.0, DBL_MAX, 0.0, 
                     CHECK_AND_TERMINATE, REQUIRED, 
                     "Second lattice constant.\n", 
                     "b_length must be a positive number. Terminating.\n");

    If.RegisterInputKey("c_length", &celldm[2], 0.0, DBL_MAX, 0.0, 
                     CHECK_AND_TERMINATE, REQUIRED, 
                     "Third lattice constant.\n", 
                     "c_length must be a positive number. Terminating.\n");

    If.RegisterInputKey("grid_spacing", &grid_spacing, 0.0, DBL_MAX, 0.35, 
                     CHECK_AND_TERMINATE, OPTIONAL, 
                     "Approximate grid spacing (bohr).\n", 
                     "grid_spacing must be a positive number. Terminating.\n");

    // Deault of zero is OK because this means to try to set it automatically later on.
    // The value of 64 covers any possible hardware scenario I can imagine currently but might
    // need to be adjusted at some point in the future.
    If.RegisterInputKey("threads_per_node", &lc.THREADS_PER_NODE, 0, 64, 0, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "Number of threads each MPI process will use. A value of 0 selects automatic setting.\n", 
                     "threads_per_node cannnot be a negative number and must be less than 64.\n");

    If.RegisterInputKey("potential_grid_refinement", &lc.FG_RATIO, 0, 5, 2, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "Ratio of the fine grid to the wavefunction grid.", 
                     "potential_grid_refinement must be in the range (1 <= ratio <= 4). Resetting to the default value of 2.\n");

    If.RegisterInputKey("potential_acceleration_constant_step", &lc.potential_acceleration_constant_step, 0.0, 2.0, 0.0, 
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

    If.RegisterInputKey("ionic_time_step_increase", &lc.iondt_inc, 1.0, 3.0, 1.1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Factor by which iondt is increased when dynamic timesteps are enabled.\n",
                     "ionic_time_step_increase must lie in the range (1.0,1.1). Resetting to the default value of 1.1.\n");
    
    If.RegisterInputKey("ionic_time_step_decrease", &lc.iondt_dec, 0.0, 1.0, 0.5,
                     CHECK_AND_FIX, OPTIONAL,
                     "Factor by which iondt is decreased when dynamic timesteps are enabled.\n",
                     "ionic_time_step_decrease must lie in the range (0.0,1.0). Resetting to the default value of 0.5.\n");

    If.RegisterInputKey("max_ionic_time_step", &lc.iondt_max, 0.0, 150.0, 150.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "Maximum ionic time step to use.\n",
                     "max_ionic_time_step must lie in the range (0.0,150.0). Resetting to the default value of 150.0.\n");

    If.RegisterInputKey("unoccupied_states_per_kpoint", &lc.num_unocc_states, 10, INT_MAX, 10, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "The number of unoccupied orbitals.\n", 
                     "unoccupied_states_per_kpoint must be at least 10. Fixing.\n");

    If.RegisterInputKey("state_block_size", &lc.state_block_size, 1, INT_MAX, 64, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "state_block used in nlforce.\n", 
                     "it is better to be 2^n.\n");


    If.RegisterInputKey("extra_random_lcao_states", &lc.extra_random_lcao_states, 0, INT_MAX, 0, 
                     CHECK_AND_TERMINATE, OPTIONAL, 
                     "Extra random wavefunctions to use for LCAO starts.\n", 
                     "extra_random_lcao_states must be greater than 0. Terminating.\n");

    If.RegisterInputKey("system_charge", &lc.background_charge, -DBL_MAX, DBL_MAX, 0.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of excess holes in the system (useful for doped systems). Example: 2 means system is missing two electrons\n",
                     "system_charge must be a real number.\n");

    If.RegisterInputKey("occupation_electron_temperature_eV", &lc.occ_width, 0.0, 2.0, 0.04,
                     CHECK_AND_FIX, OPTIONAL,
                     "Target electron temperature when not using fixed occupations.\n ",
                     "occupation_electron_temperature_eV must lie in the range (0.0,2.0). Resetting to the default value of 0.04.\n");

    If.RegisterInputKey("occupation_number_mixing", &lc.occ_mix, 0.0, 1.0, 1.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "Mixing parameter for orbital occupations when not using fixed occupations.\n",
                     "occupation_number_mixing must lie in the range (0.0,1.0). Resetting to the default value of 0.3.\n");
    If.RegisterInputKey("MP_order", &lc.mp_order, 0, 5, 2, 
                     CHECK_AND_FIX, OPTIONAL, 
                     "order of Methefessel Paxton occupation\n", 
                     "0 means simple error function as distribution\n");

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
    If.RegisterInputKey("tddft_steps", &lc.tddft_steps, 0, INT_MAX, 2000,
                     CHECK_AND_FIX, OPTIONAL, 
                     "Maximum number of tddft steps to perform.\n", 
                     "tddft steps must be greater than 0. Resetting to the default value of 2000\n");

    If.RegisterInputKey("charge_pulay_order", &lc.charge_pulay_order, 1, 10, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_pulay_scale", &lc.charge_pulay_scale, 0.0, 1.0, 0.50,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "charge_pulay_scale must lie in the range (0.0,1.0). Resetting to the default value of 0.50\n");

    If.RegisterInputKey("unoccupied_tol_factor", &lc.unoccupied_tol_factor, 0.000001, 10000.0, 1000.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "unoccupied_tol_factor must lie in the range (0.000001,10000.0). Resetting to the default value of 1000.0\n");

    If.RegisterInputKey("charge_pulay_refresh", &lc.charge_pulay_refresh, 1, INT_MAX, 100,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_broyden_order", &lc.charge_broyden_order, 1, 10, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_broyden_scale", &lc.charge_broyden_scale, 0.0, 1.0, 0.50,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "charge_broyden_scale must lie in the range (0.0,1.0). Resetting to the default value of 0.50\n");

    If.RegisterInputKey("projector_expansion_factor", &lc.projector_expansion_factor, 0.5, 3.0, 1.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "projector_expansion_factor must lie in the range (0.5,3.0). Resetting to the default value of 1.0\n");

    If.RegisterInputKey("write_data_period", &lc.checkpoint, 5, 50, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("write_eigvals_period", &lc.write_eigvals_period, 1, 100, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "How often to output eigenvalues.",
                     "write_eigvals_period must lie in the range (1,100). Resetting to the default value of 5.\n");

    If.RegisterInputKey("max_md_steps", &lc.max_md_steps, 0, INT_MAX, 100,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "Maximum number of molecular dynamics steps to perform.",
                     "max_md_steps must be a positive value. Terminating.\n");

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

    If.RegisterInputKey("kohn_sham_mucycles", &lc.eig_parm.mucycles, 1, 6, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of mu (also known as W) cycles to use in the kohn-sham multigrid preconditioner.\n",
                     "kohn_sham_mucycles must lie in the range (1,3). Resetting to the default value of 1.\n");

    If.RegisterInputKey("kohn_sham_fd_order", &lc.kohn_sham_fd_order, 4, 12, 8,
                     CHECK_AND_FIX, OPTIONAL,
                     "Order of the global grid finite difference operators to be used in the kohn-sham multigrid preconditioner.\n ",
                     "kohn_sham_fd_order must lie in the range (4,12). Resetting to the default value of 8.\n");

    If.RegisterInputKey("force_grad_order", &lc.force_grad_order, 0, 12, 8,
                     CHECK_AND_FIX, OPTIONAL,
                     "Order of the global grid finite difference operators to be used in the kohn-sham multigrid preconditioner.\n ",
                     "kohn_sham_fd_order must lie in the range (4,12). Resetting to the default value of 8.\n");

    If.RegisterInputKey("kohn_sham_coarse_time_step", &lc.eig_parm.sb_step, 0.0, 1.2, 1.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "Time step to use in the kohn-sham multigrid solver on the coarse levels.\n",
                     "kohn_sham_coarse_time_step must lie in the range (0.5,1.2). Resetting to the default value of 1.0.\n");

    If.RegisterInputKey("kohn_sham_time_step", &lc.eig_parm.gl_step, 0.0, 2.0, 0.66,
                     CHECK_AND_FIX, OPTIONAL,
                     "Smoothing timestep to use on the fine grid in the the kohn-sham multigrid preconditioner.\n",
                     "kohn_sham_time_step must lie in the range (0.4,2.0). Resetting to the default value of 0.66.\n");

    If.RegisterInputKey("kohn_sham_mg_timestep", &lc.eig_parm.mg_timestep, 0.0, 2.0, 0.6666666666666,
                     CHECK_AND_FIX, OPTIONAL,
                     "timestep for multigrid correction.\n",
                     "kohn_sham_mg_step must lie in the range (0.0,2.0). Resetting to the default value of 0.66\n");

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

    If.RegisterInputKey("poisson_finest_time_step", &lc.poi_parm.gl_step, 0.4, 0.8, 0.6,
                     CHECK_AND_FIX, OPTIONAL,
                     "Time step to use in the poisson multigrid solver on the finest level.\n",
                     "poisson_finest_time_step must lie in the range (0.4,0.8). Resetting to the default value of 0.6.\n");

    If.RegisterInputKey("poisson_coarse_time_step", &lc.poi_parm.sb_step, 0.4, 0.8, 0.6,
                     CHECK_AND_FIX, OPTIONAL,
                     "Time step to use in the poisson multigrid solver on the coarse levels.\n",
                     "poisson_coarse_time_step must lie in the range (0.4,0.8). Resetting to the default value of 0.6.\n");

    If.RegisterInputKey("poisson_coarsest_steps", &lc.poi_parm.coarsest_steps, 10, 100, 25,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of smoothing steps to use on the coarsest level in the hartree multigrid solver.\n",
                     "poisson_coarsest_steps must lie in the range (10,100). Resetting to the default value of 25.\n");

    If.RegisterInputKey("kohn_sham_mg_levels", &lc.eig_parm.levels, -1, 6, -1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of multigrid levels to use in the kohn-sham multigrid preconditioner.\n",
                     "kohn_sham_mg_levels must lie in the range (-1,6) where -1=automatic. Resetting to the default value of automatic (-1).\n");

    If.RegisterInputKey("poisson_mg_levels", &lc.poi_parm.levels, -1, 6, -1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of multigrid levels to use in the hartree multigrid solver.\n",
                     "poisson_mg_levels must lie in the range (-1,6) where -1=automatic. Resetting to the default value of automatic (-1).\n");

    If.RegisterInputKey("fine_grid_non_local_pp", &lc.nxfgrid, 1, 8, 4,
                     CHECK_AND_FIX, OPTIONAL,
                     "Fine grid for non-local pseudopotential.\n",
                     "fine_grid_non_local_pp must lie in the range (1,4). Resetting to the default value of 4.\n");

    If.RegisterInputKey("scalapack_block_factor", &lc.scalapack_block_factor, 4, 512,64,
                     CHECK_AND_FIX, OPTIONAL,
                     "Block size to use with scalapack. Optimal value is dependent on matrix size and system hardware.\n",
                     "scalapack_block_factor must lie in the range (4,512). Resetting to the default value of 64.\n");

    If.RegisterInputKey("non_local_block_size", &lc.non_local_block_size, 64, 2048, 512,
                     CHECK_AND_FIX, OPTIONAL,
                     "Block size to use when applying the non-local and S operators.\n",
                     "non_local_block_size must lie in the range (64,2048). Resetting to the default value of 512.\n");

    If.RegisterInputKey("cublasxt_block_size", &lc.cublasxt_block_size, 512, 10240, 2048,
                     CHECK_AND_FIX, OPTIONAL,
                     "Block size for cublasxt routines.\n",
                     "cublasxt_block_size must lie in the range (512,10240). Resetting to the default value of 512.\n");

    If.RegisterInputKey("E_POINTS", &lc.E_POINTS, 201, 201, 201,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("max_rmg_steps", &lc.max_rmg_steps, 0, 1000, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Number of rmg \"restart like\" (NEB/exchange/ARTs) steps to perform",
                     "max_rmg_steps must lie in the range (0,1000). Resetting to the default value of 1.");

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

    If.RegisterInputKey("b_spline_order", &lc.interp_order, 0, 7, 5,
                     CHECK_AND_FIX, OPTIONAL,
                     "Order of interpolation to use if b-spline is the selected method.\n",
                     "b_spline_order must lie in the range (0,7). Resetting to the default value of 5.\n");

    If.RegisterInputKey("b_spline_trade_order", &lc.interp_trade, 3, 3, 3,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("charge_density_mixing", &lc.mix, 0.0, 1.0, 0.5,
                     CHECK_AND_FIX, OPTIONAL,
                     "Proportion of the current charge density to replace with the new density after each scf step when linear mixing is used.\n",
                     "charge_density_mixing must lie in the range (0.0, 1.0) Resetting to the default value of 0.5.\n");

    If.RegisterInputKey("folded_spectrum_width", &lc.folded_spectrum_width, 0.15, 1.0, 0.3,
                     CHECK_AND_FIX, OPTIONAL,
                     "Submatrix width to use as a fraction of the full spectrum.\n",
                     "folded_spectrum_width must lie in the range (0.15,1.0). Resetting to the default value of 0.3.\n");

    If.RegisterInputKey("charge_pulay_special_metrics_weight", &lc.charge_pulay_special_metrics_weight, -DBL_MAX, DBL_MAX, 100.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "charge_pulay_special_metrics_weight must be a real number.");

    //RMG2BGW options
    If.RegisterInputKey("use_symmetry", &lc.is_use_symmetry, false, 
                        "For non-gamma point, always true, for gamma point, optional");
    If.RegisterInputKey("rmg2bgw", &lc.rmg2bgw, false, 
                        "Write wavefunction in G-space to BerkeleyGW WFN file.");
    If.RegisterInputKey("ecutrho", &lc.ecutrho, 0.0, 10000.0, -100.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "ecut for rho in unit of Ry.\n",
                     "\n");
    If.RegisterInputKey("ecutwfc", &lc.ecutwfc, 0.0, 10000.0, 40.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "ecut for rho in unit of Ry.\n",
                     "\n");

    If.RegisterInputKey("vxc_diag_nmin", &lc.vxc_diag_nmin, 1, 10000, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Minumum band index for diagonal Vxc matrix elements.\n ",
                     "vxc_diag_nmin must lie in the range (1, 10000). Resetting to the default value of 1.\n");

    If.RegisterInputKey("vxc_diag_nmax", &lc.vxc_diag_nmax, 1, 10000, 1,
                     CHECK_AND_FIX, OPTIONAL,
                     "Maximum band index for diagonal Vxc matrix elements.\n ",
                     "vxc_diag_nmax must lie in the range (1, 10000). Resetting to the default value of 1.\n");

    // Booleans next. Booleans are never required.
    If.RegisterInputKey("coalesce_states", &lc.coalesce_states, false,
                        "Flag indicating whether or not to coalesce states.");

    If.RegisterInputKey("localize_projectors", &lc.localize_projectors, true,
                        "Flag indicating whether or not the non-local projectors should be localized.");

    If.RegisterInputKey("localize_localpp", &lc.localize_localpp, true,
                        "Flag indicating whether or not the local pseudopotential should be localized.");

    If.RegisterInputKey("charge_pulay_special_metrics", &lc.charge_pulay_special_metrics, false,
                        "Flag to test whether or not the modified metrics should be used in Pulay mixing.");

    If.RegisterInputKey("write_pseudopotential_plots", NULL, false,
                        "Flag to indicate whether or not to write pseudopotential plots.\n");

    If.RegisterInputKey("equal_initial_density", &lc.init_equal_density_flag, false,
                        "Specifies whether to set initial up and down density to be equal.");

    If.RegisterInputKey("write_pdos", &lc.pdos_flag, false,
                        "Flag to write partial density of states.");

    If.RegisterInputKey("sort_wavefunctions", &lc.sortflag, false, 
                        "Sort wavefunctions by eigenvalue. Not needed if using subspace diagonalization.");

    If.RegisterInputKey("initial_diagonalization", &lc.initdiag, true, 
                        "Perform initial subspace diagonalization.");
    
    If.RegisterInputKey("verbose", &lc.verbose, false,
                        "Flag for writing out extra information\n");

    If.RegisterInputKey("folded_spectrum", &lc.use_folded_spectrum, false, 
                         "Use folded spectrum.");

    If.RegisterInputKey("use_numa", &lc.use_numa, true, 
                         "Use internal numa setup if available.");

    If.RegisterInputKey("mpi_queue_mode", &lc.mpi_queue_mode, false, 
                         "Use mpi queue mode.");

    If.RegisterInputKey("require_huge_pages", &lc.require_huge_pages, false, 
                         "If set RMG assumes that sufficient huge pages are available. Bad results may occur if that is not true.");

    If.RegisterInputKey("relax_dynamic_timestep", NULL, false,
                        "Flag indicating whether or not to use dynamic timesteps in relaxation mode.\n");

    If.RegisterInputKey("freeze_occupied", NULL, false,
                        "Flag indicating whether or not to freeze the density and occupied orbitals after a restart.\n");

    If.RegisterInputKey("relax_max_force", &lc.thr_frc, 0.0, DBL_MAX, 2.5E-3,
                     CHECK_AND_FIX, OPTIONAL,
                     "Force value at which an ionic relaxation is considered to be converged.\n",
                     "relax_max_force must be a positive value. Resetting to default value of 2.5e-03.\n");

    If.RegisterInputKey("md_randomize_velocity", &lc.nose.randomvel, true,
                        "");

    If.RegisterInputKey("scalapack_global_sums", &lc.scalapack_global_sums, true,
                        "");

    If.RegisterInputKey("output_rho_xsf", NULL, false,
                        "Generate xsf format for electronic density.");

    If.RegisterInputKey("rms_convergence_criterion", &lc.thr_rms, 0.0, 1.0e-3, 1.0e-7,
                     CHECK_AND_FIX, OPTIONAL,
                     "The RMS value of the change in the total potential where we assume self consistency has been achieved.\n",
                     "rms_convergence_criterion must lie in the range (1.0e-04,1.0e-14). Resetting to default value of 1.0e-7.\n");

    If.RegisterInputKey("energy_convergence_criterion", &lc.thr_energy, 1.0e-20, 1.0e-8, 1.0e-13,
                     CHECK_AND_FIX, OPTIONAL,
                     "The RMS value of the change in the total potential where we assume self consistency has been achieved.\n",
                     "rms_convergence_criterion must lie in the range (1.0e-04,1.0e-14). Resetting to default value of 1.0e-7.\n");

    If.RegisterInputKey("preconditioner_threshold", &lc.preconditioner_thr, 1.0e-9, 1.0e-1, 1.0e-1,
                     CHECK_AND_FIX, OPTIONAL,
                     "The RMS value of the change in the total potential where we switch the preconditioner from single to double precision.\n",
                     "preconditioner_threshold must lie in the range (1.0e-9,1.0e-1). Resetting to default value of 1.0e-1.\n");

    If.RegisterInputKey("gw_residual_convergence_criterion", &lc.gw_threshold, 1.0e-14, 4.0e-4, 1.0e-6,
                     CHECK_AND_FIX, OPTIONAL,
                     "The max value of the residual for unoccupied orbitals when performing a GW calculation.\n",
                     "gw_residual_convergence_criterion must lie in the range (4.0e-04,1.0e-14). Resetting to default value of 4.0e-04.\n");

    If.RegisterInputKey("gw_residual_fraction", &lc.gw_residual_fraction, 0.0, 1.0, 0.90,
                     CHECK_AND_FIX, OPTIONAL,
                     "The residual value specified by gw_residual_convergence_criterion is applied to this fraction of the total spectrum.\n",
                     "gw_residual_fraction must lie in the range (0.0,1.0). Resetting to default value of 0.90.\n");

    If.RegisterInputKey("hartree_rms_ratio", &lc.hartree_rms_ratio, 1000.0, DBL_MAX, 100000.0,
                     CHECK_AND_FIX, OPTIONAL,
                     "Ratio between target RMS for get_vh and RMS total potential.\n",
                     "hartree_rms_ratio must be in the range (1000.0, 1000000.0). Resetting to default value of 100000.0.\n");

    If.RegisterInputKey("electric_field_magnitude", &lc.e_field, 0.0, DBL_MAX, 0.0,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "Magnitude of external electric field.\n",
                     "electric_field_magnitude must be a postive value.\n");

    Ri::ReadVector<double> def_electric_field({{0.0,0.0,1.0}});
    Ri::ReadVector<double> electric_field;
    If.RegisterInputKey("electric_field_vector", &electric_field, &def_electric_field, 3, OPTIONAL,
                     "Components of the electric field.\n",
                     "You must specify a triplet of (X,Y,Z) dimensions for the electric field vector.\n");

    If.RegisterInputKey("Emin", &lc.Emin, -100.0, 100.0, -6.0,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("Emax", &lc.Emax, -100.0, 100.0, 0.0,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "",
                     "");

    If.RegisterInputKey("neb_spring_constant", &lc.neb_spring_constant, 0.05, 3.0, 0.5,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "",
                     "neb_spring_constant must be in the range (0.05, 3.0).\n");

    If.RegisterInputKey("energy_cutoff_parameter", &lc.cparm, 0.6, 1.0, 0.8,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "energy_cutoff_parameter must be in the range (0.6,1.0). Resetting to default value of 0.8.\n");


    std::string Occup, Occdown;
    std::string Occ;
//    if(lc.spin_flag) {

        If.RegisterInputKey("states_count_and_occupation_spin_up", &Occup, "",
                         CHECK_AND_FIX, OPTIONAL,
                         "Occupation string for spin up states.\n",
                         "");

        If.RegisterInputKey("states_count_and_occupation_spin_down", &Occdown, "",
                         CHECK_AND_FIX, OPTIONAL,
                         "Occupation string for spin down states.\n",
                         "");

        
//    }
//    else {

        If.RegisterInputKey("states_count_and_occupation", &Occ, "",
                         CHECK_AND_FIX, OPTIONAL,
                         "Occupation string for states.\n",
                         "");

//    }

    If.RegisterInputKey("kpoint_distribution", &pelc.pe_kpoint, -INT_MAX, INT_MAX, -1,
                     CHECK_AND_FIX, OPTIONAL,
                     "",
                     "");

    If.LoadInputKeys();

    // Check items that require custom handling
    // Some hacks here to deal with code branches that are still in C
    if(!Description.length()) Description = "RMG electronic structure calculation.";
    std::strncpy(lc.description, Description.c_str(), sizeof(lc.description));

    if(!Infile.length()) Infile = "Waves/wave.out";
    std::strncpy(lc.infile, Infile.c_str(), sizeof(lc.infile));

    if(!Outfile.length()) Outfile = "Waves/wave.out";
    std::strncpy(lc.outfile, Outfile.c_str(), sizeof(lc.outfile));

    if(lc.outfile[0] !='/') 
    {
        char *temp = new char[255];
        snprintf(temp, 255, "%s%s", pct.image_path[pct.thisimg], lc.outfile);
        std::strncpy(lc.outfile, temp, sizeof(lc.outfile));
        delete [] temp;
    }
    if(lc.infile[0] !='/') 
    {
        char *temp = new char[255];
        snprintf(temp, 255, "%s%s", pct.image_path[pct.thisimg], lc.infile);
        std::strncpy(lc.infile, temp, sizeof(lc.infile));
        delete [] temp;
    }

    if(!Infile_tddft.length()) Infile = "Waves/wave_tddft.out";
    std::strncpy(lc.infile_tddft, Infile_tddft.c_str(), sizeof(lc.infile_tddft));

    if(!Outfile_tddft.length()) Outfile = "Waves/wave_tddft.out";
    std::strncpy(lc.outfile_tddft, Outfile_tddft.c_str(), sizeof(lc.outfile_tddft));

    if(lc.outfile_tddft[0] !='/') 
    {
        char *temp = new char[255];
        snprintf(temp, 255, "%s%s", pct.image_path[pct.thisimg], lc.outfile_tddft);
        std::strncpy(lc.outfile_tddft, temp, sizeof(lc.outfile_tddft));
        delete [] temp;
    }
    if(lc.infile[0] !='/') 
    {
        char *temp = new char[255];
        snprintf(temp, 255, "%s%s", pct.image_path[pct.thisimg], lc.infile_tddft);
        std::strncpy(lc.infile_tddft, temp, sizeof(lc.infile_tddft));
        delete [] temp;
    }




    if(lc.spin_flag) {

        std::strncpy(lc.occupation_str_spin_up, Occup.c_str(), sizeof(lc.occupation_str_spin_up));
        std::strncpy(lc.occupation_str_spin_down, Occdown.c_str(), sizeof(lc.occupation_str_spin_down));

    }
    else {

        std::strncpy(lc.occupation_str, Occ.c_str(), sizeof(lc.occupation_str));

    }


    try {
       for(int ix = 0;ix < 3;ix++) {
           lc.kpoint_mesh[ix] = kpoint_mesh.vals.at(ix);
           lc.kpoint_is_shift[ix] = kpoint_is_shift.vals.at(ix);
       }
    }
    catch (const std::out_of_range& oor) {
        throw RmgFatalException() << "You must specify a triplet of (X,Y,Z) dimensions for kpoint_mesh and kpoint_is_shift.\n";
    }


    // Currently, fine grid has to be the same in each direction
    lc.nzfgrid = lc.nyfgrid = lc.nxfgrid;

    // This is an ugly hack but it's necessary for now since we want to print the
    // original options back out for reuse and the code modifies these
    static double orig_celldm[3];
    orig_celldm[0] = celldm[0];
    orig_celldm[1] = celldm[1];
    orig_celldm[2] = celldm[2];
    InputKey *ik = InputMap["a_length"];
    ik->Readdoubleval = &orig_celldm[0];
    ik = InputMap["b_length"];
    ik->Readdoubleval = &orig_celldm[1];
    ik = InputMap["c_length"];
    ik->Readdoubleval = &orig_celldm[2];


    // Transform to atomic units, which are used internally if input is in angstrom 
    if (Verify ("crds_units", "Angstrom", InputMap))
    {
        celldm[0] *= A_a0;
        celldm[1] *= A_a0;
        celldm[2] *= A_a0;
    }


    // Here we read celldm as a,b,c but for most lattice types code uses a, b/a, c/a 
    // Every lattice type uses a, b/a, c/a except CUBIC_PRIMITIVE, CUBIC_FC and CUBIC_BC 
    if (!Verify ("bravais_lattice_type", "Cubic Face Centered", InputMap) &&
        !Verify ("bravais_lattice_type", "Cubic Body Centered", InputMap))
    {
        celldm[1] /= celldm[0];
        celldm[2] /= celldm[0];
    }

    // Lattice vectors are orthogonal except for Hex which is setup inside latgen
    celldm[3] = 0.0;
    celldm[4] = 0.0;
    celldm[5] = 0.0;

    // Set up the lattice vectors
    Rmg_L.set_ibrav_type(ibrav);
    int flag = 0;
    Rmg_L.latgen(celldm, &omega, a0, a1, a2, &flag);


    int NX_GRID = WavefunctionGrid.vals.at(0);
    int NY_GRID = WavefunctionGrid.vals.at(1);
    int NZ_GRID = WavefunctionGrid.vals.at(2);

    lc.dipole_corr[0] = DipoleCorrection.vals.at(0);
    lc.dipole_corr[1] = DipoleCorrection.vals.at(1);
    lc.dipole_corr[2] = DipoleCorrection.vals.at(2);

    CheckAndTerminate(NX_GRID, 1, INT_MAX, "The value given for the global wavefunction grid X dimension is " + boost::lexical_cast<std::string>(NX_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(NY_GRID, 1, INT_MAX, "The value given for the global wavefunction grid Y dimension is " + boost::lexical_cast<std::string>(NY_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(NZ_GRID, 1, INT_MAX, "The value given for the global wavefunction grid Z dimension is " + boost::lexical_cast<std::string>(NZ_GRID) + " and only postive values are allowed.");


    /* read the electric field vector */
    try {
        ct.x_field_0 = electric_field.vals.at(0);
        ct.y_field_0 = electric_field.vals.at(1);
        ct.z_field_0 = electric_field.vals.at(2);
    }
    catch (const std::out_of_range& oor) {
        throw RmgFatalException() << "You must specify a triplet of (X,Y,Z) values for the electric field vector.\n";
    }


    if (lc.iondt_max < lc.iondt)
        throw RmgFatalException() << "max_ionic_time_step " << lc.iondt_max << " has to be >= than ionic_time_step " << ct.iondt << "\n";

    // Constraints for Neb relax. What does this magic number mean? (set constraint type for switch in Common/constrain.c)
    if(Verify("calculation_mode", "NEB Relax", InputMap)) lc.constrainforces = 5;

    // Background charge is defined to be the opposite of system charge
    lc.background_charge *= -1.0;

    lc.occ_width *= eV_Ha;
    lc.e_field *= eV_Ha;

    // Potential acceleration must be disabled if freeze_occupied is true
    if(Verify ("freeze_occupied", true, InputMap) ) {
        lc.potential_acceleration_constant_step = 0.0;
        lc.potential_acceleration_poisson_step = 0.0;
        std::cout << "You have set freeze_occupied=true so potential acceleration is disabled." << std::endl;
    }

    // Debug code
    //for(auto it = InputMap.begin();it != InputMap.end(); ++it) {
    //    std::pair<std::string, InputKey*> Check = *it;
    //    InputKey *CheckKey = it->second;
        //std::cout << Check.first << " = " << CheckKey->Print() << std::endl;
    //}

    // Set up energy output units
    static char *Hartree_str = "Ha";
    static char *Rydberg_str = "Ry";
    lc.energy_output_conversion[0] = 1.0;
    lc.energy_output_conversion[1] = 2.0;
    lc.energy_output_string[0] = Hartree_str;
    lc.energy_output_string[1] = Rydberg_str;

    ct.use_vdwdf_finegrid = Verify ("vdwdf_grid_type", "Fine", InputMap);
}
