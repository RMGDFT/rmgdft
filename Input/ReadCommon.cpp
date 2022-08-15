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
#include "BaseGrid.h"
#include "transition.h"

#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "CheckValue.h"
#include "RmgException.h"
#include "RmgInputFile.h"
#include "InputOpts.h"
#include "rmg_error.h"



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

// Auxiliary function used to convert input file names into a full path
void MakeFullPath(char *fullpath, PE_CONTROL& pelc)
{
    if(fullpath[0] !='/')
    {
        char temp[MAX_PATH];
        snprintf(temp, sizeof(temp) - 1, "%s%s", pelc.image_path[pelc.thisimg], fullpath);
        std::strncpy(fullpath, temp, MAX_PATH);
    }
}
void MakeFullPath(std::string &fullpath, PE_CONTROL& pelc)
{
    if(fullpath[0] != '/')
    {
        fullpath = std::string(pelc.image_path[pelc.thisimg]) + fullpath;
    }
}


namespace Ri = RmgInput;

void ReadCommon(char *cfile, CONTROL& lc, PE_CONTROL& pelc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    RmgInputFile If(cfile, InputMap, pelc.img_comm);
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
    std::string Weightsfile;
    std::string Qfunctionfile;
    std::string Workfile;
    std::string Orbitalfile;
    std::string ExxIntfile;
    std::string PseudoPath;
    std::string VdwKernelfile;
 
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

    static Ri::ReadVector<int> Cell_movable;
    Ri::ReadVector<int> DefCell_movable({{0,0,0,0,0,0,0,0,0}});

    double celldm[6] = {1.0,1.0,1.0,0.0,0.0,0.0};
    static double grid_spacing;
    double a0[3], a1[3], a2[3], omega;

 
    If.RegisterInputKey("description", &Description, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Description of the run. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("pseudopotential", NULL , "",
                     CHECK_AND_FIX, OPTIONAL,
                     "External pseudopotentials may be specfied with this input key. The format uses the atomic symbol followed by the pseudopotential file name.    pseudopotential = \"Ni Ni.UPF  O O.UPF\"", 
                     "", PSEUDO_OPTIONS);

    If.RegisterInputKey("Hubbard_U", NULL , "",
                     CHECK_AND_FIX, OPTIONAL,
"Hubbard U parameter for each atomic species using the format &"
"Hubbard_U=\"Ni  6.5  3d  0.0  0.0  0.0\"",
                     "LDA+U prameters: elements, U value, orbital label, J[0-3]", LDAU_OPTIONS);

    If.RegisterInputKey("input_wave_function_file", &Infile, "Waves/wave.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Input file/path to  read wavefunctions and other binary data from on a restart. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("input_tddft_file", &Infile_tddft, "Waves/wave_tddft.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Input file/path to  read wavefunctions and other binary data from on a restart. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("nvme_weights_filepath", &Weightsfile, "Weights/",
                     CHECK_AND_FIX, OPTIONAL,
                     "File/path for disk storage of projector weights. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("nvme_work_filepath", &Workfile, "Work/",
                     CHECK_AND_FIX, OPTIONAL,
                     "File/path for disk storage of workspace. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("nvme_orbitals_filepath", &Orbitalfile, "Orbitals/",
                     CHECK_AND_FIX, OPTIONAL,
                     "File/path for runtime disk storage of orbitals. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("qfunction_filepath", &Qfunctionfile, "Qfunctions/",
                     CHECK_AND_FIX, OPTIONAL,
                     "File/path for runtime disk storage of qfunctions. ",
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("exx_integrals_filepath", &ExxIntfile, "afqmc_rmg",
                     CHECK_AND_FIX, OPTIONAL,
                     "File/path for exact exchange integrals. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("vdwdf_kernel_filepath", &VdwKernelfile, "vdW_kernel_table",
                     CHECK_AND_FIX, OPTIONAL,
                     "File/path for vdW_kernel_table data. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("pseudo_dir", &PseudoPath, ".",
                     CHECK_AND_FIX, OPTIONAL,
                     "Directory where pseudopotentials are stored. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("output_wave_function_file", &Outfile, "Waves/wave.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Output file/path to store wavefunctions and other binary data. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("output_tddft_file", &Outfile_tddft, "Waves/wave_tddft.out",
                     CHECK_AND_FIX, OPTIONAL,
                     "Output file/path to store wavefunctions and other binary data. ", 
                     "", CONTROL_OPTIONS);

    If.RegisterInputKey("tddft_mode", NULL, &lc.tddft_mode, "electric field",
                     CHECK_AND_TERMINATE, OPTIONAL,tddft_mode, 
                     "TDDFT mode ", 
                     "either electric field kick or point charge kick", TDDFT_OPTIONS);

    If.RegisterInputKey("restart_tddft", &lc.restart_tddft, false, 
                        "restart TDDFT", TDDFT_OPTIONS);

    If.RegisterInputKey("stress", &lc.stress, false, 
                        "flag to control stress cacluation", CONTROL_OPTIONS);

    If.RegisterInputKey("cell_relax", &lc.cell_relax, false, 
                        "flag to control unit cell relaxation", CONTROL_OPTIONS);

    If.RegisterInputKey("processor_grid", &ProcessorGrid, &DefProcessorGrid, 3, OPTIONAL, 
                     "Three-D (x,y,z) layout of the MPI processes. ", 
                     "You must specify a triplet of (X,Y,Z) dimensions for the processor grid. ", CELL_OPTIONS);

    If.RegisterInputKey("wavefunction_grid", &WavefunctionGrid, &DefWavefunctionGrid, 3, OPTIONAL, 
                     "Three-D (x,y,z) dimensions of the grid the wavefunctions are defined on. ", 
                     "You must specify a triplet of (X,Y,Z) dimensions for the wavefunction grid. ", CELL_OPTIONS);

    If.RegisterInputKey("dipole_correction", &DipoleCorrection, &DefDipoleCorrection, 3, OPTIONAL, 
                     "(1,1,1) for molecule, dipole correction in all directions.  ", 
                     "(0,0,0) means no correction by default, (1,0,0) or others have not programed ");

    If.RegisterInputKey("cell_movable", &Cell_movable, &DefCell_movable, 9, OPTIONAL, 
                     "9 numbers to control cell relaxation  ", 
                     "0 0 0 0 0 0 0 0 0 by default, no cell relax ", CELL_OPTIONS);

    If.RegisterInputKey("kpoint_mesh", &kpoint_mesh, &def_kpoint_mesh, 3, OPTIONAL, 
                     "Three-D layout of the kpoint mesh. ", 
                     "You must specify a triplet of coordinate dimensions for the kpoint_mesh. ", CELL_OPTIONS);

    If.RegisterInputKey("kpoint_is_shift", &kpoint_is_shift, &def_kpoint_is_shift, 3, OPTIONAL, 
                     "Three-D layout of the kpoint shift. ", 
                     "You must specify a triplet of coordinate dimensions for kpoint_is_shift. ", CELL_OPTIONS);

    int ibrav;
    If.RegisterInputKey("bravais_lattice_type", NULL, &ibrav, "None",
                     CHECK_AND_TERMINATE, OPTIONAL, bravais_lattice_type,
                     "Bravais Lattice Type. ", 
                     "bravais_lattice_type not found. ", CELL_OPTIONS);

    If.RegisterInputKey("vdw_corr", NULL, &lc.vdw_corr, "None",
                     CHECK_AND_TERMINATE, OPTIONAL, vdw_corr,
                     "Type of vdw correction  ", 
                     "Only Grimme D2 is implemented now ", CONTROL_OPTIONS);

    If.RegisterInputKey("dftd3_version", &lc.dftd3_version, 2, 6, 3, 
            CHECK_AND_FIX, OPTIONAL, 
            "Grimme's DFT-D3 versions,", 
            "version 2: DFT-D2, version 3: DFT-D3 zero damp, version 4: DFT-D3 BJ damp. ", MISC_OPTIONS);

    If.RegisterInputKey("start_mode", NULL, &lc.runflag, "LCAO Start",
                     CHECK_AND_TERMINATE, OPTIONAL, start_mode,
                     "Type of run. ", 
                     "start_mode must be one of  \"Random Start\", \"Restart From File\", or \"LCAO Start\". Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("atomic_orbital_type", NULL, &ct.atomic_orbital_type, "delocalized",
                     CHECK_AND_TERMINATE, OPTIONAL, atomic_orbital_type,
                     "Atomic Orbital Type. Choices are localized and delocalized. ", 
                     "atomic_orbital_type not found. ", PSEUDO_OPTIONS);

    If.RegisterInputKey("internal_pseudo_type", NULL, &ct.internal_pseudo_type, "nc_accuracy",
                     CHECK_AND_TERMINATE, OPTIONAL, internal_pseudo_type,
                     "Internal pseudopotential type. Choices are ultrasoft, nc_accuracy ", 
                     "internal pseudopotential type not found. ", PSEUDO_OPTIONS);

    If.RegisterInputKey("subdiag_driver", NULL, &lc.subdiag_driver, "auto",
                     CHECK_AND_FIX, OPTIONAL, subdiag_driver,
                     "Driver type used for subspace diagonalization of the eigenvectors. ", 
                     "subdiag_driver must be lapack, scalapack, cusolver or auto. Resetting to auto. ", DIAG_OPTIONS);

    If.RegisterInputKey("kohn_sham_solver", NULL, &lc.kohn_sham_solver, "davidson",
                     CHECK_AND_FIX, OPTIONAL, kohn_sham_solver,
"RMG supports a pure multigrid Kohn-Sham solver as well as "
"a multigrid preconditioned davidson solver. The davidson "
"solver is usually better for smaller problems with the pure "
"multigrid solver often being a better choice for very large "
"problems.",
                     "kohn_sham_solver must be multigrid or davidson. Resetting to multigrid. ", KS_SOLVER_OPTIONS);

    If.RegisterInputKey("poisson_solver", NULL, &lc.poisson_solver, "pfft",
                     CHECK_AND_FIX, OPTIONAL, poisson_solver,
                     "poisson solver. ", 
                     "poisson_solver must be multigrid or pfft. Resetting to pfft. ", POISSON_OPTIONS);


    If.RegisterInputKey("crds_units", NULL, NULL, "Bohr",
                     CHECK_AND_FIX, OPTIONAL, crds_units,
                     "Units for the atomic coordinates. ", 
                     "Coordinates must be specified in either Bohr or Angstrom. ", CELL_OPTIONS);

    If.RegisterInputKey("lattice_units", NULL, NULL, "Bohr",
                     CHECK_AND_FIX, OPTIONAL, lattice_units,
                     "Units for the lattice vectors ", 
                     "lattice vectors' unit  must be specified in either Bohr or Angstrom, or Alat. ", CELL_OPTIONS);

    If.RegisterInputKey("charge_mixing_type", NULL, &lc.charge_mixing_type, "Pulay",
                     CHECK_AND_TERMINATE, OPTIONAL, charge_mixing_type,
"RMG supports Broyden, Pulay and Linear mixing "
"When the davidson Kohn-Sham solver is selected Broyden or "
"Pulay are preferred. For the multigrid solver Linear with "
"potential acceleration is often (but not always) the best "
"choice.",
                     "charge_mixing_type must be either \"Broyden\", \"Linear\" or \"Pulay\". Terminating. ", MIXING_OPTIONS);
    
    If.RegisterInputKey("ldau_mixing_type", NULL, &lc.ldau_mixing_type, "Linear",
                     CHECK_AND_TERMINATE, OPTIONAL, charge_mixing_type,
"RMG supports Pulay and Linear mixing for DFT+U occupation mixing ",
                     "ldau_mixing_type must be either \"Linear\" or \"Pulay\". Terminating. ", MIXING_OPTIONS);

    If.RegisterInputKey("charge_analysis", NULL, &lc.charge_analysis_type, "Voronoi",
                     CHECK_AND_TERMINATE, OPTIONAL, charge_analysis,
                     "Type of charge analysis to use. Only Voronoi deformation density is currently available. ", 
                     "charge_analysis must be either \"Voronoi\" or \"None\". Terminating. ");
    
    If.RegisterInputKey("charge_analysis_period", &lc.charge_analysis_period, 0, 500, 0,
                     CHECK_AND_FIX, OPTIONAL,
                     "How often to  perform and write out charge analysis.",
                     "charge_analysis_write_period must lie in the range (1,500). Resetting to the default value of 0. ", OUTPUT_OPTIONS);
    
    If.RegisterInputKey("dipole_moment", &lc.dipole_moment, false, 
                        "Turns on calculation of dipole moment for the entire cell.");

    If.RegisterInputKey("vdwdf_grid_type", NULL, NULL, "Coarse",
                     CHECK_AND_TERMINATE, OPTIONAL, vdwdf_grid_type,
                     "Type of grid to use when computing vdw-df correlation. ", 
                     "vdwdf_grid_type be either \"Coarse\" or \"Fine\". Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("relax_mass", NULL, &lc.relax_mass, "Atomic",
                     CHECK_AND_TERMINATE, OPTIONAL, relax_mass,
"Mass to use for structural relaxation, either atomic masses, or the mass of carbon for all atoms. ", 
                     "relax_mass must be either \"Atomic\" or \"Equal\". Terminating. ", MD_OPTIONS);

    If.RegisterInputKey("md_integration_order", NULL, &lc.mdorder, "5th Beeman-Velocity Verlet",
                     CHECK_AND_TERMINATE, OPTIONAL, md_integration_order,
                     "Integration order for molecular dynamics. ", 
                     "md_integration_order must be either \"2nd Velocity Verlet\", \"3rd Beeman-Velocity Verlet\" or \"5th Beeman-Velocity Verlet\". Terminating. ", MD_OPTIONS);

    If.RegisterInputKey("z_average_output_mode", NULL, &lc.zaverage, "None",
                     CHECK_AND_TERMINATE, OPTIONAL, z_average_output_mode,
                     "z_average_output_mode. ", 
                     "z_average_output_mode not supported. Terminating. ");

    If.RegisterInputKey("atomic_coordinate_type", NULL, &lc.crd_flag, "Absolute",
                     CHECK_AND_TERMINATE, OPTIONAL, atomic_coordinate_type,
                     "Flag indicated whether or not atomic coordinates are absolute or cell relative. ", 
                     "atomic_coordinate_type must be either \"Absolute\" or \"Cell Relative\". Terminating. ", CELL_OPTIONS);

    if(lc.forceflag != NEB_RELAX) 
        If.RegisterInputKey("calculation_mode", NULL, &lc.forceflag, "Quench Electrons",
                CHECK_AND_TERMINATE, OPTIONAL, calculation_mode,
                "Type of calculation to perform. ", 
                "calculation_mode not available. ", CONTROL_OPTIONS);

    If.RegisterInputKey("ldaU_mode", NULL, &lc.ldaU_mode, "None",
            CHECK_AND_TERMINATE, OPTIONAL, ldaU_mode,
            "Type of lda+u implementation. ", 
            "lda+u type not available. ", LDAU_OPTIONS);

    If.RegisterInputKey("relax_method", NULL, &lc.relax_method, "Fast Relax",
            CHECK_AND_TERMINATE, OPTIONAL, relax_method,
            "Type of relaxation method to use for structural optimizations. ", 
            "relax_method not supported. ", MD_OPTIONS);

    If.RegisterInputKey("md_temperature_control", NULL, &lc.tcontrol, "Nose Hoover Chains",
            CHECK_AND_TERMINATE, OPTIONAL, md_temperature_control,
            "Type of temperature control method to use in molecular dynamics. ", 
            "md_temperature_control type not supported. ", MD_OPTIONS);

    If.RegisterInputKey("md_temperature", &lc.nose.temp, 0.0, DBL_MAX, 300.0,
            CHECK_AND_FIX, OPTIONAL,
            "Target MD Temperature. ",
            "md_temperature must be a positive number. ", MD_OPTIONS);

    If.RegisterInputKey("md_nose_oscillation_frequency_THz", &lc.nose.fNose, 0.0, DBL_MAX, 15.59,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "md_nose_oscillation_frequency_THz must be a positive real number.", MD_OPTIONS);

    If.RegisterInputKey("energy_output_units", &EnergyOutputType, &lc.energy_output_units, "Hartrees",
            CHECK_AND_FIX, OPTIONAL, energy_output_units,
            "Units to be used when writing energy values to the output file. "
            " Hartrees or Rydbergs are available.", 
            "energy_output_units must be either \"Hartrees\" or \"Rydbergs\". Setting to \"Hartrees\". ", CONTROL_OPTIONS);

    If.RegisterInputKey("boundary_condition_type", NULL, &lc.boundaryflag, "Periodic",
            CHECK_AND_TERMINATE, OPTIONAL, boundary_condition_type,
            "Boundary condition type Only periodic is currently implemented. ", 
            "discretization_type must be Periodic. ");

    If.RegisterInputKey("exchange_correlation_type", NULL, &lc.xctype, "AUTO_XC",
            CHECK_AND_TERMINATE, OPTIONAL, exchange_correlation_type,
            "Most pseudopotentials specify the exchange correlation type they "
            "were generated with and the default value of AUTO_XC means that "
            "the type specified in the pseudopotial is what RMG will use. That "
            "can be overridden by specifying a value here.",
            "exchange_correlation_type not supported. Terminating. ", XC_OPTIONS);

    If.RegisterInputKey("dos_broading", &lc.gaus_broad, 0.0, 1.0, 0.1,
            CHECK_AND_FIX, OPTIONAL,
            "For DOS with Gaussian broading method",
            "in unit of eV ", OCCUPATION_OPTIONS);

    If.RegisterInputKey("dos_method", NULL, &lc.dos_flag, "tetrahedra",
            CHECK_AND_TERMINATE, OPTIONAL, dos_method,
            "tetrahedra or gauss smearing method for DOS calculation ",
            "dos_method not supported. Terminating. ", OCCUPATION_OPTIONS);

    If.RegisterInputKey("occupations_type", NULL, &lc.occ_flag, "Fermi Dirac",
            CHECK_AND_TERMINATE, OPTIONAL, occupations_type,
            "RMG supports several different ways of specifying orbital occupations. "
            "For a spin polarized system one may specify the occupations for up and "
            "down separately. In the case of a non-zero electronic temperature these "
            "will be adjusted as the calculation proceeds based on this setting. ",
            "occupations_type not supported. Terminating. ", OCCUPATION_OPTIONS);

    If.RegisterInputKey("interpolation_type", NULL, &lc.interp_flag, "FFT",
            CHECK_AND_TERMINATE, OPTIONAL, interpolation_type,
            "Interpolation method for transferring data between the potential grid "
            "and the wavefunction grid. Mostly for diagnostic purposes.", 
            "interpolation_type not supported. Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("exx_mode", NULL, &lc.exx_mode, "Local fft",
            CHECK_AND_TERMINATE, OPTIONAL, exx_mode,
            "FFT mode for exact exchange computations.",
            "exx mode not supported. Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("ExxIntCholosky", &lc.ExxIntChol, true, 
            "if set true, Exx integrals are Cholesky factorized to 3-index ");

    If.RegisterInputKey("exxdiv_treatment", NULL, &lc.exxdiv_treatment, "gygi-baldereschi",
            CHECK_AND_TERMINATE, OPTIONAL, exxdiv_treatment,
            "Exact exchange method for handling exx divergence at G=0. ",
            "Method not supported. Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("x_gamma_extrapolation", &lc.gamma_extrapolation, true, 
            "if set true, use exx extrapolation to gamma ", XC_OPTIONS);

    If.RegisterInputKey("ExxCholMax", &lc.exxchol_max, 1, 64, 8, 
            CHECK_AND_FIX, OPTIONAL, 
            "maximum number of Exx integral cholesky vectors ", 
            "exxchol_max * num_states ");

    If.RegisterInputKey("exx_fraction", &lc.exx_fraction, -1.0, 1.0, -1.0,
            CHECK_AND_FIX, OPTIONAL,
            "when hybrid functional is used, the fraction of Exx",
            "when the value is negative, it is defined from Functionals ", XC_OPTIONS);

    If.RegisterInputKey("charge_pulay_Gspace", &lc.charge_pulay_Gspace, false, 
            "if set true, charge density mixing the residual in G space ", MIXING_OPTIONS);

    If.RegisterInputKey("drho_precond", &lc.drho_precond, false, 
            "if set true, charge density residual is preconded with q^2/(q^2+q0^2) ", MIXING_OPTIONS);

    If.RegisterInputKey("cube_rho", &lc.cube_rho, true, 
            "if set true, charge density is printed out in cube format ");

    If.RegisterInputKey("cube_vh", &lc.cube_vh, false, 
            "if set true, hatree potential is printed out in cube format ");

    If.RegisterInputKey("cube_pot", &lc.cube_vh, false, 
            "if set true, total potential is printed out in cube format ");

    std::string states_list;

    If.RegisterInputKey("cube_states_list", &states_list, "",
            CHECK_AND_FIX, OPTIONAL,
            "plot the states listed here", "");

    If.RegisterInputKey("exx_int_flag", &lc.exx_int_flag, false, 
            "if set true, calculate the exact exchange integrals ");

    If.RegisterInputKey("noncollinear", &lc.noncoll, false, 
            "if set true, calculate noncollinear ", CONTROL_OPTIONS);

    If.RegisterInputKey("spinorbit", &lc.spinorbit, false, 
            "if set true, calculate with spinorbit coupling ", CONTROL_OPTIONS);

    If.RegisterInputKey("AFM", &lc.AFM, false, 
            "if set true, anti-feromagnetic will be forced by symmetry operation if exist ", CONTROL_OPTIONS);
   

    If.RegisterInputKey("a_length", &celldm[0], 0.0, DBL_MAX, 0.0, 
            CHECK_AND_TERMINATE, OPTIONAL, 
            "First lattice constant. ", 
            "a_length must be a positive number. Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("b_length", &celldm[1], 0.0, DBL_MAX, 0.0, 
            CHECK_AND_TERMINATE, OPTIONAL, 
            "Second lattice constant. ", 
            "b_length must be a positive number. Terminating. ", CONTROL_OPTIONS);

    If.RegisterInputKey("c_length", &celldm[2], 0.0, DBL_MAX, 0.0, 
            CHECK_AND_TERMINATE, OPTIONAL, 
            "Third lattice constant. ", 
            "c_length must be a positive number. Terminating. ", CONTROL_OPTIONS);


    Ri::ReadVector<double> def_lattice_vector({{0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0}});
    Ri::ReadVector<double> lattice_vector;
    If.RegisterInputKey("lattice_vector", &lattice_vector, &def_lattice_vector, 9, OPTIONAL,
            "The simulation cell may be specified using either lattice vectors, "
            "a0, a1, a2 or by lattice constants and a bravais lattice type. "
            "If lattice vectors are used they should be entered as a 3x3 matrix.",
            "Optional Lattice vector input as 3x3 matrix or 9 numbers", CELL_OPTIONS);

    If.RegisterInputKey("grid_spacing", &grid_spacing, 0.0, DBL_MAX, 0.35, 
            CHECK_AND_TERMINATE, OPTIONAL, 
            "Approximate grid spacing (bohr). ", 
            "grid_spacing must be a positive number. Terminating. ", CELL_OPTIONS);

    If.RegisterInputKey("filter_factor", &lc.filter_factor, 0.06, 1.0, 1.0, 
            CHECK_AND_TERMINATE, OPTIONAL, 
            "Filtering factor. ", 
            "filter_factor must lie in the range (0.1, 1.0). Terminating. ", PSEUDO_OPTIONS|EXPERT_OPTION);

    // Default of zero is OK because this means to try to set it automatically later on.
    // The max value of 128 covers any possible hardware scenario I can imagine currently but might
    // need to be adjusted at some point in the future.
    If.RegisterInputKey("omp_threads_per_node", &lc.OMP_THREADS_PER_NODE, 0, 64, 0, 
            CHECK_AND_FIX, OPTIONAL, 
            "Number of Open MP threads each MPI process will use. A value of 0 selects automatic setting. ", 
            "threads_per_node cannnot be a negative number and must be less than 64. ", CONTROL_OPTIONS);

    If.RegisterInputKey("fd_allocation_limit", &lc.fd_allocation_limit, 1024, 262144, 65536, 
            CHECK_AND_FIX, OPTIONAL, 
            "Allocation sizes in finite difference routines less than this value are stack "
            "rather than heap based. ", 
            "fd_allocation_limit must lie in the range 1024 to 262144. ", PERF_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("rmg_threads_per_node", &lc.MG_THREADS_PER_NODE, 0, 64, 0, 
            CHECK_AND_FIX, OPTIONAL, 
            "Number of Multigrid/Davidson threads each MPI process will use. A value of 0 means set automatically.", 
            "threads_per_node cannnot be a negative number and must be less than 64. ", PERF_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("potential_grid_refinement", &lc.FG_RATIO, 0, 4, 0, 
            CHECK_AND_FIX, OPTIONAL, 
            "Ratio of the potential grid density to the wavefunction grid "
            "density. For example if the wavefunction grid is (72,72,72) and "
            "potential_grid_refinement = \"2\" then the potential grid would be "
            "(144,144,144). The default value is 2 but it may sometimes be "
            "beneficial to adjust this. (For USPP the minimum value is also 2 "
            "and it cannot be set lower. NCPP can be set to 1).",
            "potential_grid_refinement must be in the range (0 <= ratio <= 4) where 0 means autoset. ", CELL_OPTIONS);

    If.RegisterInputKey("davidson_multiplier", &lc.davidx, 0, 6, 0, 
            CHECK_AND_FIX, OPTIONAL, 
            "The davidson solver expands the eigenspace with the maximum expansion "
            "factor being set by the value of davidson_multiplier. Larger values "
            "often lead to faster convergence but because the computational cost "
            "of the davidson diagonalization step scales as the cube of the number of "
            "eigenvectors the optimal value based on the fastest time to solution "
            "depends on the number of orbitals. If not specified explicitly or set "
            "to 0 RMG uses the following algorithm to set the value. "
            " & &"
            "Number of orbitals <= 600 davidson_multiplier= \"4\"&"
            "600 < Number of orbitals <= 900    davidson_multiplier = \"3\"&"
            "Number of orbitals > 900           davidson_multiplier = \"2\"&"
            " &"
            "For very large problems the N^3 scaling makes even a factor of 2 "
            "prohibitively costly and the multigrid solver is a better choice. ",
            "davidson_multiplier must be in the range (2 <= davidson_multiplier <= 6).", KS_SOLVER_OPTIONS);

    If.RegisterInputKey("davidson_max_steps", &lc.david_max_steps, 5, 20, 8, 
            CHECK_AND_FIX, OPTIONAL, 
            "Maximum number of iterations for davidson diagonalization.", 
            "davidson_max_steps must be in the range (5 <= davidson_max_steps <= 20). ", KS_SOLVER_OPTIONS);

    If.RegisterInputKey("davidson_premg", &lc.davidson_premg, 0, 8, 4, 
            CHECK_AND_FIX, OPTIONAL, 
            "If the davidson solver is selected this parameter controls the number of multigrid steps to use before enabling davidson.", 
            "davidson_premg must be in the range (0 <= davidson_premg <= 8). ", KS_SOLVER_OPTIONS);

    If.RegisterInputKey("ldaU_radius", &lc.ldaU_radius, 1.0, 12.0, 9.0, 
            CHECK_AND_FIX, OPTIONAL, 
            "Max radius of atomic orbitals to be used in LDA+U projectors. ",
            "ldaU_range must lie in the range (1.0, 12.0). Resetting to the default value of 9.0. ", LDAU_OPTIONS);

    If.RegisterInputKey("potential_acceleration_constant_step", &lc.potential_acceleration_constant_step, 0.0, 4.0, 0.0, 
            CHECK_AND_FIX, OPTIONAL, 
            "When set to a non-zero value this parameter causes RMG to "
            "perform a band by band update of the self-consistent potential "
            "during the course of an SCF step when the multigrid kohn_sham_solver "
            "is chosen. This means that updates to the lower energy orbitals "
            "are incorporated into the SCF potential seen by the higher energy orbitals "
            "as soon as they are computed. This can lead to faster convergence "
            "and better stability for many systems. The option should only be used "
            "with Linear mixing. Even when the davidson solver is chosen this parameter "
            "may be used since the first few steps with davidson usually uses the "
            "multigrid solver.",
            "potential_acceleration_constant_step must lie in the range (0.0, 4.0). Resetting to the default value of 0.0. ", MIXING_OPTIONS);

    If.RegisterInputKey("tddft_time_step", &lc.tddft_time_step, 0.0, DBL_MAX, 0.2, 
            CHECK_AND_TERMINATE, OPTIONAL,
            "TDDFT time step for use in TDDFT mode ",
            "tddft_time_step is in atomic unit ", MD_OPTIONS);

    If.RegisterInputKey("ionic_time_step", &lc.iondt, 0.0, DBL_MAX, 50.0, 
            CHECK_AND_TERMINATE, OPTIONAL,
            "Ionic time step for use in molecular dynamics and structure optimizations. ",
            "ionic_time_step must be greater than 0.0. ", MD_OPTIONS);

    If.RegisterInputKey("ionic_time_step_increase", &lc.iondt_inc, 1.0, 3.0, 1.1,
            CHECK_AND_FIX, OPTIONAL,
            "Factor by which ionic timestep is increased when dynamic timesteps are enabled. ",
            "ionic_time_step_increase must lie in the range (1.0,1.1). Resetting to the default value of 1.1. ", MD_OPTIONS);

    If.RegisterInputKey("ionic_time_step_decrease", &lc.iondt_dec, 0.0, 1.0, 0.5,
            CHECK_AND_FIX, OPTIONAL,
            "Factor by which ionic timestep is decreased when dynamic timesteps are enabled. ",
            "ionic_time_step_decrease must lie in the range (0.0,1.0). Resetting to the default value of 0.5. ", MD_OPTIONS);

    If.RegisterInputKey("max_ionic_time_step", &lc.iondt_max, 0.0, 150.0, 150.0,
            CHECK_AND_FIX, OPTIONAL,
            "Maximum ionic time step to use for molecular dynamics or structural optimizations. ",
            "max_ionic_time_step must lie in the range (0.0,150.0). Resetting to the default value of 150.0. ", MD_OPTIONS);

    If.RegisterInputKey("qmc_nband", &lc.qmc_nband, 0, INT_MAX, 0, 
            CHECK_AND_FIX, OPTIONAL, 
            "The number of band used in rmg-qmcpack interface. ", 
            "By default without this input, it will be set to number of states. ", CONTROL_OPTIONS);

    If.RegisterInputKey("unoccupied_states_per_kpoint", &lc.num_unocc_states, 0, INT_MAX, 10, 
            CHECK_AND_FIX, OPTIONAL, 
            "The number of unoccupied orbitals. A value that is 15-20% of the number of occupied orbitals generally works well.", 
            "Unoccupied_states_per_kpoint must be greater than 0. Setting to default value of 10. ", OCCUPATION_OPTIONS);

    If.RegisterInputKey("state_block_size", &lc.state_block_size, 1, INT_MAX, 64, 
            CHECK_AND_FIX, OPTIONAL, 
            "state_block used in nlforce. ", 
            "it is better to be 2^n. ", PERF_OPTIONS|EXPERT_OPTION);


    If.RegisterInputKey("extra_random_lcao_states", &lc.extra_random_lcao_states, 0, INT_MAX, 0, 
            CHECK_AND_TERMINATE, OPTIONAL, 
            "LCAO (Linear Combination of Atomic Orbitals) is the default startup method "
            "for RMG. The atomic orbitals are obtained from the pseudpotentials but in some "
            "cases better convergence may be obtained by adding extra random wavefunctions "
            "in addition to the atomic orbitals.", 
            "extra_random_lcao_states must be greater than 0. Terminating. ", DIAG_OPTIONS);

    If.RegisterInputKey("system_charge", &lc.background_charge, -DBL_MAX, DBL_MAX, 0.0,
            CHECK_AND_FIX, OPTIONAL,
            "Number of excess holes in the system (useful for doped systems). Example: 2 means system is missing two electrons ",
            "system_charge must be a real number. ", CONTROL_OPTIONS);

    If.RegisterInputKey("occupation_electron_temperature_eV", &lc.occ_width, 0.0, 2.0, 0.04,
            CHECK_AND_FIX, OPTIONAL,
            "Target electron temperature when not using fixed occupations.  ",
            "occupation_electron_temperature_eV must lie in the range (0.0,2.0). Resetting to the default value of 0.04. ", OCCUPATION_OPTIONS);

    std::string stm_bias_str;
    If.RegisterInputKey("STM_bias", &stm_bias_str, "-1.0 1.0",
            CHECK_AND_FIX, OPTIONAL,
            "bias (in unit of Volt)  for STM calculation  ",
            "integrate density in energy window EF->EF+bias ", CONTROL_OPTIONS);
    std::string stm_height_str;
    If.RegisterInputKey("STM_height", &stm_height_str, "2.0 4.0",
            CHECK_AND_FIX, OPTIONAL,
            "height range for STM calculation  ",
            "vaccum must be along z direction, height is above the maximum z coordinate ", CONTROL_OPTIONS);

    If.RegisterInputKey("occupation_number_mixing", &lc.occ_mix, 0.0, 1.0, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "Mixing parameter for orbital occupations when not using fixed occupations. ",
            "occupation_number_mixing must lie in the range (0.0,1.0). Resetting to the default value of 0.3. ", OCCUPATION_OPTIONS);

    If.RegisterInputKey("MP_order", &lc.mp_order, 0, 5, 2, 
            CHECK_AND_FIX, OPTIONAL, 
            "order of Methefessel Paxton occupation.", 
            "0 means simple error function as distribution ", OCCUPATION_OPTIONS);

    If.RegisterInputKey("period_of_diagonalization", &lc.diag, 0, INT_MAX, 1, 
            CHECK_AND_FIX, OPTIONAL, 
            "Diagonalization period (per scf step). Mainly for debugging and should not be changed for production.", 
            "Diagonalization period must be greater than 0. Resetting to the default value of 1. ", DIAG_OPTIONS);

    If.RegisterInputKey("max_scf_steps", &lc.max_scf_steps, 0, INT_MAX, 500,
            CHECK_AND_FIX, OPTIONAL, 
            "Maximum number of self consistent steps to perform. Inner loop for hybrid functionals. ", 
            "max_scf_steps must be greater than 0. Resetting to the default value of 500 ", CONTROL_OPTIONS);

    If.RegisterInputKey("max_exx_steps", &lc.max_exx_steps, 1, INT_MAX, 100,
            CHECK_AND_FIX, OPTIONAL, 
            "Maximum number of self consistent steps to perform with hybrid functionals. ", 
            "max_exx_steps must be greater than 0. Resetting to the default value of 500 ", CONTROL_OPTIONS);

    If.RegisterInputKey("tddft_steps", &lc.tddft_steps, 0, INT_MAX, 2000,
            CHECK_AND_FIX, OPTIONAL, 
            "Maximum number of tddft steps to perform. ", 
            "tddft steps must be greater than 0. Resetting to the default value of 2000 ", TDDFT_OPTIONS);

    If.RegisterInputKey("charge_pulay_order", &lc.charge_pulay_order, 1, 10, 5,
            CHECK_AND_FIX, OPTIONAL,
            "Number of previous steps to use when Pulay mixing is used to update the charge density.",
            "", MIXING_OPTIONS);

    If.RegisterInputKey("charge_pulay_scale", &lc.charge_pulay_scale, 0.0, 1.0, 0.50,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "charge_pulay_scale must lie in the range (0.0,1.0). Resetting to the default value of 0.50 ", MIXING_OPTIONS);


    If.RegisterInputKey("ldau_pulay_order", &lc.ldau_pulay_order, 1, 10, 5,
            CHECK_AND_FIX, OPTIONAL,
            "Number of previous steps to use when Pulay mixing is used to update the ldau occupation .",
            "", MIXING_OPTIONS);

    If.RegisterInputKey("ldau_pulay_scale", &lc.ldau_pulay_scale, 0.0, 1.0, 1.00,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "ldau_pulay_scale must lie in the range (0.0,1.0). Resetting to the default value of 1.00 ", MIXING_OPTIONS);

    If.RegisterInputKey("unoccupied_tol_factor", &lc.unoccupied_tol_factor, 1.0, 100000.0, 1000.0,
            CHECK_AND_FIX, OPTIONAL,
            "When using the Davidson Kohn-Sham solver unoccupied states are converged to a less stringent tolerance than occupied orbitals with the ratio set by this parameter.",
            "unoccupied_tol_factor must lie in the range (0.000001,100000.0). Resetting to the default value of 1000.0 ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("charge_pulay_refresh", &lc.charge_pulay_refresh, 1, INT_MAX, 100,
            CHECK_AND_FIX, OPTIONAL,
            "charge Pulay mixing reset steps.",
            "", MIXING_OPTIONS);
    If.RegisterInputKey("ldau_pulay_refresh", &lc.ldau_pulay_refresh, 1, INT_MAX, 100,
            CHECK_AND_FIX, OPTIONAL,
            "ldau pulay mixing reset steps",
            "", MIXING_OPTIONS);

    If.RegisterInputKey("charge_broyden_order", &lc.charge_broyden_order, 1, 10, 5,
            CHECK_AND_FIX, OPTIONAL,
            "Number of previous steps to use when Broyden mixing is used to update the charge density.",
            "", MIXING_OPTIONS);

    If.RegisterInputKey("charge_broyden_scale", &lc.charge_broyden_scale, 0.0, 1.0, 0.50,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "charge_broyden_scale must lie in the range (0.0,1.0). Resetting to the default value of 0.50 ", MIXING_OPTIONS);

    If.RegisterInputKey("projector_expansion_factor", &lc.projector_expansion_factor, 0.5, 3.0, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "When using localized projectors the radius can be adjusted with this parameter.",
            "projector_expansion_factor must lie in the range (0.5,3.0). Resetting to the default value of 1.0 ", PSEUDO_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("write_data_period", &lc.checkpoint, 5, 50, 5,
            CHECK_AND_FIX, OPTIONAL,
            "How often to write checkpoint files during the initial quench in units of SCF steps. During structural relaxations of molecular dynamics checkpoints are written each ionic step.",
            "", CONTROL_OPTIONS);

    If.RegisterInputKey("write_eigvals_period", &lc.write_eigvals_period, 1, 100, 5,
            CHECK_AND_FIX, OPTIONAL,
            "How often to output eigenvalues in units of scf steps.",
            "write_eigvals_period must lie in the range (1,100). Resetting to the default value of 5. ", OUTPUT_OPTIONS);

    If.RegisterInputKey("max_md_steps", &lc.max_md_steps, 0, INT_MAX, 100,
            CHECK_AND_TERMINATE, OPTIONAL,
            "Maximum number of molecular dynamics steps to perform.",
            "max_md_steps must be a positive value. Terminating. ", MD_OPTIONS);

    If.RegisterInputKey("hartree_max_sweeps", &lc.hartree_max_sweeps, 5, 100, 10,
            CHECK_AND_FIX, OPTIONAL,
            "Maximum number of hartree iterations to perform per scf step. ",
            "hartree_max_sweeps must lie in the range (5,100). Resetting to the default value of 10. ", POISSON_OPTIONS);

    If.RegisterInputKey("hartree_min_sweeps", &lc.hartree_min_sweeps, 0, 5, 5,
            CHECK_AND_FIX, OPTIONAL,
            "Minimum number of hartree iterations to perform per scf step. ",
            "hartree_min_sweeps must lie in the range (0.5). Resetting to the default value of 5. ", POISSON_OPTIONS);

    If.RegisterInputKey("kohn_sham_pre_smoothing", &lc.eig_parm.gl_pre, 1, 5, 2,
            CHECK_AND_FIX, OPTIONAL,
            "Number of global grid pre-smoothing steps to perform before a "
            "multigrid preconditioner iteration. ",
            "kohn_sham_pre_smoothing must lie in the range (1,5). Resetting to the default value of 2. ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("kohn_sham_post_smoothing", &lc.eig_parm.gl_pst, 1, 5, 2,
            CHECK_AND_FIX, OPTIONAL,
            "Number of global grid post-smoothing steps to perform after a "
            "multigrid preconditioner iteration. ",
            "kohn_sham_post_smoothing must lie in the range (1,5). Resetting to the default value of 2. ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("kohn_sham_mucycles", &lc.eig_parm.mucycles, 1, 6, 2,
            CHECK_AND_FIX, OPTIONAL,
            "Number of mu (also known as W) cycles to use in the kohn-sham multigrid preconditioner. ",
            "kohn_sham_mucycles must lie in the range (1,6). Resetting to the default value of 2. ", KS_SOLVER_OPTIONS);

    If.RegisterInputKey("kohn_sham_fd_order", &lc.kohn_sham_fd_order, 6, 10, 8,
            CHECK_AND_FIX, OPTIONAL,
            "RMG uses finite differencing to represent the kinetic energy operator "
            "and the accuracy of the representation is controllable by the "
            "kohn_sham_fd_order parameter. The default is 8 and is fine for most "
            "purposes but higher accuracy is obtainable with 10th order at the cost "
            "of some additional computational expense.",
            "kohn_sham_fd_order must lie in the range (6,10). Resetting to the default value of 8. ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("use_gpu_fd", &lc.use_gpu_fd, false, 
            "Use gpus for kohn-sham orbital finite differencing. Depending on the "
            "balance of hardware characteristics this can provide a significant "
            "speedup but individual testing is required. Experimental.", MISC_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("force_grad_order", &lc.force_grad_order, 0, 12, 8,
            CHECK_AND_FIX, OPTIONAL,
            "Atomic forces may be computed to varying degrees of accuracy depending "
            "on the requirements of a specific problem. A value of 0 implies highest "
            "accuracy which is obtained by using FFTs in place of finite differencing. ",
            "kohn_sham_fd_order must lie in the range (4,12). Resetting to the default value of 8. ", MD_OPTIONS);

    If.RegisterInputKey("kohn_sham_coarse_time_step", &lc.eig_parm.sb_step, 0.0, 1.2, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "Time step to use in the kohn-sham multigrid solver on the coarse levels. ",
            "kohn_sham_coarse_time_step must lie in the range (0.5,1.2). Resetting to the default value of 1.0. ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("kohn_sham_time_step", &lc.eig_parm.gl_step, 0.0, 2.0, 0.66,
            CHECK_AND_FIX, OPTIONAL,
            "Smoothing timestep to use on the fine grid in the the kohn-sham multigrid preconditioner. ",
            "kohn_sham_time_step must lie in the range (0.4,2.0). Resetting to the default value of 0.66. ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("kohn_sham_mg_timestep", &lc.eig_parm.mg_timestep, 0.0, 2.0, 0.6666666666666,
            CHECK_AND_FIX, OPTIONAL,
            "timestep for multigrid correction. ",
            "kohn_sham_mg_step must lie in the range (0.0,2.0). Resetting to the default value of 0.66 ", KS_SOLVER_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("poisson_pre_smoothing", &lc.poi_parm.gl_pre, 1, 6, 2,
            CHECK_AND_FIX, OPTIONAL,
            "Number of global hartree grid pre-smoothing steps to perform before a multigrid iteration. ",
            "poisson_pre_smoothing must lie in the range (1,6). Resetting to the default value of 2. ", POISSON_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("poisson_post_smoothing", &lc.poi_parm.gl_pst, 1, 6, 1,
            CHECK_AND_FIX, OPTIONAL,
            "Number of global hartree grid post-smoothing steps to perform after a multigrid iteration. ",
            "", POISSON_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("poisson_mucycles", &lc.poi_parm.mucycles, 1, 4, 3,
            CHECK_AND_FIX, OPTIONAL,
            "Number of mu (also known as W) cycles to use in the hartree multigrid solver. ",
            "poisson_mucycles must lie in the range (1,4). Resetting to the default value of 3. ", POISSON_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("poisson_finest_time_step", &lc.poi_parm.gl_step, 0.4, 1.0, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "Time step to use in the poisson multigrid solver on the finest level. ",
            "poisson_finest_time_step must lie in the range (0.4,1.0). Resetting to the default value of 1.0. ", POISSON_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("poisson_coarse_time_step", &lc.poi_parm.sb_step, 0.4, 1.0, 0.8,
            CHECK_AND_FIX, OPTIONAL,
            "Time step to use in the poisson multigrid solver on the coarse levels. ",
            "poisson_coarse_time_step must lie in the range (0.4,1.0). Resetting to the default value of 0.8. ", POISSON_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("poisson_coarsest_steps", &lc.poi_parm.coarsest_steps, 10, 100, 25,
            CHECK_AND_FIX, OPTIONAL,
            "Number of smoothing steps to use on the coarsest level in the hartree multigrid solver. ",
            "poisson_coarsest_steps must lie in the range (10,100). Resetting to the default value of 25. ", POISSON_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("kohn_sham_mg_levels", &lc.eig_parm.levels, -1, 6, -1,
            CHECK_AND_FIX, OPTIONAL,
            "Number of multigrid levels to use in the kohn-sham multigrid preconditioner. ",
            "kohn_sham_mg_levels must lie in the range (-1,6) where -1=automatic. Resetting to the default value of automatic (-1). ", KS_SOLVER_OPTIONS);

    If.RegisterInputKey("poisson_mg_levels", &lc.poi_parm.levels, -1, 6, -1,
            CHECK_AND_FIX, OPTIONAL,
            "Number of multigrid levels to use in the hartree multigrid solver. ",
            "poisson_mg_levels must lie in the range (-1,6) where -1=automatic. Resetting to the default value of automatic (-1). ",
            POISSON_OPTIONS);

    If.RegisterInputKey("scalapack_block_factor", &lc.scalapack_block_factor, 4, 512,32,
            CHECK_AND_FIX, OPTIONAL,
            "Block size to use with scalapack. Optimal value is dependent on matrix size and system hardware. ",
            "scalapack_block_factor must lie in the range (4,512). Resetting to the default value of 32. ", DIAG_OPTIONS);

    If.RegisterInputKey("non_local_block_size", &lc.non_local_block_size, 64, 40000, 512,
            CHECK_AND_FIX, OPTIONAL,
            "Block size to use when applying the non-local and S operators. ",
            "non_local_block_size must lie in the range (64,40000). Resetting to the default value of 512. ", PERF_OPTIONS);

    If.RegisterInputKey("E_POINTS", &lc.E_POINTS, 201, 201, 201,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "");

    If.RegisterInputKey("md_number_of_nose_thermostats", &lc.nose.m, 5, 5, 5,
            CHECK_AND_FIX, OPTIONAL,
            "Number of Nose thermostats to use during Constant Volume and Temperature MD.",
            "", MD_OPTIONS);

    If.RegisterInputKey("dynamic_time_delay", &lc.relax_steps_delay, 5, 5, 5,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "", MD_OPTIONS);

    If.RegisterInputKey("dynamic_time_counter", &lc.relax_steps_counter, 0, 0 , 0,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "", MD_OPTIONS);

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

    If.RegisterInputKey("coalesce_factor", &pelc.coalesce_factor, 1, 16, 4,
            CHECK_AND_FIX, OPTIONAL,
            "Grid coalescing factor.",
            "coalesce_factor must lie in the range (1,8). Resetting to default value of 4.", CONTROL_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("charge_density_mixing", &lc.mix, 0.0, 1.0, 0.5,
            CHECK_AND_FIX, OPTIONAL,
            "Proportion of the current charge density to replace with the new density after each scf step when linear mixing is used. ",
            "charge_density_mixing must lie in the range (0.0, 1.0) Resetting to the default value of 0.5. ", MIXING_OPTIONS);

    If.RegisterInputKey("ldau_mixing", &lc.ldau_mix, 0.0, 1.0, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "Proportion of the current ldau occupation to replace with the new ones after each scf step when linear mixing is used. ",
            "charge_density_mixing must lie in the range (0.0, 1.0) Resetting to the default value of 1.0. ", MIXING_OPTIONS);

    If.RegisterInputKey("drho_precond_q0", &lc.drho_q0, 0.0, 10.0, 0.5,
            CHECK_AND_FIX, OPTIONAL,
            "preconding the charge density residual by q^2/(q^2+q0^2) ",
            "See Kresse and Furthmueller,  Computational Materials Science 6 (1996) 15-50  ", MIXING_OPTIONS);
    If.RegisterInputKey("folded_spectrum_width", &lc.folded_spectrum_width, 0.10, 1.0, 0.3,
            CHECK_AND_FIX, OPTIONAL,
            "Submatrix width to use as a fraction of the full spectrum. "
            "The folded spectrum width ranges from 0.10 to 1.0. For insulators and "
            "semiconductors a value of 0.3 is appropriate. For metals values between "
            "0.15 to 0.2 tend to be better. The default value is 0.3 ",
            "folded_spectrum_width must lie in the range (0.10,1.0). Resetting to the default value of 0.3. ", DIAG_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("folded_spectrum_iterations", &lc.folded_spectrum_iterations, 0, 20, 2,
            CHECK_AND_FIX, OPTIONAL,
            "Number of folded spectrum iterations to perform. ",
            "folded_spectrum_iterations must lie in the range (0,20). Resetting to the default value of 2. ", DIAG_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("laplacian_offdiag", &lc.laplacian_offdiag, false, 
            "if set to true, we use LaplacianCoeff.cpp to generate coeff", MISC_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("laplacian_autocoeff", &lc.laplacian_autocoeff, false, 
            "if set to true, we use LaplacianCoeff.cpp to generate coeff", MISC_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("use_cpdgemr2d", &lc.use_cpdgemr2d, true, 
            "if set to true, we use Cpdgemr2d to change matrix distribution");

    //RMG2BGW options
    If.RegisterInputKey("use_symmetry", &lc.is_use_symmetry, true, 
            "For non-gamma point, always true, for gamma point, optional", CELL_OPTIONS);

    If.RegisterInputKey("frac_symmetry", &lc.frac_symm, true, 
            "For supercell calculation, one can disable the fractional translation symmetry", CELL_OPTIONS);

    If.RegisterInputKey("rmg2bgw", &lc.rmg2bgw, false, 
            "Write wavefunction in G-space to BerkeleyGW WFN file.", MISC_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("ecutrho", &lc.ecutrho, 0.0, 10000.0, 0.0,
            CHECK_AND_FIX, OPTIONAL,
            "ecut for rho in unit of Ry. ",
            " ");

    If.RegisterInputKey("ecutwfc", &lc.ecutwfc, 0.0, 10000.0, 0.0,
            CHECK_AND_FIX, OPTIONAL,
            "ecut for wavefunctions in unit of Ry. ",
            " ");

    If.RegisterInputKey("vxc_diag_nmin", &lc.vxc_diag_nmin, 1, 10000, 1,
            CHECK_AND_FIX, OPTIONAL,
            "Minimum band index for diagonal Vxc matrix elements. ",
            "vxc_diag_nmin must lie in the range (1, 10000). Resetting to the default value of 1. ");

    If.RegisterInputKey("vxc_diag_nmax", &lc.vxc_diag_nmax, 1, 10000, 1,
            CHECK_AND_FIX, OPTIONAL,
            "Maximum band index for diagonal Vxc matrix elements. ",
            "vxc_diag_nmax must lie in the range (1, 10000). Resetting to the default value of 1. ");

    If.RegisterInputKey("max_nlradius", &lc.max_nlradius, 2.0, 10000.0, 10000.0,
            CHECK_AND_FIX, OPTIONAL,
            "maximum radius for non-local projectors ",
            " ", PSEUDO_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("min_nlradius", &lc.min_nlradius, 1.0, 10000.0, 2.0,
            CHECK_AND_FIX, OPTIONAL,
            "minimum radius for non-local projectors ",
            " ", PSEUDO_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("max_qradius", &lc.max_qradius, 2.0, 10000.0, 10000.0,
            CHECK_AND_FIX, OPTIONAL,
            "maximum radius for qfunc in ultra-pseudopotential ",
            " ", PSEUDO_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("min_qradius", &lc.min_qradius, 1.0, 10000.0, 2.0,
            CHECK_AND_FIX, OPTIONAL,
            "minimum radius for qfunc in ultra-pseudopotential ",
            " ", PSEUDO_OPTIONS|EXPERT_OPTION);

    // Testing options
    If.RegisterInputKey("test_energy", &lc.test_energy, -1.0e9, 1.0e9, (double)NAN,
            CHECK_AND_FIX, OPTIONAL,
            "Expected final energy for testing. ",
            "test_energy must lie in the range (-1.0e9,1.0e9). Ignoring. ", TESTING_OPTIONS);

    If.RegisterInputKey("test_energy_tolerance", &lc.test_energy_tolerance, 1.0e-8, 1.0e-4, 1.0e-7,
            CHECK_AND_FIX, OPTIONAL,
            "Test final energy tolerance. ",
            "test_energy_tolerance must lie in the range (1.0e-8,1.0e-4). Resetting to the default value of 1.0e-7. ", TESTING_OPTIONS);

    If.RegisterInputKey("test_bond_length", &lc.test_bond_length, 0.0 , 20.0, (double)NAN,
            CHECK_AND_FIX, OPTIONAL,
            "Expected dimer bond length for testing. ",
            "bond length must lie in the range (0.0, 20.0). Ignoring. ", TESTING_OPTIONS);

    If.RegisterInputKey("test_bond_length_tolerance", &lc.test_bond_length_tolerance, 1.0e-4, 1.0e-1, 1.0e-3,
            CHECK_AND_FIX, OPTIONAL,
            "Test bond length tolerance. ",
            "test_bond_length_tolerance must lie in the range (1.0e-4,1.0e-1). Resetting to the default value of 1.0e-3. ", TESTING_OPTIONS);


    // Booleans next. Booleans are never required.
#if CUDA_ENABLED
    // If GPU memory is constrained this one should be set to true.
    If.RegisterInputKey("pin_nonlocal_weights", &lc.pin_nonlocal_weights, false,
            "Flag indicating whether or not nonlocal weights should use pinned instead of managed memory.", PERF_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("use_cublasxt", &lc.use_cublasxt, false,
            "This flag enables the use of the cublasxt library in place of the standard cublas library with cuda enabled builds. Intended for use when GPU memory is constrained.", CONTROL_OPTIONS);

#endif

    If.RegisterInputKey("use_bessel_projectors", &lc.use_bessel_projectors, false,
            "When a semi-local pseudopotential is being used projectors will be generated using Bloechl's procedure with Bessel functions as the basis set if this is true.", PSEUDO_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("write_orbital_overlaps", &lc.write_orbital_overlaps, false,
            "If true the orbital overlap matrix from successive MD steps is written.", OUTPUT_OPTIONS);

    If.RegisterInputKey("kohn_sham_ke_fft", &lc.kohn_sham_ke_fft, false,
            "Special purpose flag which will force use of an FFT for the kinetic energy operator.", MISC_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("fast_density", &lc.fast_density, true,
            "Use a faster but less accurate method to generate the charge density from the electronic wavefunctions. "
            "As the cutoff (grid-density) increases this method improves in accuracy. This option should be set to "
            "false if you receive warnings about negative charge densities after interpolation.");

    If.RegisterInputKey("lcao_use_empty_orbitals", &lc.lcao_use_empty_orbitals, false,
            "Some pseudopotentials contain unbound atomic orbitals and this flag indicates "
            "whether or not they should be used for LCAO starts.", MISC_OPTIONS);

    If.RegisterInputKey("write_serial_restart", &lc.write_serial_restart, false,
            "RMG normally writes parallel restart files. These require that restarts have the "
            "same processor topology. If write_serial_restart = \"true\" then RMG will also "
            "write a serial restart file that can be used with a different processor topology ", CONTROL_OPTIONS);

    If.RegisterInputKey("read_serial_restart", &lc.read_serial_restart, false,
            "Directs RMG to read from serial restart files. Normally used when changing "
            "the sprocessor topology used during a restart run ", CONTROL_OPTIONS);

    If.RegisterInputKey("write_qmcpack_restart", &lc.write_qmcpack_restart, false,
            "If true then a QMCPACK restart file is written as well as a serial restart file.", CONTROL_OPTIONS);
    If.RegisterInputKey("write_qmcpack_restart_localized", &lc.write_qmcpack_restart_localized, false,
            "If true then a QMCPACK restart file for localized orbitals", CONTROL_OPTIONS);

    If.RegisterInputKey("compressed_infile", &lc.compressed_infile, true,
            "Flag indicating whether or not parallel restart wavefunction file uses compressed format.", CONTROL_OPTIONS);

    If.RegisterInputKey("compressed_outfile", &lc.compressed_outfile, true,
            "Flag indicating whether or not  parallel output wavefunction file uses compressed format.", CONTROL_OPTIONS);

    If.RegisterInputKey("nvme_weights", &lc.nvme_weights, false,
            "Flag indicating whether or not projector weights should be mapped to disk.", CONTROL_OPTIONS);

    If.RegisterInputKey("nvme_work", &lc.nvme_work, false,
            "Flag indicating whether or not work arrays should be mapped to disk.", CONTROL_OPTIONS);

    If.RegisterInputKey("nvme_orbitals", &lc.nvme_orbitals, false,
            "Flag indicating whether or not orbitals should be mapped to disk.", CONTROL_OPTIONS);

    If.RegisterInputKey("alt_laplacian", &lc.alt_laplacian, true,
            "Flag indicating whether or not to use alternate laplacian weights for some operators.", MISC_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("use_alt_zgemm", &lc.use_alt_zgemm, false,
            "Flag indicating whether or not to use alternate zgemm implementation.", PERF_OPTIONS);

    If.RegisterInputKey("filter_dpot", &lc.filter_dpot, false,
            "Flag indicating whether or not to filter density dependent potentials.", PSEUDO_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("sqrt_interpolation", &lc.sqrt_interpolation, false,
            "Flag indicating whether or not to use square root technique for density interpolation.", MISC_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("renormalize_forces", &lc.renormalize_forces, true,
            "Flag indicating whether or not to renormalize forces.", MD_OPTIONS);

    If.RegisterInputKey("coalesce_states", &lc.coalesce_states, false,
            "Flag indicating whether or not to coalesce states.", CONTROL_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("localize_projectors", &lc.localize_projectors, true,
            "The Beta function projectors for a particular ion decay rapidly "
            "in real-space with increasing r. For large cells truncating the "
            "real-space representation of the projector can lead to "
            "significant computational savings with a small loss of accuracy. "
            "For smaller cells the computational cost is the same for localized "
            "or delocalized projectors so it is better to set localize_projectors "
            "to false.", PSEUDO_OPTIONS);

    If.RegisterInputKey("localize_localpp", &lc.localize_localpp, true,
            "The local potential associated with a particular ion also decays "
            "rapidly in real-space with increasing r. As with beta projectors "
            "truncating the real-space representation for large cells can lead "
            "to significant computational savings with a small loss of accuracy "
            "but it should be set to false for small cells.", PSEUDO_OPTIONS);

    If.RegisterInputKey("write_pseudopotential_plots", &lc.write_pp_flag, false,
            "Flag to indicate whether or not to write pseudopotential plots. ", OUTPUT_OPTIONS);

    If.RegisterInputKey("equal_initial_density", &lc.init_equal_density_flag, false,
            "Specifies whether to set initial up and down density to be equal.");

    If.RegisterInputKey("write_pdos", &lc.pdos_flag, false,
            "Flag to write partial density of states.", OUTPUT_OPTIONS);

    If.RegisterInputKey("initial_diagonalization", &lc.initdiag, true, 
            "Perform initial subspace diagonalization.", DIAG_OPTIONS);

    If.RegisterInputKey("verbose", &lc.verbose, false,
            "Flag for writing out extra information ", CONTROL_OPTIONS);

    If.RegisterInputKey("folded_spectrum", &lc.use_folded_spectrum, false, 
            "When the number of eigenvectors is large using folded_spectrum is "
            "substantially faster than standard diagonalization. It also tends "
            "to converge better for metallic systems. It works with the "
            "multigrid kohn_sham_solver but not the davidson solver. ", DIAG_OPTIONS|EXPERT_OPTION);


    If.RegisterInputKey("use_numa", &lc.use_numa, true, 
            "Numa stands for Non Uniform Memory Access and means that the main "
            "memory of a computer is organized into seperate distinct banks. "
            "Each bank is then attached to a CPU core or group of cores and "
            "while all cores can normally access all banks the access speed "
            "is faster for directly attached banks. Ensuring that individual "
            "CPU cores mostly access memory in banks they are directly attached "
            "to can have a large impact on performance. Process mapping that "
            "does this can normally be done when jobs are submitted and run "
            "via arguments to mpirun/mpiexec but if this is not done RMG will "
            "attempt to provide an optimal mapping if use_numa is set to true. ",
            PERF_OPTIONS);

    If.RegisterInputKey("use_hwloc", &lc.use_hwloc, false, 
            "Use internal hwloc setup if available. If both this and use_numa are true hwloc takes precedence.", PERF_OPTIONS);

    If.RegisterInputKey("use_async_allreduce", &lc.use_async_allreduce, true, 
            "RMG uses MPI_Allreduce function calls in several places and for "
            "large problems these can account for a significant fraction of the "
            "total run time. In most cases using the asynchronous MPI versions "
            "of the functions is faster but this is not true for all platforms "
            "and in that casesetting this flag to false can improve performance.", PERF_OPTIONS);

    If.RegisterInputKey("mpi_queue_mode", &lc.mpi_queue_mode, true, 
            "Use mpi queue mode.", PERF_OPTIONS);

    If.RegisterInputKey("spin_manager_thread", &lc.spin_manager_thread, true, 
            "When mpi_queue_mode is enabled the manager thread spins instead of sleeping.", PERF_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("spin_worker_threads", &lc.spin_worker_threads, true, 
            "When mpi_queue_mode is enabled the worker threads spin instead of sleeping.", PERF_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("require_huge_pages", &lc.require_huge_pages, false, 
            "If set RMG assumes that sufficient huge pages are available. "
            "Bad things may happen if this is not true.", PERF_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("relax_dynamic_timestep", NULL, false,
            "Flag indicating whether or not to use dynamic timesteps in relaxation mode. ", MD_OPTIONS);

    If.RegisterInputKey("freeze_occupied", NULL, false,
            "Flag indicating whether or not to freeze the density and occupied orbitals after a restart. ", MISC_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("relax_max_force", &lc.thr_frc, 0.0, DBL_MAX, 2.5E-3,
            CHECK_AND_FIX, OPTIONAL,
            "Force value at which an ionic relaxation is considered to be converged. ",
            "relax_max_force must be a positive value. Resetting to default value of 2.5e-03. ", MD_OPTIONS);

    If.RegisterInputKey("md_randomize_velocity", &lc.nose.randomvel, true,
            "The initial ionic velocities for a molecular dyanamics run are randomly initialized to the target temperature.", MD_OPTIONS);

    If.RegisterInputKey("output_rho_xsf", NULL, false,
            "Generate xsf format for electronic density.", OUTPUT_OPTIONS);

    If.RegisterInputKey("rms_convergence_criterion", &lc.thr_rms, 0.0, 1.0e-3, 1.0e-7,
            CHECK_AND_FIX, OPTIONAL,
            "The RMS value of the change in the total potential from step to step "
            "where we assume self consistency has been achieved.",
            "rms_convergence_criterion must lie in the range (1.0e-04,1.0e-14). Resetting to default value of 1.0e-7. ", CONTROL_OPTIONS);

    If.RegisterInputKey("stress_convergence_criterion", &lc.thr_stress, 0.0, 5.0e+1, 0.5,
            CHECK_AND_FIX, OPTIONAL,
            "The stress criteria ",
            "stress_convergence_criterion must lie in the range (1.0e-04,50). Resetting to default value of 0.5. ", CONTROL_OPTIONS);

    If.RegisterInputKey("energy_convergence_criterion", &lc.thr_energy, 1.0e-20, 1.0e-7, 1.0e-10,
            CHECK_AND_FIX, OPTIONAL,
            "The RMS value of the estimated change in the total energy per step where we assume self "
            "consistency has been achieved. ",
            "rms_convergence_criterion must lie in the range (1.0e-07,1.0e-20). Resetting to default value of 1.0e-10. ", CONTROL_OPTIONS);

    If.RegisterInputKey("exx_convergence_criterion", &lc.exx_convergence_criterion, 1.0e-12, 1.0e-6, 1.0e-9,
            CHECK_AND_FIX, OPTIONAL,
            "Convergence criterion for the EXX delta from step to step where we assume EXX "
            "consistency has been achieved. ",
            "exx_convergence_criterion must lie in the range (1.0e-12,1.0e-6). Resetting to default value of 1.0e-9. ", XC_OPTIONS);

    If.RegisterInputKey("vexx_fft_threshold", &lc.vexx_fft_threshold, 1.0e-14, 1.0e-1, 1.0e-14,
            CHECK_AND_FIX, OPTIONAL,
            "The value for the EXX delta where we switch from single to double precision ffts. "
            "Single precision ffts are generally accurate enough. ",
            "vexx_fft_threshold must lie in the range (1.0e-14,1.0e-1). Resetting to default value of 1.0e-14. ", XC_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("preconditioner_threshold", &lc.preconditioner_thr, 1.0e-9, 1.0e-1, 1.0e-1,
            CHECK_AND_FIX, OPTIONAL,
            "The RMS value of the change in the total potential where we switch "
            "the preconditioner from single to double precision.",
            "preconditioner_threshold must lie in the range (1.0e-9,1.0e-1). Resetting to default value of 1.0e-1. ", PERF_OPTIONS|EXPERT_OPTION);

    If.RegisterInputKey("gw_residual_convergence_criterion", &lc.gw_threshold, 1.0e-14, 4.0e-4, 1.0e-6,
            CHECK_AND_FIX, OPTIONAL,
            "The max value of the residual for unoccupied orbitals when performing a GW calculation. ",
            "gw_residual_convergence_criterion must lie in the range (4.0e-04,1.0e-14). Resetting to default value of 4.0e-04. ", MISC_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("gw_residual_fraction", &lc.gw_residual_fraction, 0.0, 1.0, 0.90,
            CHECK_AND_FIX, OPTIONAL,
            "The residual value specified by gw_residual_convergence_criterion is applied "
            "to this fraction of the total spectrum. ",
            "gw_residual_fraction must lie in the range (0.0,1.0). Resetting to default value of 0.90. ", MISC_OPTIONS|EXPERIMENTAL_OPTION);

    If.RegisterInputKey("hartree_rms_ratio", &lc.hartree_rms_ratio, 1000.0, DBL_MAX, 100000.0,
            CHECK_AND_FIX, OPTIONAL,
            "Ratio between target RMS for get_vh and RMS total potential. ",
            "hartree_rms_ratio must be in the range (1000.0, 1000000.0). Resetting to default value of 100000.0. ", POISSON_OPTIONS);

    If.RegisterInputKey("electric_field_magnitude", &lc.e_field, 0.0, DBL_MAX, 0.0,
            CHECK_AND_TERMINATE, OPTIONAL,
            "Magnitude of external electric field. ",
            "electric_field_magnitude must be a postive value. ");

    Ri::ReadVector<double> def_tddft_qpos({{0.0,0.0,0.0}});
    Ri::ReadVector<double> tddft_qpos;
    If.RegisterInputKey("tddft_qpos", &tddft_qpos, &def_tddft_qpos, 3, OPTIONAL,
            "cartesian coordinate of the point charge for tddft",
            "You must specify a triplet of (X,Y,Z) dimensions. ", TDDFT_OPTIONS);

    If.RegisterInputKey("tddft_qgau", &lc.tddft_qgau, 0.0, DBL_MAX, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "Gaussian parameter for point charge to Gaussian charge",
            "", TDDFT_OPTIONS);

    Ri::ReadVector<double> def_electric_field({{0.0,0.0,1.0}});
    Ri::ReadVector<double> electric_field;
    If.RegisterInputKey("electric_field_vector", &electric_field, &def_electric_field, 3, OPTIONAL,
            "Components of the electric field. ",
            "You must specify a triplet of (X,Y,Z) dimensions for the electric field vector. ");

    If.RegisterInputKey("Emin", &lc.Emin, -100.0, 100.0, -6.0,
            CHECK_AND_TERMINATE, OPTIONAL,
            "",
            "");

    If.RegisterInputKey("Emax", &lc.Emax, -100.0, 100.0, 0.0,
            CHECK_AND_TERMINATE, OPTIONAL,
            "",
            "");

    If.RegisterInputKey("energy_cutoff_parameter", &lc.cparm, 0.6, 1.0, 0.8,
            CHECK_AND_FIX, OPTIONAL,
            "",
            "energy_cutoff_parameter must be in the range (0.6,1.0). Resetting to default value of 0.8. ", PSEUDO_OPTIONS|EXPERT_OPTION);

    std::string Occup, Occdown;
    std::string Occ;

    If.RegisterInputKey("states_count_and_occupation_spin_up", &Occup, "",
            CHECK_AND_FIX, OPTIONAL,
            "Occupation string for spin up states. Format is the same as for states_count_and_occupation. Total number of states must match spin down occupation string.",
            "", OCCUPATION_OPTIONS);

    If.RegisterInputKey("states_count_and_occupation_spin_down", &Occdown, "",
            CHECK_AND_FIX, OPTIONAL,
            "Occupation string for spin down states. Format is the same as for states_count_and_occupation. Total number of states must match spin up occupation string.",
            "", OCCUPATION_OPTIONS);


    If.RegisterInputKey("states_count_and_occupation", &Occ, "",
            CHECK_AND_FIX, OPTIONAL,
            "Occupation string for states. Format for a system with 240 electrons and 20 unoccupied states would be. \"120 2.0 20 0.0\" ",
            "", OCCUPATION_OPTIONS);

    If.RegisterInputKey("kpoint_distribution", &pelc.pe_kpoint, -INT_MAX, INT_MAX, -1,
            CHECK_AND_FIX, OPTIONAL,
"This option affects kpoint parallelization. If there are M MPI procs then N = M/kpoint_distribution procs "
" are assigned to each kpoint. M must be evenly divisible by kpoint_distribution.", 
            "", CELL_OPTIONS);

    If.RegisterInputKey("time_reversal", &lc.time_reversal, true,
            "if false, no k -> -k symmetry", CONTROL_OPTIONS);
    If.RegisterInputKey("wannier90", &lc.wannier90, false,
            "set up informations for wannier90 interface", CONTROL_OPTIONS);
    If.RegisterInputKey("num_wanniers", &lc.num_wanniers, 0, INT_MAX, 0, 
            CHECK_AND_FIX, OPTIONAL, 
            "number of wannier functions to be used in wannier90 ", 
            "wannier90 is an open source package ");

    If.RegisterInputKey("wannier90_scdm_mu", &lc.wannier90_scdm_mu, -DBL_MAX, DBL_MAX, 0.0,
            CHECK_AND_FIX, OPTIONAL,
            "when wannier90 is used to build wannier functions, the energy window parameter",
            "it usually is around Fermi energy ", CONTROL_OPTIONS);
    If.RegisterInputKey("wannier90_scdm_sigma", &lc.wannier90_scdm_sigma, 0.0, DBL_MAX, 1.0,
            CHECK_AND_FIX, OPTIONAL,
            "when wannier90 is used to build wannier functions, the energy window parameter",
            "it is used for Gaussian or erfc filtering the band ");
    If.RegisterInputKey("wannier90_scdm", &lc.wannier90_scdm, -INT_MAX, 2, 0,
            CHECK_AND_FIX, OPTIONAL,
            "use scdm method to set the trial wannier functions ",
            "0: isolated, 1: gaussian, 2: erfc ");

    // Command line help request?
    bool cmdline = (std::find(ct.argv.begin(), ct.argv.end(), std::string("--help")) != ct.argv.end());
    bool markdown = (std::find(ct.argv.begin(), ct.argv.end(), std::string("--markdown")) != ct.argv.end());
    if(cmdline || markdown)
    {

        // These keys are parsed in another file but we want to register them here just so they show up in the
        // documentation. If you change the documentation in the other location make sure to change it here as well.
        std::string KpointArray;
        If.RegisterInputKey("kpoints", &KpointArray, "",
                         CHECK_AND_FIX, REQUIRED,
                         "Normally kpoints are specified using the kpoint_mesh and kpoint_is_shift options but one can also enter a list of kpoints and their weights with this option. If kpoint_mesh is not specified or this is a bandstructure calculation this is required otherwise it is optional. \n", "");

        If.RegisterInputKey("kpoints_bandstructure", &KpointArray, "",
                         CHECK_AND_FIX, OPTIONAL,
                         "List of kpoints to use in a bandstructure calculation. For more detailed information look at the github wiki page on kpoint calculations.\n", "");

        if(pct.imgpe == 0 && cmdline) WriteInputOptions(InputMap, std::string("cmdline"));
        if(pct.imgpe == 0 && markdown) WriteInputOptions(InputMap, std::string("markdown"));
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    If.LoadInputKeys();

    // Check items that require custom handling
    // Some hacks here to deal with code branches that are still in C
    if(!Description.length()) Description = "RMG electronic structure calculation.";
    lc.description = Description;

    if(!Infile.length()) Infile = "Waves/wave.out";
    std::strncpy(lc.infile, Infile.c_str(), sizeof(lc.infile)-1);
    MakeFullPath(lc.infile, pelc);

    if(!Outfile.length()) Outfile = "Waves/wave.out";
    std::strncpy(lc.outfile, Outfile.c_str(), sizeof(lc.outfile)-1);
    MakeFullPath(lc.outfile, pelc);

    lc.exx_int_file = ExxIntfile;
    MakeFullPath(lc.exx_int_file, pelc);

    lc.vdW_kernel_file = VdwKernelfile;
    MakeFullPath(lc.vdW_kernel_file, pelc);

    lc.nvme_weights_path = Weightsfile;
    if(lc.nvme_weights_path.length()) lc.nvme_weights_path.append("/");
    MakeFullPath(lc.nvme_weights_path, pelc);

    lc.qfunction_path = Qfunctionfile;
    if(lc.qfunction_path.length()) lc.qfunction_path.append("/");
    MakeFullPath(lc.qfunction_path, pelc);

    lc.nvme_work_path = Workfile;
    if(lc.nvme_work_path.length()) lc.nvme_work_path.append("/");
    MakeFullPath(lc.nvme_work_path, pelc);

    lc.nvme_orbitals_path = Orbitalfile;
    if(lc.nvme_orbitals_path.length()) lc.nvme_orbitals_path.append("/");
    MakeFullPath(lc.nvme_orbitals_path, pelc);

    if(!Infile_tddft.length()) Infile = "Waves/wave_tddft.out";
    std::strncpy(lc.infile_tddft, Infile_tddft.c_str(), sizeof(lc.infile_tddft)-1);
    MakeFullPath(lc.infile_tddft, pelc);

    if(!Outfile_tddft.length()) Outfile = "Waves/wave_tddft.out";
    std::strncpy(lc.outfile_tddft, Outfile_tddft.c_str(), sizeof(lc.outfile_tddft)-1);
    MakeFullPath(lc.outfile_tddft, pelc);


    lc.pseudo_dir = PseudoPath;
    if(lc.pseudo_dir.length()) lc.pseudo_dir.append("/");

    if((Occup.length() != 0) && (Occdown.length() != 0) && (Occ.length() == 0)) lc.spin_flag = 1;

    if((Occup.length() != 0) && (Occdown.length() != 0) && (Occ.length() != 0))
    {
        rmg_error_handler (__FILE__, __LINE__, "You have specified occupations for spin-up spin-down and non-spin cases which is ambiguous. Terminating.");
    }

    if(lc.spinorbit) ct.noncoll = true;
    lc.nspin = 1;
    lc.noncoll_factor = 1;
    if(lc.noncoll) 
    {
        lc.nspin = 4;
        lc.noncoll_factor = 2;
    }
    if(lc.spin_flag) {

        lc.nspin = 2;
        std::strncpy(lc.occupation_str_spin_up, Occup.c_str(), sizeof(lc.occupation_str_spin_up)-1);
        std::strncpy(lc.occupation_str_spin_down, Occdown.c_str(), sizeof(lc.occupation_str_spin_down)-1);

    }
    else {

        std::strncpy(lc.occupation_str, Occ.c_str(), sizeof(lc.occupation_str)-1);

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

    if (Verify ("lattice_units", "Angstrom", InputMap) && Verify ("crds_units", "Bohr", InputMap))
    {
        std::cout <<"WARNING:: lattice_unit Angstrom and  crds_unit Bohr"<<std::endl; 
    }
    if (Verify ("lattice_units", "Bohr", InputMap) && Verify ("crds_units", "Angstrom", InputMap))
    {
        std::cout <<"WARNING:: lattice_unit Bohr and  crds_unit Angstrom"<<std::endl; 
    }

    // Check if a lattice vector was specified and if not 
    if(lattice_vector.vals == def_lattice_vector.vals && ibrav == None)
        rmg_error_handler(__FILE__,__LINE__,"\nNeither a lattice_vector or a lattice type was specified. Terminating.\n");

    if(ibrav == None && lc.stress)
        rmg_error_handler(__FILE__,__LINE__, "\nStress not supported for arbitary lattice vectors. You need to specify bravais_lattice_type.: Terminating.\n");

    if(ibrav == None)
    {
        double Lunit = 1.0;
        if (Verify ("lattice_units", "Angstrom", InputMap))
            Lunit = A_a0;
        if (Verify ("lattice_units", "Alat", InputMap))
            Lunit = celldm[0];

        if(Lunit < 1.0e-10)    
            throw RmgFatalException() << "lattice units = Alat, but Alat(a_length) = 0.0" << "\n";

        for (int i = 0; i < 3; i++)
        {
            a0[i] = Lunit * lattice_vector.vals.at(i);
            a1[i] = Lunit * lattice_vector.vals.at(i+3);
            a2[i] = Lunit * lattice_vector.vals.at(i+6);
        }

        // Save copy of scaled input vectors
        Rmg_L.save_vectors(a0, a1, a2);

//        Rmg_L.lat2abc(a0, a1, a2);
        // Detects ibrav and generates a,b,c,cosab,cosac,cosbc
        int ibb = Rmg_L.lat2ibrav (a0, a1, a2);
        Rmg_L.rotate_vectors(a0, a1, a2);
        if(pct.imgpe==0) printf("Detected ibrav %d from lattice vectors.\n",ibb);

        // Sets up celldm
        Rmg_L.abc2celldm();

        for(int i=0;i < 6;i++) {celldm[i] = Rmg_L.get_celldm(i);}
        ibrav = Rmg_L.get_ibrav_type();
        Rmg_L.latgen(celldm, &omega, a0, a1, a2, true);

    }
    else
    {
        // Transform to atomic units, which are used internally if input is in angstrom 
        if (Verify ("lattice_units", "Angstrom", InputMap))
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

        if (Verify ("bravais_lattice_type", "Cubic Primitive", InputMap) )
        {
            if( std::abs(celldm[1] - 1.0 ) > 1.0e-5 || std::abs(celldm[2] - 1.0) > 1.0e-5)

                throw RmgFatalException() << "a, b, c is not consistent with Cubic Primitive\n";
        }
    
        // Lattice vectors are orthogonal except for Hex which is setup inside latgen
        celldm[3] = 0.0;
        celldm[4] = 0.0;
        celldm[5] = 0.0;

        // Set up the lattice vectors
        Rmg_L.set_ibrav_type(ibrav);
        Rmg_L.latgen(celldm, &omega, a0, a1, a2, false);
        Rmg_L.save_vectors(Rmg_L.a0, Rmg_L.a1, Rmg_L.a2);
    }

    int NX_GRID = WavefunctionGrid.vals.at(0);
    int NY_GRID = WavefunctionGrid.vals.at(1);
    int NZ_GRID = WavefunctionGrid.vals.at(2);

    lc.dipole_corr[0] = DipoleCorrection.vals.at(0);
    lc.dipole_corr[1] = DipoleCorrection.vals.at(1);
    lc.dipole_corr[2] = DipoleCorrection.vals.at(2);

    for(int i = 0; i < 9; i++) lc.cell_movable[i] = Cell_movable.vals.at(i);

    CheckAndTerminate(NX_GRID, 1, INT_MAX, "The value given for the global wavefunction grid X dimension is " + boost::lexical_cast<std::string>(NX_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(NY_GRID, 1, INT_MAX, "The value given for the global wavefunction grid Y dimension is " + boost::lexical_cast<std::string>(NY_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(NZ_GRID, 1, INT_MAX, "The value given for the global wavefunction grid Z dimension is " + boost::lexical_cast<std::string>(NZ_GRID) + " and only postive values are allowed.");

    pelc.pe_x = ProcessorGrid.vals.at(0);
    pelc.pe_y = ProcessorGrid.vals.at(1);
    pelc.pe_z = ProcessorGrid.vals.at(2);

    CheckAndTerminate(pelc.pe_x, 1, INT_MAX, "The value given for the global processor grid X dimension is " + boost::lexical_cast<std::string>(NX_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_y, 1, INT_MAX, "The value given for the global processor grid Y dimension is " + boost::lexical_cast<std::string>(NY_GRID) + " and only postive values are allowed.");
    CheckAndTerminate(pelc.pe_z, 1, INT_MAX, "The value given for the global processor grid Z dimension is " + boost::lexical_cast<std::string>(NZ_GRID) + " and only postive values are allowed.");

    ct.tddft_qpos[0] = tddft_qpos.vals.at(0);
    ct.tddft_qpos[1] = tddft_qpos.vals.at(1);
    ct.tddft_qpos[2] = tddft_qpos.vals.at(2);
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
        std::cout << "You have set freeze_occupied=true so potential acceleration is disabled." << std::endl;
    }

    // Debug code
    //for(auto it = InputMap.begin();it != InputMap.end(); ++it) {
    //    std::pair<std::string, InputKey*> Check = *it;
    //    InputKey *CheckKey = it->second;
    //std::cout << Check.first << " = " << CheckKey->Print() << std::endl;
    //}

    // Set up energy output units
    lc.energy_output_conversion[0] = 1.0;
    lc.energy_output_conversion[1] = 2.0;
    lc.energy_output_string[0] = "Ha";
    lc.energy_output_string[1] = "Ry";

    ct.use_vdwdf_finegrid = Verify ("vdwdf_grid_type", "Fine", InputMap);

    if(states_list.length() != 0 )
    {
        size_t pos = 0;
        while (states_list.size() > 0)
        {
            pos = states_list.find(",");
            pos = std::min(pos, states_list.size());
            std::string stcount = states_list.substr(0, pos);
            states_list.erase(0, pos + 1);

            pos = stcount.find("-");
            if(pos > stcount.length() )
            {
                
                ct.cube_states_list.push_back(std::stoi(stcount));
            }
            else
            {
                int st_start = std::stoi(stcount.substr(0,pos));
                int st_end = std::stoi(stcount.substr(pos+1));
                
                for(int st = st_start; st <= st_end; st++)
                    ct.cube_states_list.push_back(st);
            }
        }
    }

    {
        const char *buf = stm_bias_str.c_str();
        char *end;
        for (double bias = strtod(buf, &end); buf != end; bias = strtod(buf, &end) )
        {
            buf = end;
            ct.stm_bias_list.push_back(bias);
        }
    }
    {
        const char *buf = stm_height_str.c_str();
        char *end;
        for (double height = strtod(buf, &end); buf != end; height = strtod(buf, &end) )
        {
            buf = end;
            ct.stm_height_list.push_back(height);
        }
    }
    if(ct.verbose && pct.imgpe == 0)
    {
        for(size_t st = 0; st < ct.cube_states_list.size(); st++)
            printf("\n state list %d", ct.cube_states_list[st]);
    }

    if((ct.kohn_sham_solver == DAVIDSON_SOLVER) && Verify("charge_mixing_type","Linear", InputMap))
    {
        rmg_error_handler (__FILE__, __LINE__, "\nError. You have selected Linear Mixing with the Davidson kohn-sham solver\nwhich is not valid. Please change to Broyden or Pulay mixing. Terminating.\n\n");
    }

}
