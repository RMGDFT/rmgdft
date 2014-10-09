

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;
#include <iostream> 
#include <fstream>
#include <sstream>
#include <iterator>
#include <string> 
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

 


std::string preprocess_config_file(char *cfile)
{
    std::string config_file(cfile);
    std::string outbuf;
    std::string tbuf;

    std::ifstream ifs(cfile);

    // First pass to clean it up
    for(std::string line; std::getline(ifs, line); ) {

        // Strip leading and trailing whitespace
        boost::trim(line);

        // Filter lines and only include non-empty and non-comment lines
        std::size_t len = line.size();
        if(len > 0) {
            if( line.compare(0, 1, "#") ) {
                // Chop trailing comments. Might need to fix this for # characters embedded in strings
                std::size_t found;
                while(std::string::npos != (found = line.find_last_of("#"))) line.resize(found);
                boost::trim(line);
                tbuf = tbuf + line + "\n";
            }
        }

    }
    // Split into individual lines again
    std::string delims = "\n";
    std::vector<std::string> lines;
    std::vector<std::string> lines1;
    std::vector<std::string>::iterator it;
    boost::algorithm::split( lines, tbuf, boost::is_any_of(delims), boost::token_compress_on );
    tbuf.clear();
    int idx = -1;
    bool join = false;
    for (it = lines.begin(); it != lines.end(); ++it) {
        std::string line = *it;

        if(join) {
            lines1.at(idx) = lines1.at(idx) + line;
            join = false;
        }

        // If line contains "=" then it's a new key-value pair
        if(std::string::npos != line.find("=")) {
            idx++;
            // new key-value pair
            lines1.push_back(line);
            join = false;
        }

        // Forward line
        std::vector<std::string>::iterator it1;
        it1 = it + 1;

        // Any quotes in this line?
        std::size_t f1 = line.find("\"", 1);
        if(f1 != std::string::npos) {
            std::size_t f2 = line.find("\"", f1 + 1);
            // If two quotes then we are done, if one quote then join with next line unless it contains an =
            if(f2 != std::string::npos) {
                join = false;
            }
            else {
                join = true;
                if(it1 != lines.end()) {
                    std::string fline = *it1;
                    f1 = fline.find("=", 1);
                    if(f1 != std::string::npos) join = false;
                }
            }
        }
        else {
            // no quotes so join unless next line contains an equal sign
            join = true;
            if(it1 != lines.end()) {
                std::string fline = *it1;
                f1 = fline.find("=", 1);
                if(f1 != std::string::npos) join = false;
            }
        }
    }

    for (it = lines1.begin(); it != lines1.end(); ++it) {
        outbuf = outbuf + *it + "\n";
    }

    //std::cout << outbuf << std::endl;
    return outbuf;

}


// Custom validator for boost program options. 
namespace RmgInput { 

    template <typename VectorType>
    class ReadVector
    {
        public:
            std::vector<VectorType> vals;
    };

    void validate(boost::any& v, const std::vector<std::string>& values, ReadVector<int>*, int)
    {
        ReadVector<int> A;
        po::validators::check_first_occurrence(v);
        const std::string& s = po::validators::get_single_string(values);
        std::string t1 = s;
        boost::trim_if(t1, boost::algorithm::is_any_of("\" \t"));

        while(t1.size() > 0) {
            A.vals.push_back( std::atoi(t1.c_str()));
            size_t f1 = t1.find_first_of(" \t");
            if(f1 == std::string::npos) break;
            t1 = t1.substr(f1, std::string::npos);
            boost::trim(t1);
        }
        v = A;
    }

    void validate(boost::any& v, const std::vector<std::string>& values, ReadVector<double>*, double)
    {
        ReadVector<double> A;
        po::validators::check_first_occurrence(v);
        const std::string& s = po::validators::get_single_string(values);
        std::string t1 = s;
        boost::trim_if(t1, boost::algorithm::is_any_of("\" \t"));

        while(t1.size() > 0) {
            A.vals.push_back( std::atof(t1.c_str()));
            size_t f1 = t1.find_first_of(" \t");
            if(f1 == std::string::npos) break;
            t1 = t1.substr(f1, std::string::npos);
            boost::trim(t1);
        }
        v = A;
    }

}

namespace Ri = RmgInput;

void ReadCommon(int argc, char *argv[], char *cfile)
{

    std::string sfile = preprocess_config_file(cfile);
    std::string pa_constant_step, pa_poisson_step;
   
    Ri::ReadVector<int> ProcessorGrid;
    Ri::ReadVector<int> CoarseGrid;
    int FG_RATIO;
 
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "print version string")
        ("help", "produce help message")
    ;

    po::options_description control("RMG control options");
    control.add_options()
       //("description", po::value<std::string>(&ii)->default_value("RMG RUN"), "description of run")

        ("processor_grid", po::value<Ri::ReadVector<int> >(&ProcessorGrid), "Three-D (x,y,z) layout of the MPI processes.")
        ("coarse_grid", po::value<Ri::ReadVector<int> >(&CoarseGrid), "Three-D (x,y,z) dimensions of the global wavefunction grid.")
        ("potential_grid_refinement=", po::value<int>(&FG_RATIO)->default_value(2), "Ratio of the fine grid to the coarse grid.")
        ("potential_acceleration_constant_step", po::value<double>(&ct.potential_acceleration_constant_step)->default_value(0.0), "Time step used for constant potential acceleration.")
        ("potential_acceleration_poisson_step", po::value<double>(&ct.potential_acceleration_poisson_step)->default_value(0.0), "Time step used for poisson potential acceleration.")
        ("ionic_time_step", po::value<double>(&ct.iondt)->default_value(50.0), "Ionic time step for use in molecular dynamics and structure optimizations.")
        ("unoccupied_states_per_kpoint", po::value<int>(&ct.num_unocc_states)->default_value(10), "The number of unoccupied orbitals.")
        ("folded_spectrum", po::value<bool>(&ct.use_folded_spectrum)->default_value(false), "Use folded spectrum.")
        ("sort_wavefunctions", po::value<bool>(&ct.sortflag)->default_value(false), "Sort wavefunctions by eigenvalue. Not needed if using subspace diagonalization.")
        ("period_of_diagonalization", po::value<int>(&ct.diag)->default_value(1), "Diagonalization period (per scf step).")
        ("end_diagonalization_step", po::value<int>(&ct.end_diag)->default_value(9999999), "Stop diagonalizing after end_diag steps.")
        ("initial_diagonalization", po::value<bool>(&ct.initdiag)->default_value(false), "Perform initial subspace diagonalization.")
        ("charge_pulay_special_metrics", po::value<bool>(&ct.charge_pulay_special_metrics)->default_value(false), "Flag to test whether or not the modified metrics should be used in Pulay mixing.")
        ("write_pdos", po::value<bool>(&ct.pdos_flag)->default_value(false), "Flag to write partial density of states.")
        ("equal_initial_density", po::value<bool>(&ct.init_equal_density_flag)->default_value(false), "Specifies whether to set initial up and down density to be equal.")
        ("threads_per_node", po::value<int>(&ct.THREADS_PER_NODE)->default_value(-1),"")
        ("charge_pulay_order", po::value<int>(&ct.charge_pulay_order)->default_value(5),"")
        ("charge_pulay_refresh", po::value<int>(&ct.charge_pulay_refresh)->default_value(0),"")
        ("max_scf_steps", po::value<int>(&ct.max_scf_steps)->default_value(500),"")
        ("write_data_period", po::value<int>(&ct.checkpoint)->default_value(5),"")
        ("write_eigvals_period", po::value<int>(&ct.write_eigvals_period)->default_value(5),"")
        ("max_md_steps", po::value<int>(&ct.max_md_steps)->default_value(100),"")
        ("hartree_max_sweeps", po::value<int>(&ct.hartree_max_sweeps)->default_value(100),"")
        ("hartree_min_sweeps", po::value<int>(&ct.hartree_min_sweeps)->default_value(5),"")
        ("kohn_sham_pre_smoothing", po::value<int>(&ct.eig_parm.gl_pre)->default_value(2),"")
        ("kohn_sham_post_smoothing", po::value<int>(&ct.eig_parm.gl_pst)->default_value(1),"")
        ("kohn_sham_mucycles", po::value<int>(&ct.eig_parm.mucycles)->default_value(1),"")
        ("kohn_sham_fd_order", po::value<int>(&ct.kohn_sham_fd_order)->default_value(6),"")
        ("poisson_pre_smoothing", po::value<int>(&ct.poi_parm.gl_pre)->default_value(3),"")
        ("poisson_post_smoothing", po::value<int>(&ct.poi_parm.gl_pst)->default_value(3),"")
        ("poisson_mucycles", po::value<int>(&ct.poi_parm.mucycles)->default_value(1),"")
        ("poisson_coarsest_steps", po::value<int>(&ct.poi_parm.coarsest_steps)->default_value(80),"")
        ("kohn_sham_mg_levels", po::value<int>(&ct.eig_parm.levels)->default_value(1),"")
        ("poisson_mg_levels", po::value<int>(&ct.poi_parm.levels)->default_value(-1),"")
        ("fine_grid_non_local_pp", po::value<int>(&ct.nxfgrid)->default_value(4),"")
        ("scalapack_block_factor", po::value<int>(&ct.scalapack_block_factor)->default_value(32),"")
        ("E_POINTS", po::value<int>(&ct.E_POINTS)->default_value(201),"")
        ("max_rmg_steps", po::value<int>(&ct.max_rmg_steps)->default_value(1),"")
        ("md_number_of_nose_thermostats", po::value<int>(&ct.nose.m)->default_value(5),"")
        ("dynamic_time_delay", po::value<int>(&ct.relax_steps_delay)->default_value(5),"")
        ("dynamic_time_counter", po::value<int>(&ct.relax_steps_counter)->default_value(0),"")
        ("scf_steps_offset", po::value<int>(&ct.scf_steps)->default_value(0),"")
        ("total_scf_steps_offset", po::value<int>(&ct.total_scf_steps)->default_value(0),"")
        ("md_steps_offset", po::value<int>(&ct.md_steps)->default_value(0),"")
        ("b_spline_order", po::value<int>(&ct.interp_order)->default_value(5),"")
        ("b_spline_trade_order", po::value<int>(&ct.interp_trade)->default_value(3),"")

    ;

    po::options_description config_file_options;
    config_file_options.add(control);

    po::variables_map vm;
    std::stringstream ss;
    ss << sfile;

    auto parsedOptions = po::parse_config_file(ss, config_file_options, true);
    po::store(parsedOptions, vm);
    po::notify(vm);
    auto unregistered = po::collect_unrecognized(parsedOptions.options, po::include_positional);

    try {
        pct.pe_x = ProcessorGrid.vals.at(0);
        pct.pe_y = ProcessorGrid.vals.at(1);
        pct.pe_z = ProcessorGrid.vals.at(2);
    }
    catch (const std::out_of_range& oor) {
        std::cerr << "You must specify a triplet of (X,Y,Z) dimensions for the processor grid." << std::endl;
        rmg_error_handler(__FILE__, __LINE__, "Terminating.\n");
    }

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

    // Set grid info up
    Rmg_G = new BaseGrid(NX_GRID, NY_GRID, NZ_GRID, pct.pe_x, pct.pe_y, pct.pe_z, 0, FG_RATIO);



    if (vm.count("ionic_time_step")) {
        std::cout << "AAionic_time_step=" << ct.iondt << std::endl;
    }
    if (vm.count("potential_acceleration_constant_step")) {
        std::cout << "AApotential_acceleration_constant_step=" << ct.potential_acceleration_constant_step << std::endl;
    }
    if (vm.count("unoccupied_states_per_kpoint")) {
        std::cout << "AAunoccupied_states_per_kpoint=" << ct.num_unocc_states << std::endl;
    }
    if (vm.count("description")) {
        //std::cout << "AAdescription=" << ii << std::endl;
    }



}
