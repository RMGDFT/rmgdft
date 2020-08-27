
#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>



#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "MapElements.h"
#include "RmgException.h"
#include "transition.h"
#include "InternalPseudo.h"

static std::unordered_map<std::string, int> atomic_map = {
        {"s", 0},
        {"p", 1},
        {"d", 2},
        {"f", 3},
        {"g", 4}};

using boost::property_tree::ptree;

// Reads a pseudopotential stored in XML format into our internal data structures.
void LoadXmlPseudo(SPECIES *sp)
{

    ptree xml_tree;
    char *pp_buffer = NULL;
    int pp_buffer_len;
    int max_nlprojectors = 0;
    std::stringstream ss; 
    double_2d_array ddd0;  // Used to read in the PP_DIJ
    double_2d_array qqq;   // Used to read in the norms of the augmentation functions (PP_Q)
    std::string Msg;
    sp->max_l = 0;

    // PP format places occupations in the xml attributes of the potentials but we need
    // them in the wavefunctions so we use this map to associate them.
    std::unordered_map<int, double> occupation_map;

    // Open on one pe and read entire file into a character buffer
    if(pct.imgpe == 0) {

        // Check for file existence
        boost::filesystem::path pp_filepath(sp->pseudo_filename);
        if( !boost::filesystem::exists(pp_filepath) ) {

            Msg = "Pseudopotential file " + boost::lexical_cast<std::string>(sp->pseudo_filename) + " does not exist.\n";

        }
        else {

            try {
                std::ifstream pp_file;
                pp_file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
                pp_file.open(sp->pseudo_filename);
                pp_file.seekg(0, std::ios::end);
                pp_buffer_len = pp_file.tellg();
                pp_file.seekg(0, std::ios::beg);
                pp_buffer = new char[pp_buffer_len + 1]();
                pp_file.read(pp_buffer, pp_buffer_len);       // read the whole file into the buffer
                pp_file.close();
            }
            // Catch any file io errors and rethrow later with our own error message
            catch (std::exception &e) {
                Msg = "Unable to read from pseudopotential file " + boost::lexical_cast<std::string>(sp->pseudo_filename) + "\n";
            }
        }

    }

    int openfail = Msg.length();
    MPI_Bcast(&openfail, 1, MPI_INT, 0, pct.img_comm);
    if(openfail) 
        throw RmgFatalException() << Msg << " in " << __FILE__ << " at line " << __LINE__ << " file name: " + boost::lexical_cast<std::string>(sp->pseudo_filename) + "\n";


    // Send it to everyone else
    MPI_Bcast (&pp_buffer_len, 1, MPI_INT, 0, pct.img_comm);
    if(pct.imgpe != 0) {
        pp_buffer = new char[pp_buffer_len + 1]();
    }
    MPI_Bcast (pp_buffer, pp_buffer_len, MPI_CHAR, 0, pct.img_comm);
    std::string pp_string(pp_buffer);
    ss << pp_string;


    // Get the compulsory stuff first
    read_xml(ss, xml_tree);
    sp->INFO = xml_tree.get<std::string>("pseudo.header"); 

    // Check for relativistic and thrown an error if found
    std::string s_is_relativistic = xml_tree.get<std::string>("pseudo.header.<xmlattr>.relativistic");
    boost::to_upper(s_is_relativistic);
    if(s_is_relativistic == "YES")
        throw RmgFatalException() << "RMG does not support fully relativistic pseudopotentials. Terminating.\n";

    // Atomic symbol, mass, number and zvalence and mesh size
    std::string atomic_symbol = xml_tree.get<std::string>("pseudo.header.<xmlattr>.symbol");
    boost::trim(atomic_symbol);

    // Maybe check symbols here
    sp->atomic_mass = GetAtomicMass(atomic_symbol);
    sp->atomic_number = GetAtomicNumber(atomic_symbol);
    sp->zvalence = xml_tree.get<double>("pseudo.header.<xmlattr>.zval");
    if(sp->zvalence > ct.max_zvalence) ct.max_zvalence = sp->zvalence;

    // Store functional information for later processing
    sp->functional = xml_tree.get<std::string>("pseudo.header.<xmlattr>.xc-functional-type");
    std::string xctype = xml_tree.get<std::string>("pseudo.header.<xmlattr>.xc-functional-parametrization");
    if(xctype == "Perdew-Burke-Ernzerhof")
    {
       sp->functional = "PBE";
    }

    // Setup the radial mesh
    std::string grid_type = xml_tree.get<std::string>("pseudo.grid.<xmlattr>.type");
    int r_total = xml_tree.get<int>("pseudo.grid.<xmlattr>.npts");
    double r_i = xml_tree.get<double>("pseudo.grid.<xmlattr>.ri");
    double r_f = xml_tree.get<double>("pseudo.grid.<xmlattr>.rf");

    if(grid_type == std::string("linear"))
    {
        sp->rg_points = r_total - 1;
        sp->gtype = LINEAR_GRID;
        double delta_r = (r_f - r_i) / (double)sp->rg_points;
        sp->r = new double[sp->rg_points];
        sp->rab = new double[sp->rg_points];
        for(int i = 0;i < sp->rg_points;i++) {
            sp->r[i] = delta_r * (double)(i + 1);
            sp->rab[i] = delta_r;
        }
    }

    // Determine log mesh parameters directly from the mesh
    sp->aa = (sp->r[0] * sp->r[0]) / (sp->r[1] - 2 * sp->r[0]);
    sp->bb = log (sp->r[1] / sp->r[0] - 1);

    // Only norm-conserving supported
    sp->is_norm_conserving = true;

    // Core correction flag
    std::string s_core_correction = xml_tree.get<std::string>("pseudo.header.<xmlattr>.core-corrections");
    boost::to_upper(s_core_correction);
    if(!s_core_correction.compare(0,2,"NO")) sp->nlccflag = false;
    if(!s_core_correction.compare(0,3,"YES")) sp->nlccflag = true;

    // First pass for semi-local find the local potential

    BOOST_FOREACH( ptree::value_type const& s, xml_tree.get_child("pseudo.semilocal") ) 
    {
        int lval = 0;
        sp->local = xml_tree.get<int>("pseudo.semilocal.<xmlattr>.l-local", -3);
        sp->is_semi_local = true;
        if( s.first == "vps" ) 
        {
            std::string l = s.second.get("<xmlattr>.l", "s");
            lval = atomic_map[l];

            // Get occupation here needed for wavefunctions below
            occupation_map[lval] = s.second.get<double>("<xmlattr>.occupation", 2.0);
            if(lval == sp->local)
            {
                std::string vdata = s.second.get<std::string>("radfunc.data");
                sp->vloc0 = UPF_str_to_double_array(vdata, r_total, 1);
                for(int idx = 0;idx < sp->rg_points;idx++) sp->vloc0[idx] /= sp->r[idx];
            } 
        }
    }

    // Second pass for semi-local to take care of the rest
    if(sp->is_semi_local)
    {
        sp->kkbeta = sp->rg_points;
        sp->nbeta = 0;
        sp->is_ddd_diagonal = true;
        BOOST_FOREACH( ptree::value_type const& s, xml_tree.get_child("pseudo.semilocal") )
        {
            int lval = 0;
            if( s.first == "vps" )
            {
                std::string l = s.second.get("<xmlattr>.l", "s");
                lval = atomic_map[l];
                if(lval != sp->local)
                {

                    // Generate the difference potentials.
                    std::string vdata = s.second.get<std::string>("radfunc.data");
                    sp->dVl[sp->nbeta] = UPF_str_to_double_array(vdata, r_total, 1);
                    sp->dVl_l[sp->nbeta] = lval;
                    for(int idx=0;idx < sp->rg_points;idx++)
                    {
                        sp->dVl[sp->nbeta][idx] = sp->dVl[sp->nbeta][idx] / sp->r[idx] - sp->vloc0[idx];
                    }
                    sp->nbeta++;
                }
            }
        }
    }

    int iwf = 0;
    sp->num_atomic_waves_m = 0;
    sp->atomic_wave = new double *[MAX_INITWF];
    sp->aradius = new double [MAX_INITWF];
    sp->awave_lig = new double *[MAX_INITWF];
    sp->atomic_wave_l = new int [MAX_INITWF];
    sp->atomic_wave_j = new double [MAX_INITWF]();
    sp->atomic_wave_oc.resize(MAX_INITWF);
    sp->atomic_wave_oc.assign(MAX_INITWF, 0.0);
    sp->atomic_wave_energy.resize(MAX_INITWF);
    sp->atomic_wave_energy.assign(MAX_INITWF, 0.0);
    sp->atomic_rho = new double[sp->rg_points]();

    BOOST_FOREACH( ptree::value_type const& s, xml_tree.get_child("pseudo.pseudowave-functions") )
    {
        int lval=0;
        if( s.first == "pswf" )
        {
            std::string l = s.second.get<std::string>("<xmlattr>.l", "s");
            lval = atomic_map[l];
            sp->atomic_wave_l[iwf] = lval;
            sp->num_atomic_waves_m += 2*sp->atomic_wave_l[iwf] + 1;
            std::string vdata = s.second.get<std::string>("radfunc.data");
            sp->atomic_wave[iwf] = UPF_str_to_double_array(vdata, r_total, 1);
            sp->awave_lig[iwf] = new double[MAX_LOCAL_LIG]();

            // Remove extra factors of r
            for(int ix = 0;ix < sp->rg_points;ix++) sp->atomic_wave[iwf][ix] /= sp->r[ix];

            sp->atomic_wave_oc[iwf] = occupation_map[lval];
            sp->aradius[iwf] = 12.0;

            // Accumulate charge for atomic rho
            for(int idx=0;idx < sp->rg_points;idx++)
            {
                sp->atomic_rho[idx] += sp->atomic_wave[iwf][idx]*sp->atomic_wave[iwf][idx] * sp->atomic_wave_oc[iwf]/4.0/PI;
            }
            iwf++;
        }
    }
    sp->num_atomic_waves = iwf;

    // Next generate the Kleinman-Bylander projectors.
    // Diagonals of ddd0 array are the KB normalization coefficients.
    ddd0.resize(boost::extents[sp->nbeta][sp->nbeta]);
    for (int j = 0; j < sp->nbeta; j++)
    {
        for (int k = 0; k < sp->nbeta; k++) ddd0[j][k] = 0.0;
    }

    int nb = 0;
    for(int ip=0;ip < sp->num_atomic_waves;ip++)
    {
        if(sp->atomic_wave_l[ip] > ct.max_l) ct.max_l = sp->atomic_wave_l[ip];
        if(sp->atomic_wave_l[ip] == sp->local) continue; 

        int vl = 0;
        for(vl=0;vl <= MAX_L;vl++) if(sp->atomic_wave_l[ip] == sp->dVl_l[vl]) break;
        if(vl > MAX_L)
            throw RmgFatalException() << "Problem detected with pseudopotential file. MAX_L too large. Terminating.\n";

        sp->beta[ip] = new double[sp->rg_points];

        for(int idx=0;idx < sp->rg_points;idx++) sp->beta[ip][idx] = sp->atomic_wave[ip][idx] * sp->dVl[vl][idx];
        sp->llbeta[ip] = sp->atomic_wave_l[ip];
        if(sp->llbeta[ip] > sp->max_l) sp->max_l = sp->llbeta[ip];

        // Evaluate the normalization constant
        double *work = new double[sp->rg_points]();
        for(int idx=0;idx<sp->rg_points;idx++) 
            work[idx] = sp->beta[ip][idx]*sp->atomic_wave[ip][idx];

        double sum = FOUR * PI * radint1 (work, sp->r, sp->rab, sp->rg_points);
        ddd0[nb][nb] = FOUR * PI / sum;
        delete [] work;
        nb++;
    }

#if 0
   
    if(sp->nlccflag) {
        std::string PP_NLCC = xml_tree.get<std::string>("UPF.PP_NLCC");
        sp->rspsco = UPF_read_mesh_array(PP_NLCC, r_total, ibegin);
        for(int ix = 0;ix < sp->rg_points;ix++) sp->rspsco[ix] = sp->rspsco[ix] * 4.0 * PI;
    }
#endif


    sp->nqf=0;
    sp->nlc=0;

    int ihmax = 0;
    for (int j = 0; j < sp->nbeta; j++)
    {
        int l = sp->llbeta[j];
        for (int k = 0; k < 2 * l + 1; k++)
        {
            ++ihmax;
        }
    }

    for (int j = 0; j < ihmax; j++)
    {
        sp->nhtol.push_back(0);
        sp->nhtom.push_back(0);
        sp->indv.push_back(0);
        sp->nh_l2m.push_back(0);
    }

    int ih = 0;
    for (int j = 0; j < sp->nbeta; j++)
    {
        int l = sp->llbeta[j];
        for (int k = 0; k < 2 * l + 1; k++)
        {
            sp->nhtol[ih] = l;
            sp->nhtom[ih] = k;
            sp->indv[ih] = j;
            sp->nh_l2m[ih] = l*l + k;
            ++ih;
        }
    }
    sp->nh = ih;
    if (ih > max_nlprojectors)
        max_nlprojectors = ih;

    sp->ddd0.resize(boost::extents[ih][ih]);
    sp->qqq.resize(boost::extents[ih][ih]);
    for (int j = 0; j < ih; j++)
    {
        for (int k = 0; k < ih; k++)
        {
            if ((sp->nhtol[j] == sp->nhtol[k]) && (sp->nhtom[j] == sp->nhtom[k]))
            {
                sp->ddd0[j][k] = ddd0[sp->indv[j]][sp->indv[k]];
                sp->qqq[j][k] = qqq[sp->indv[j]][sp->indv[k]];
                /*                          sp->ddd[j][k]=ddd[indv[j]][indv[k]]; */
            }
            else
            {
                sp->ddd0[j][k] = 0.0;
                sp->qqq[j][k] = 0.0;
                /*                          sp->ddd[j][k]=0.0; */
            }
        }
    }


    // Set the maximum number of non-local projecters needed
    if(max_nlprojectors > ct.max_nl) 
        ct.max_nl = max_nlprojectors;

    // Optional stuff next

    sp->description = xml_tree.get<std::string>("pseudo.header.<xmlattr>.comment", "Pseudopotential");

    // Stuff not present in the UPF format that RMG requires. 
    // We need to find a consistent way of automatically setting these.
    sp->rc = fabs(2.0 * sp->zvalence / sqrt(PI) / sp->vloc0[0]);
    sp->rc = std::max(sp->rc, 0.5);
    sp->rc = std::min(sp->rc, 1.5);
    //sp->rc = 1.0;
    sp->lradius = 8.5;
    sp->gwidth = 8.0;
    sp->rwidth = 15.0; 
    sp->agwidth = 10.0;
    sp->arwidth = 25.0;

    // Leftover initializations
    sp->mill_radius = 9.0;

    if(pp_buffer) delete [] pp_buffer;
}


