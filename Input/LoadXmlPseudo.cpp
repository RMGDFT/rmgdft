
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


#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "MapElements.h"
#include "RmgException.h"
#include "transition.h"
#include "InternalPseudo.h"


using boost::property_tree::ptree;

// Reads a pseudopotential stored in XML format into our internal data structures.
void LoadXmlPseudo(SPECIES *sp)
{

    int ibegin=0, iend=0;
    ptree xml_tree;
    char *pp_buffer = NULL;
    int pp_buffer_len;
    int max_nlprojectors = 0;
    int l_max;
    std::stringstream ss; 
    double  ddd0[MAX_NL][MAX_NL];  // Used to read in the PP_DIJ
    double qqq[MAX_NL][MAX_NL];    // Used to read in the norms of the augmentation functions (PP_Q)

    std::string Msg;


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
    std::string PP_INFO = xml_tree.get<std::string>("pseudo.header"); 
    sp->INFO = new char[PP_INFO.size() + 1]();
    std::strncpy(sp->INFO, PP_INFO.c_str(), PP_INFO.size());

    // Atomic symbol, mass, number and zvalence and mesh size
    std::string atomic_symbol = xml_tree.get<std::string>("pseudo.header.<xmlattr>.symbol");
    boost::trim(atomic_symbol);

    // Maybe check symbols here
    sp->atomic_mass = GetAtomicMass(atomic_symbol);
    sp->atomic_number = GetAtomicNumber(atomic_symbol);
    sp->zvalence = xml_tree.get<double>("pseudo.header.<xmlattr>.zval");
    if(sp->zvalence > ct.max_zvalence) ct.max_zvalence = sp->zvalence;

    // Store functional information for later processing
    std::string PP_FUNCTYPE = xml_tree.get<std::string>("pseudo.header.<xmlattr>.xc-functional-type");
    std::strncpy(sp->functional_type, PP_FUNCTYPE.c_str(), sizeof(sp->functional_type)-1);
    std::string PP_FUNC = xml_tree.get<std::string>("pseudo.header.<xmlattr>.xc-functional-parametrization");
    std::strncpy(sp->functional, PP_FUNC.c_str(), sizeof(sp->functional)-1);


    // Read in the radial mesh and keep values between 1.0e-05 < r < 50.0
    // adjusting mesh size accordingly
    std::string grid_type = xml_tree.get<std::string>("pseudo.grid.<xmlattr>.type");
    int r_total = xml_tree.get<int>("pseudo.grid.<xmlattr>.npts");
    double r_i = xml_tree.get<double>("pseudo.grid.<xmlattr>.ri");
    double r_f = xml_tree.get<double>("pseudo.grid.<xmlattr>.rf");

    if(grid_type == std::string("linear"))
    {
        sp->rg_points = r_total;
        double delta_r = (r_f - r_i) / (double)(r_total - 1);
        sp->r = new double[sp->rg_points];
        for(int i = 0;i < r_total;i++) {
            sp->r[i] = delta_r * (double)i;
        }
    }

    // Only norm-conserving supported
    sp->is_norm_conserving = true;

    // Core correction flag
    std::string s_core_correction = xml_tree.get<std::string>("pseudo.header.<xmlattr>.core-corrections");
    boost::to_upper(s_core_correction);
    if(!s_core_correction.compare(0,2,"NO")) sp->nlccflag = false;
    if(!s_core_correction.compare(0,3,"YES")) sp->nlccflag = true;

#if 0
    l_max = xml_tree.get<int>("pseudo.header.<xmlattr>.l_max");
    sp->kkbeta = sp->rg_points;
   
    // Check for full relativistic and thrown an error if found
    std::string s_is_relativistic = xml_tree.get<std::string>("pseudo.header.<xmlattr>.relativistic");
    boost::to_upper(s_is_relativistic);
    if(!s_is_relativistic.compare(0,4,"FULL")) {
        throw RmgFatalException() << "RMG does not support fully relativistic pseudopotentials. Terminating.\n";
    }


    // Determine log mesh parameters directly from the mesh
    sp->aa = (sp->r[0] * sp->r[0]) / (sp->r[1] - 2 * sp->r[0]);
    sp->bb = log (sp->r[1] / sp->r[0] - 1);

    // Read in rab and convert it into a C style array
    std::string PP_RAB = xml_tree.get<std::string>("UPF.PP_MESH.PP_RAB");
    sp->rab = UPF_read_mesh_array(PP_RAB, r_total, ibegin);

    // Local potential
    std::string PP_LOCAL = xml_tree.get<std::string>("UPF.PP_LOCAL");
    sp->vloc0 = UPF_read_mesh_array(PP_LOCAL, r_total, ibegin);

    // Get into our internal units
    for(int ix = 0;ix < sp->rg_points;ix++) sp->vloc0[ix] /= 2.0;

    // Get the l-value for the local potential if present
    sp->local = xml_tree.get<int>("pseudo.header.<xmlattr>.l_local", -3);

    // Atomic charge density
    std::string PP_RHOATOM = xml_tree.get<std::string>("UPF.PP_RHOATOM");
    sp->atomic_rho = UPF_read_mesh_array(PP_RHOATOM, r_total, ibegin);

    // UPF stores rhoatom * r^2 so rescale
    for(int ix = 0;ix < sp->rg_points;ix++) sp->atomic_rho[ix] = sp->atomic_rho[ix] / (4.0 * PI * sp->r[ix] * sp->r[ix]);

    if(sp->nlccflag) {
        std::string PP_NLCC = xml_tree.get<std::string>("UPF.PP_NLCC");
        sp->rspsco = UPF_read_mesh_array(PP_NLCC, r_total, ibegin);
        for(int ix = 0;ix < sp->rg_points;ix++) sp->rspsco[ix] = sp->rspsco[ix] * 4.0 * PI;
    }

    // Number of atomic orbitals
    sp->num_atomic_waves = xml_tree.get<double>("pseudo.header.<xmlattr>.number_of_wfc", 0);
    sp->num_atomic_waves_m = 0;
    if(sp->num_atomic_waves  > 0) {

        sp->atomic_wave = new double *[MAX_INITWF];
        sp->awave_lig = new double *[MAX_INITWF];
        sp->atomic_wave_l = new int [MAX_INITWF];
        sp->atomic_wave_oc = new double [MAX_INITWF]();

        for(int iwf = 0;iwf < sp->num_atomic_waves;iwf++) {
            // Ugh. UPF format has embedded .s so use / as a separator
            typedef ptree::path_type path;
            std::string chi = "UPF/PP_PSWFC/PP_CHI." + boost::lexical_cast<std::string>(iwf + 1);
            std::string PP_CHI = xml_tree.get<std::string>(path(chi, '/'));
            sp->atomic_wave[iwf] = UPF_read_mesh_array(PP_CHI, r_total, ibegin);

            sp->atomic_wave_l[iwf] = xml_tree.get<int>(path(chi + "/<xmlattr>/l", '/'));
            if(sp->atomic_wave_l[iwf] == 0) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 1;
            if(sp->atomic_wave_l[iwf] == 1) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 3;
            if(sp->atomic_wave_l[iwf] == 2) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 5;
            if(sp->atomic_wave_l[iwf] == 3) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 7;

            //sp->atomic_wave_label[j][0] =
            sp->atomic_wave_oc[iwf] = xml_tree.get<double>(path(chi + "/<xmlattr>/occupation", '/'));

            // UPF stores atomic wavefunctions * r so divide through
            for(int ix = 0;ix < sp->rg_points;ix++) sp->atomic_wave[iwf][ix] /= sp->r[ix];
            sp->awave_lig[iwf] = new double[MAX_LOCAL_LIG]();
            
        }

    }
    else {
       throw RmgFatalException() << "RMG requires pseudopotentials with pseudo atomic wfs. Terminating.\n";  
    }

    // Number of projectors
    sp->nbeta = xml_tree.get<double>("pseudo.header.<xmlattr>.number_of_proj");
    if(sp->nbeta > MAX_NB)
        throw RmgFatalException() << "Pseudopotential has " << sp->nbeta << " projectors but this version of RMG only supports s,p and d proejectors with a limit of " << MAX_NB << " projectors.\n";

    sp->is_ddd_diagonal = true;
    if(sp->nbeta > 0) {

        for(int ip = 0;ip < sp->nbeta;ip++) {
            // Ugh. UPF format has embedded .s so use / as a separator
            typedef ptree::path_type path;
            std::string betapath = "UPF/PP_NONLOCAL/PP_BETA." + boost::lexical_cast<std::string>(ip + 1);
            std::string PP_BETA = xml_tree.get<std::string>(path(betapath, '/'));
            sp->beta[ip] = UPF_read_mesh_array(PP_BETA, r_total, ibegin);

            for(int ix = 0;ix < sp->rg_points;ix++) sp->beta[ip][ix] /= sp->r[ix];
            sp->llbeta[ip] =  xml_tree.get<int>(path(betapath + "/<xmlattr>/angular_momentum", '/'));
            if(sp->llbeta[ip] > ct.max_l) ct.max_l = sp->llbeta[ip];  // For all species
            if(sp->llbeta[ip] > l_max) l_max = sp->llbeta[ip];        // For this species
//               double cutoff_radius = xml_tree.get<int>(betapath + ".<xmlattr>.cutoff_radius");
            
        }

        /*read in the Matrix ddd0(nbeta,nbeta) */
        std::string PP_DIJ = xml_tree.get<std::string>("UPF.PP_NONLOCAL.PP_DIJ");
        double *tmatrix = UPF_str_to_double_array(PP_DIJ, sp->nbeta*sp->nbeta, 0);
        double offd_sum = 0.0;
        for (int j = 0; j < sp->nbeta; j++)
        {
            for (int k = 0; k < sp->nbeta; k++)
            {
                ddd0[j][k] = tmatrix[j*sp->nbeta + k] / 2.0;
                if(j != k) offd_sum += ddd0[j][k]*ddd0[j][k];
            }
        }
        delete [] tmatrix;
        if(offd_sum > 1.0e-20) sp->is_ddd_diagonal = false; 
    }

    sp->nqf=0;
    sp->nlc=0;


    for (int j = 0; j < MAX_NL; j++)
    {
        sp->nhtol[j] = 0;
        sp->nhtom[j] = 0;
        sp->indv[j] = 0;
        sp->nh_l2m[j] = 0;
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

    if (max_nlprojectors > MAX_NL)
        throw RmgFatalException() << "Error in " << __FILE__ << " at line " << __LINE__ << " too many nonlocal projectors.\n"; 

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

    std::string description = xml_tree.get<std::string>("pseudo.header.<xmlattr>.comment", "Pseudopotential");
    sp->description = new char[description.length() + 1]();
    std::strcpy(sp->description, description.c_str());

    // Stuff not present in the UPF format that RMG requires. 
    // We need to find a consistent way of automatically setting these.
    sp->rc = fabs(2.0 * sp->zvalence / sqrt(PI) / sp->vloc0[0]);
    sp->rc = std::max(sp->rc, 0.5);
    sp->rc = std::min(sp->rc, 1.5);
    //sp->rc = 1.0;
    sp->lradius = 8.5;
    sp->gwidth = 8.0;
    sp->rwidth = 15.0; 
    sp->aradius = 9.0;
    sp->acut = 7.0;
    sp->agwidth = 10.0;
    sp->arwidth = 25.0;

    // Leftover initializations
    sp->mill_radius = 9.0;

    // Finally adjust sp->rg_points to skip the high end
    sp->rg_points = iend - ibegin + 1;
    if(pp_buffer) delete [] pp_buffer;
    if(!sp->is_ddd_diagonal) ct.is_ddd_non_diagonal = true;
#endif
}


