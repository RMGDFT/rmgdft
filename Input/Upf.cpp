
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

#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "MapElements.h"
#include "RmgException.h"
#include "transition.h"
#include "InternalPseudo.h"


// Converts a string containg a PP_MESH, PP_RAB, PP_LOCAL, etc into a double array and
// returns a pointer to the array. The calling function needs to have obtained the 
// number of elements expected in the array from the attributes field. The actual
// number of elements placed in the array is count-start where start elements
// are skipped. In other words if there are 100 elements and start=1 the first
// element is skipped and the remaining 99 are placed in the array starting from
// the zeroth position.
double * UPF_str_to_double_array(std::string str, int max_count, int start) {

   std::vector<std::string> strs;
   int count = 0;
   double *array = new double[max_count-start];
   std::string delims = " \t\n";
   boost::trim(str);
   boost::algorithm::split( strs, str, boost::is_any_of(delims), boost::token_compress_on );
 
   std::vector<std::string>::iterator it;
   for (it = strs.begin(); it != strs.end(); ++it) {
       std::string svalue = *it;
       boost::trim(svalue);
       if(count >= start)
           array[count-start] = std::atof(svalue.c_str());
       count++;
       if(count > max_count)
            throw RmgFatalException() << "Problem with UPF pseudopotential in " << __FILE__ << " at line " << __LINE__ << " Too many elements.\n";

       if(count == max_count) return array;
   }

   return array;
}



// Skips the first start elements in array since some UPF pseudopotentials
// have artifacts for small r which can cause us problems
double * UPF_read_mesh_array(std::string str, int count, int start)
{
    return UPF_str_to_double_array(str, count, start);
}


using boost::property_tree::ptree;

// Reads a pseudopotential stored in UPF format into our internal data structures.
void LoadUpf(SPECIES *sp)
{

    int ibegin=0, iend=0;
    ptree upf_tree;
    char *pp_buffer = NULL;
    int pp_buffer_len;
    int max_nlprojectors = 0;
    int l_max;
    std::stringstream ss; 
    double  ddd0[MAX_NL][MAX_NL];  // Used to read in the PP_DIJ
    double qqq[MAX_NL][MAX_NL];    // Used to read in the norms of the augmentation functions (PP_Q)

    std::string Msg;
    if(!std::strcmp(sp->pseudo_filename, "./@Internal") || !strlen(sp->pseudo_filename)) {
 
        std::string pp_string = GetInternalPseudo(&sp->atomic_symbol[0]);
        ss << pp_string;

    }
    else {

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

    }

    // Get the compulsory stuff first
    read_xml(ss, upf_tree);
    std::string PP_INFO = upf_tree.get<std::string>("UPF.PP_INFO"); 
    sp->INFO = new char[PP_INFO.size() + 1]();
    std::strncpy(sp->INFO, PP_INFO.c_str(), PP_INFO.size());
   
    // Atomic symbol, mass, number and zvalence and mesh size
    std::string atomic_symbol = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.element");
    boost::trim(atomic_symbol);
    // Maybe check symbols here
    sp->atomic_mass = GetAtomicMass(atomic_symbol);
    sp->atomic_number = GetAtomicNumber(atomic_symbol);
    sp->zvalence = upf_tree.get<double>("UPF.PP_HEADER.<xmlattr>.z_valence");
    if(sp->zvalence > ct.max_zvalence) ct.max_zvalence = sp->zvalence;

    // Store UPF functional string for later processing
    std::string PP_FUNC = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.functional");
    std::strncpy(sp->functional, PP_FUNC.c_str(), sizeof(sp->functional));


    // Read in the radial mesh and keep values between 1.0e-05 < r < 50.0
    // adjusting mesh size accordingly
    int r_total = upf_tree.get<int>("UPF.PP_HEADER.<xmlattr>.mesh_size");
    std::string PP_R = upf_tree.get<std::string>("UPF.PP_MESH.PP_R");
    double *t_r;
    t_r = UPF_str_to_double_array(PP_R, r_total, 0);
    for(int i = 0;i < r_total;i++) {
        ibegin = i;
        if(t_r[i] > 1.0e-5) break;
    }
    for(int i = ibegin;i < r_total;i++) {
        iend = i;
        if(t_r[i] > 50.0) break;
    }

    sp->rg_points = r_total - ibegin;
    sp->r = new double[sp->rg_points];
    for(int i = ibegin;i < r_total;i++) {
        sp->r[i - ibegin] = t_r[i];
    }
    delete [] t_r;



    l_max = upf_tree.get<int>("UPF.PP_HEADER.<xmlattr>.l_max");
    sp->kkbeta = sp->rg_points;
   
    // Get the type of pseudopotential
    std::string pp_type = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.pseudo_type");
    boost::to_upper(pp_type);
    if(!pp_type.compare(0, 2, "NC")) {
        sp->is_norm_conserving = true;
    }
    else if(!pp_type.compare(0, 2, "SL")) {  // Norm conserving with extra information
        sp->is_norm_conserving = true;
    }
    else if(!pp_type.compare(0, 2, "US")) {
        sp->is_norm_conserving = false;
    }
    else {
        throw RmgFatalException() << "RMG only supports norm conserving and ultrasoft pseudpotentials.\n";
    }


    // Kind of redundant information in the format
    std::string s_is_ultrasoft = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.is_ultrasoft"); 
    boost::to_upper(s_is_ultrasoft);
    if(!s_is_ultrasoft.compare(0,1,"F")) sp->is_norm_conserving = true;
    if(!s_is_ultrasoft.compare(0,5,"FALSE")) sp->is_norm_conserving = true;
    if(!s_is_ultrasoft.compare(0,1,"T")) sp->is_norm_conserving = false;
    if(!s_is_ultrasoft.compare(0,4,"TRUE")) sp->is_norm_conserving = false;

    // Check for full relativistic and thrown an error if found
    std::string s_is_relativistic = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.relativistic");
    boost::to_upper(s_is_relativistic);
    if(!s_is_relativistic.compare(0,4,"FULL")) {
        throw RmgFatalException() << "RMG does not support fully relativistic pseudopotentials. Terminating.\n";
    }

    // Core correction flag
    std::string s_core_correction = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.core_correction");
    if(!s_core_correction.compare(0,1,"F")) sp->nlccflag = false;
    if(!s_core_correction.compare(0,5,"FALSE")) sp->nlccflag = false;
    if(!s_core_correction.compare(0,1,"T")) sp->nlccflag = true;
    if(!s_core_correction.compare(0,4,"TRUE")) sp->nlccflag = true;

    // Determine log mesh parameters directly from the mesh
    sp->aa = (sp->r[0] * sp->r[0]) / (sp->r[1] - 2 * sp->r[0]);
    sp->bb = log (sp->r[1] / sp->r[0] - 1);

    // Read in rab and convert it into a C style array
    std::string PP_RAB = upf_tree.get<std::string>("UPF.PP_MESH.PP_RAB");
    sp->rab = UPF_read_mesh_array(PP_RAB, r_total, ibegin);

    // Local potential
    std::string PP_LOCAL = upf_tree.get<std::string>("UPF.PP_LOCAL");
    sp->vloc0 = UPF_read_mesh_array(PP_LOCAL, r_total, ibegin);

    // Get into our internal units
    for(int ix = 0;ix < sp->rg_points;ix++) sp->vloc0[ix] /= 2.0;

    // Get the l-value for the local potential if present
    sp->local = upf_tree.get<int>("UPF.PP_HEADER.<xmlattr>.l_local", -3);

    // Atomic charge density
    std::string PP_RHOATOM = upf_tree.get<std::string>("UPF.PP_RHOATOM");
    sp->atomic_rho = UPF_read_mesh_array(PP_RHOATOM, r_total, ibegin);

    // UPF stores rhoatom * r^2 so rescale
    for(int ix = 0;ix < sp->rg_points;ix++) sp->atomic_rho[ix] = sp->atomic_rho[ix] / (4.0 * PI * sp->r[ix] * sp->r[ix]);

    if(sp->nlccflag) {
        std::string PP_NLCC = upf_tree.get<std::string>("UPF.PP_NLCC");
        sp->rspsco = UPF_read_mesh_array(PP_NLCC, r_total, ibegin);
        for(int ix = 0;ix < sp->rg_points;ix++) sp->rspsco[ix] = sp->rspsco[ix] * 4.0 * PI;
    }

    // Number of atomic orbitals
    sp->num_atomic_waves = upf_tree.get<double>("UPF.PP_HEADER.<xmlattr>.number_of_wfc", 0);
    sp->num_atomic_waves_m = 0;
    if(sp->num_atomic_waves  > 0) {

        sp->atomic_wave = new double *[sp->num_atomic_waves];
        sp->awave_lig = new double *[sp->num_atomic_waves];
        sp->atomic_wave_l = new int [sp->num_atomic_waves];
        sp->atomic_wave_oc = new double [sp->num_atomic_waves]();

        for(int iwf = 0;iwf < sp->num_atomic_waves;iwf++) {
            // Ugh. UPF format has embedded .s so use / as a separator
            typedef ptree::path_type path;
            std::string chi = "UPF/PP_PSWFC/PP_CHI." + boost::lexical_cast<std::string>(iwf + 1);
            std::string PP_CHI = upf_tree.get<std::string>(path(chi, '/'));
            sp->atomic_wave[iwf] = UPF_read_mesh_array(PP_CHI, r_total, ibegin);

            sp->atomic_wave_l[iwf] = upf_tree.get<int>(path(chi + "/<xmlattr>/l", '/'));
            if(sp->atomic_wave_l[iwf] == 0) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 1;
            if(sp->atomic_wave_l[iwf] == 1) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 3;
            if(sp->atomic_wave_l[iwf] == 2) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 5;
            if(sp->atomic_wave_l[iwf] == 3) sp->num_atomic_waves_m = sp->num_atomic_waves_m + 7;

            //sp->atomic_wave_label[j][0] =
            sp->atomic_wave_oc[iwf] = upf_tree.get<double>(path(chi + "/<xmlattr>/occupation", '/'));

            // UPF stores atomic wavefunctions * r so divide through
            for(int ix = 0;ix < sp->rg_points;ix++) sp->atomic_wave[iwf][ix] /= sp->r[ix];
            sp->awave_lig[iwf] = new double[MAX_LOCAL_LIG]();
            
        }

    }
    else {
       throw RmgFatalException() << "RMG requires pseudopotentials with pseudo atomic wfs. Terminating.\n";  
    }

    // Number of projectors
    sp->nbeta = upf_tree.get<double>("UPF.PP_HEADER.<xmlattr>.number_of_proj");
    if(sp->nbeta > MAX_NB)
        throw RmgFatalException() << "Pseudopotential has " << sp->nbeta << " projectors but this version of RMG only supports s,p and d proejectors with a limit of " << MAX_NB << " projectors.\n";

    if(sp->nbeta > 0) {

        for(int ip = 0;ip < sp->nbeta;ip++) {
            // Ugh. UPF format has embedded .s so use / as a separator
            typedef ptree::path_type path;
            std::string betapath = "UPF/PP_NONLOCAL/PP_BETA." + boost::lexical_cast<std::string>(ip + 1);
            std::string PP_BETA = upf_tree.get<std::string>(path(betapath, '/'));
            sp->beta[ip] = UPF_read_mesh_array(PP_BETA, r_total, ibegin);

            for(int ix = 0;ix < sp->rg_points;ix++) sp->beta[ip][ix] /= sp->r[ix];
            sp->llbeta[ip] =  upf_tree.get<int>(path(betapath + "/<xmlattr>/angular_momentum", '/'));
            if(sp->llbeta[ip] > ct.max_l) ct.max_l = sp->llbeta[ip];  // For all species
            if(sp->llbeta[ip] > l_max) l_max = sp->llbeta[ip];        // For this species
//               double cutoff_radius = upf_tree.get<int>(betapath + ".<xmlattr>.cutoff_radius");
            
        }

        /*read in the Matrix ddd0(nbeta,nbeta) */
        std::string PP_DIJ = upf_tree.get<std::string>("UPF.PP_NONLOCAL.PP_DIJ");
        double *tmatrix = UPF_str_to_double_array(PP_DIJ, sp->nbeta*sp->nbeta, 0);
        for (int j = 0; j < sp->nbeta; j++)
        {
            for (int k = 0; k < sp->nbeta; k++)
            {
                ddd0[j][k] = tmatrix[j*sp->nbeta + k] / 2.0;
            }
        }
        delete [] tmatrix;
      
    }

    sp->nqf=0;
    sp->nlc=0;
    // Charge augmentation for US
    if(!sp->is_norm_conserving) {

       std::string q_with_l = upf_tree.get<std::string>("UPF.PP_NONLOCAL.PP_AUGMENTATION.<xmlattr>.q_with_l");
       boost::to_upper(q_with_l);
       if(!q_with_l.compare(0,1,"T") || !q_with_l.compare(0,4,"TRUE")) {
           throw RmgFatalException() << "RMG does not support pseudopotentials with an angular momentum dependence of the augmentation charges. Terminating.\n";  
       }



       sp->nqf = upf_tree.get<int>("UPF.PP_NONLOCAL.PP_AUGMENTATION.<xmlattr>.nqf", 0);
       sp->nlc = upf_tree.get<int>("UPF.PP_NONLOCAL.PP_AUGMENTATION.<xmlattr>.nqlc", 2*l_max + 1);
       int tnlc = (sp->nbeta * (sp->nbeta + 1)) / 2;

       // read in the rinner for Q_L(r) function */
       if(sp->nqf > 0) {
           std::string PP_RINNER = upf_tree.get<std::string>("UPF.PP_NONLOCAL.PP_AUGMENTATION.PP_RINNER");
           sp->rinner = new double[tnlc * sp->nlc * sp->nqf]();
           double *tmatrix = UPF_str_to_double_array(PP_RINNER, sp->nlc, 0);
           for(int i = 0;i < sp->nlc;i++) sp->rinner[i] = tmatrix[i];
           delete [] tmatrix;
       }

       /*read in the Matrix qqq(nbeta,nbeta) */
       std::string PP_Q = upf_tree.get<std::string>("UPF.PP_NONLOCAL.PP_AUGMENTATION.PP_Q");
       double *tmatrix = UPF_str_to_double_array(PP_Q, sp->nbeta*sp->nbeta, 0);
       for (int j = 0; j < sp->nbeta; j++)
       {
           for (int k = 0; k < sp->nbeta; k++)
           {
               qqq[j][k] = tmatrix[j*sp->nbeta + k];
           }
       }
       delete [] tmatrix;

       // Read in PP_QFCOEFF
       if(sp->nqf > 0) {
           // UPF format calls for full matrix so store in tmatrix
           std::string PP_QFCOEFF = upf_tree.get<std::string>("UPF.PP_NONLOCAL.PP_AUGMENTATION.PP_QFCOEF");
           double *tmatrix1 = UPF_str_to_double_array(PP_QFCOEFF, sp->nqf*sp->nlc*sp->nbeta*sp->nbeta, 0);
           // We only use upper half of the symmetric matrix
           sp->qfcoef = new double[tnlc * sp->nlc * sp->nqf];

           for (int idx1 = 0; idx1 < sp->nbeta; idx1++)
           {
               for (int idx2 = idx1; idx2 < sp->nbeta; idx2++)
               {
                   int nmb = idx2 * (idx2 + 1) / 2 + idx1;
                   for (int idx3 = 0; idx3 < sp->nlc; idx3++)
                   {
                       for (int idx4 = 0; idx4 < sp->nqf; idx4++)
                       {
                           int indx = nmb * sp->nlc * sp->nqf + idx3 * sp->nqf + idx4;
                           sp->qfcoef[indx] = tmatrix1[sp->nqf*sp->nlc*sp->nbeta*idx1 + sp->nqf*sp->nlc*idx2 + sp->nqf*idx3 + idx4];
                       }
                   }
               }
           }

           delete [] tmatrix1;  

           // PP_Q.i.j
           sp->qnm = new double[tnlc * MAX_RGRID]();
           for(int i = 0;i < sp->nbeta;i++) {
               for(int j = i;j < sp->nbeta;j++) {
                   typedef ptree::path_type path;
                   std::string pp_qij = "UPF/PP_NONLOCAL/PP_AUGMENTATION/PP_QIJ." + boost::lexical_cast<std::string>(i + 1);
                   pp_qij = pp_qij + "." + boost::lexical_cast<std::string>(j + 1);
                   std::string PP_Qij = upf_tree.get<std::string>(path(pp_qij, '/'));
                   double *tmatrix2 = UPF_read_mesh_array(PP_Qij, r_total, ibegin);


                   int nmb = j * (j + 1) / 2 + i;
                   for (int idx = 0; idx < sp->kkbeta; idx++)
                   {
                       sp->qnm[nmb * MAX_RGRID + idx] = tmatrix2[idx];
                   }

                   delete [] tmatrix2;
               }
           }
       }

    }

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
    
    std::string description = upf_tree.get<std::string>("UPF.PP_HEADER.<xmlattr>.comment", "Pseudopotential");
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

}

// C binding
extern "C" void LoadUpf_C(SPECIES *sp)
{
    LoadUpf(sp);
}

