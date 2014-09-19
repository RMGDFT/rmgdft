
#include <exception>
#include <iostream>
#include <vector>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "make_conf.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"


// Converts a string containg a PP_MESH, PP_RAB, PP_LOCAL, etc into a double array and
// returns a pointer to the array. The calling function needs to have obtained the 
// number of elements expected in the array from the attributes field.
double * UPF_str_to_double_array(std::string str, int max_count) {

   std::vector<std::string> strs;
   int count = 0;
   double *array = new double[max_count];
   std::string delims = " \t\n";
   boost::trim(str);
   boost::algorithm::split( strs, str, boost::is_any_of(delims), boost::token_compress_on );
 
   std::vector<std::string>::iterator it;
   for (it = strs.begin(); it != strs.end(); ++it) {
       std::string svalue = *it;
       array[count] = std::atof(svalue.c_str());
       count++;
       if(count > max_count)
           rmg_error_handler(__FILE__,__LINE__,"Problem with UPF pseudopotential. Too many elements.\n");
   }

   return array;
}


using boost::property_tree::ptree;

// Reads a pseudopotential stored in UPF format into our internal data structures.
void LoadUpf(char *file, SPECIES *sp)
{

 ptree upf_tree;
 std::string PP_INFO;
 std::string upf_file = file;
 
 read_xml(upf_file, upf_tree);
 PP_INFO = upf_tree.get<std::string>("UPF.PP_INFO"); 
   
 std::cout << PP_INFO << std::endl; 

 // Attributes of the mesh
 double PP_MESH_dx = upf_tree.get<double>("UPF.PP_MESH.<xmlattr>.dx");
 std::cout << "PP_MESH.dx    =  " << PP_MESH_dx << std::endl; 

 int PP_MESH_mesh = upf_tree.get<double>("UPF.PP_MESH.<xmlattr>.mesh");
 std::cout << "PP_MESH.mesh  =  " << PP_MESH_mesh << std::endl; 
 
 double PP_MESH_xmin = upf_tree.get<double>("UPF.PP_MESH.<xmlattr>.xmin");
 std::cout << "PP_MESH.xmin  =  " << PP_MESH_xmin << std::endl; 

 double PP_MESH_rmax = upf_tree.get<double>("UPF.PP_MESH.<xmlattr>.rmax");
 std::cout << "PP_MESH.rmax  =  " << PP_MESH_rmax << std::endl; 

 double PP_MESH_zmesh = upf_tree.get<double>("UPF.PP_MESH.<xmlattr>.zmesh");
 std::cout << "PP_MESH.zmesh  =  " << PP_MESH_zmesh << std::endl; 


 // Read in the radial mesh and convert it into a C style array
 std::string PP_R = upf_tree.get<std::string>("UPF.PP_MESH.PP_R");

 double *rmesh = UPF_str_to_double_array(PP_R, PP_MESH_mesh);
 for(int ix=0;ix < PP_MESH_mesh;ix++) {
     std::cout << rmesh[ix] << std::endl;
 }

}

// C binding
extern "C" void LoadUpf_C(char *file, SPECIES *sp)
{
    LoadUpf_C(file, sp);
}

