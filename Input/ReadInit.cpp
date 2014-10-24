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
#include <boost/filesystem.hpp>
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
#include "InputOpts.h"
#include "grid.h"


void ReadInit(char *meta, CONTROL& lc, PE_CONTROL& pelc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    if(!boost::filesystem::exists(meta) )
    {
        std::cout << "\n using default path and input"; 
        std::strncpy(pelc.image_path[0], "./", sizeof(pelc.image_path[0]));
        std::strncpy(pelc.image_input[0], "input", sizeof(pelc.image_input[0]));
        pelc.images = 1;
        pelc.image_npes[0] = pelc.total_npes;
        lc.images_per_node = 1;
        lc.spin_flag = 0;
        pelc.pe_kpoint = 1;
        return;
    }

    RmgInputFile If(meta, InputMap);

    /* do spin polarized calculation? */
    If.RegisterInputKey("spin_polarization", &lc.spin_flag, false,
                        "Spin polarized calculation.");

    If.RegisterInputKey("num_images", &pelc.images, 1, INT_MAX, 1,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "number of images to use",
                     "num_images must be a positive integer.");

    If.RegisterInputKey("image_per_node", &ct.images_per_node, 1, INT_MAX, 1,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "number of images per node",
                     "image_per_node must be a positive integer.");

    If.RegisterInputKey("pegroup_for_kpoint", &pelc.pe_kpoint, 1, INT_MAX, 1,
                     CHECK_AND_TERMINATE, OPTIONAL,
                     "pegroup_for_kpoint",
                     "pegroup_for_kpoint must be a positive integer.");


    std::string ImageInfo;
    If.RegisterInputKey("image_infos", &ImageInfo, "",
                     CHECK_AND_FIX, OPTIONAL,
                     "Image information.\n",
                     "You must provide information for at least 1 image.");

    If.LoadInputKeys();

    int num_image = 0;
    int tot_pe = 0;
    std::vector<std::string> Images;
    std::string line_delims = "\r\n^";
    std::string whitespace_delims = " \t";

    if(ImageInfo.size() > 0) {
        boost::algorithm::split( Images, ImageInfo, boost::is_any_of(line_delims), boost::token_compress_on );
        for (auto it = Images.begin(); it != Images.end(); ++it) {
            std::string image_line = *it;
            boost::trim_if(image_line, boost::algorithm::is_any_of(" \t"));
            std::vector<std::string> fields;
            boost::algorithm::split( fields, image_line, boost::is_any_of(whitespace_delims), boost::token_compress_on );
            std::strncpy(pelc.image_path[num_image], fields[0].c_str(), sizeof(pelc.image_path[num_image]));
            std::strncpy(pelc.image_input[num_image], fields[1].c_str(), sizeof(pelc.image_input[num_image]));
            pelc.image_npes[num_image] = std::atoi(fields[2].c_str());
            tot_pe += pelc.image_npes[num_image];
            num_image++;
        }
    }

    if(num_image == 0) 
    {
        tot_pe = pelc.total_npes;
        std::strcpy(pelc.image_path[0],  "./");
        std::strcpy(pelc.image_input[0] , meta);
        pelc.image_npes[0] = pelc.total_npes;
        num_image = 1;
    }
    
    
    if(tot_pe != pelc.total_npes) 
    {
        throw RmgFatalException() << "\n require " << tot_pe << " npes != from job/pe_kpoint " << pelc.total_npes;
    }
    if(num_image != pelc.images) 
    {
        throw RmgFatalException() << "\n number of image and image info not consistent";
    }
    if(pelc.images > MAX_IMGS)
    {
        throw RmgFatalException() << "\n number of image " << pelc.images << " > " << MAX_IMGS << " MAX_IMGS";
    }

    if(pelc.images % ct.images_per_node != 0)
    {
        throw RmgFatalException() << "\nyou are crazy to use the multiple image per node feature"; 
    }
    if(ct.images_per_node >1 )
    {
        for(int i = 1; i < pelc.images; i++) 
            if(pelc.image_npes[i] != pelc.image_npes[0])
            {
                throw RmgFatalException() << "\n image " << pelc.image_npes[i] << " has different NPES " << pelc.image_npes[0];
            }
    }
}

