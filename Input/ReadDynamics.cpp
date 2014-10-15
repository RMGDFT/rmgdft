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

    For detailed documentation on how the input system works look at
    Input/ReadCommon.cpp. This function is used to read specific
    information related to molecular dynamics and relaxations.

**********************************************************************/

namespace Ri = RmgInput;

void ReadDymamics(char *cfile, CONTROL& lc, std::unordered_map<std::string, InputKey *>& InputMap)
{

    RmgInputFile If(cfile, InputMap);

}
