#include "portability.h"
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <string>
#include <boost/algorithm/string.hpp>
#include "transition.h"
#include "RmgException.h"


// Used to increment the logfile number

int FilenameIncrement(char *pathname)
{

    int lognum = 0;

    boost::filesystem::path current_path(pathname);
    std::string dirname  = current_path.parent_path().string(); 
    std::string basename = boost::filesystem::basename(pathname);
    if(!dirname.length()) dirname = dirname + "./";
    // Does parent path exist?
    if( boost::filesystem::exists(dirname)) {

        // yes so check if it's a file 
        if(!boost::filesystem::is_directory(dirname)) {
            throw RmgFatalException() << "Found " << dirname << "  that is not a directory in " __FILE__ << " at line " << __LINE__ << ".\n";
        }

    }
    else {

        // no so need to make it
        if(!boost::filesystem::create_directory(dirname)) {
            throw RmgFatalException() << "Unable to create logfile directory " << dirname << " in " << __FILE__ << " at line " << __LINE__ << ".\n";
        }

    }

    
    char lognum_str[8]; 
    snprintf(lognum_str, 3, "%02d", lognum);
    std::string nextname = std::string(pathname) + "." + lognum_str + ".log";
    while (boost::filesystem::exists(nextname))
    {
        if (++lognum > 99)
            throw RmgFatalException() << "You have over 100 logfiles, you need to think of a better job naming scenario!" << "\n";
        nextname.erase();
        snprintf(lognum_str, 4, "%02d", lognum);
        nextname = std::string(pathname) + "." + lognum_str + ".log";
    }

    // return lognum
    return lognum;

}
