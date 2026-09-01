#include <filesystem>
#include <string>
#include <boost/algorithm/string.hpp>
#include "transition.h"
#include "RmgException.h"


// Used to increment the logfile number

int FilenameIncrement(char *pathname)
{

    int lognum = 0;

    std::filesystem::path current_path(pathname);
    std::string dirname  = current_path.parent_path().string(); 
    std::string basename = current_path.stem().string();
    if(!dirname.length()) dirname = dirname + "./";
    // Does parent path exist?
    if( std::filesystem::exists(dirname)) {

        // yes so check if it's a file 
        if(!std::filesystem::is_directory(dirname)) {
            throw RmgFatalException() << "Found " << dirname << "  that is not a directory in " __FILE__ << " at line " << __LINE__ << ".\n";
        }

    }
    else {

        // no so need to make it
        if(!std::filesystem::create_directory(dirname)) {
            throw RmgFatalException() << "Unable to create logfile directory " << dirname << " in " << __FILE__ << " at line " << __LINE__ << ".\n";
        }

    }

    
    char lognum_str[8]; 
    snprintf(lognum_str, 3, "%02d", lognum);
    std::string nextname = std::string(pathname) + "." + lognum_str + ".log";
    while (std::filesystem::exists(nextname))
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
