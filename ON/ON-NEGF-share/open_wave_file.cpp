/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/


#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <filesystem>
#include "main.h"
#include "prototypes_on.h"
#include "transition.h"



/*This opens s file for writing, returns a file handle
 * If opening the file fails, PE 0 tries to create a directory, since it is possible that
 * the reason for failure is that the directory does not exist*/

int open_wave_file (const char *filename)
{

    std::filesystem::path p(filename);
    std::filesystem::create_directories(p.parent_path());

    int amode;
    int fhand;

    amode = S_IREAD | S_IWRITE;

    fhand = open (filename, O_CREAT | O_TRUNC | O_RDWR, amode);

    if (fhand < 0)
    {
        rmg::printlog("\n cannot open file %s", filename);
        rmg::error("cannot open file ");
    }

    return fhand;

}
