/*
 *
 * Copyright 2018 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "portability.h"
#include <sys/types.h>
#include <sys/stat.h>
#if !(defined(_WIN32) || defined(_WIN64))
    #include <libgen.h>
#endif
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"

int FileOpenAndCreate(std::string &pathname, int flags, mode_t mode)
{
    // Try to open the file
    int fd = open(pathname.c_str(), flags, mode);

    // Check if it failed because the directory did not exist
    if(fd < 0) 
    {
        char tmpname[MAX_PATH];
        strncpy (tmpname, pathname.c_str(), sizeof(tmpname));

#if !(defined(_WIN32) || defined(_WIN64))
        if (mkdir (dirname (tmpname), S_IRWXU)) return -1;
#else
        char dirname[_MAX_DIR];
        _splitpath(tmpname, NULL, dirname, NULL, NULL);
        if (!_mkdir(dirname)) return -1;
#endif

    }

    // Try to open the file again
    fd = open(pathname.c_str(), flags, mode);
    if(fd < 0) return -1;

    return fd;
}
