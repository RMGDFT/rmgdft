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

#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "main.h"


// Creates an mmap array of the given length associated with
// the file descriptor specified by fd. If the file already
// exists it's length is adjusted to length.
void *CreateMmapArray(int &fd, size_t length)
{
    ftruncate(fd, length);
    lseek(fd, length-1, SEEK_SET);
    write(fd, "", 1);
    return mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
}
