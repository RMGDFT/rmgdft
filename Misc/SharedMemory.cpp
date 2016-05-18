/*
 *
 * Copyright 2014 The RMG Project Developers. See the COPYRIGHT file 
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

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/creation_tags.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <cstring>
#include <cstdlib>
#include <string>
#include <unordered_set>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "RmgShm.h"


// This collection of functions is used to manage shared memory objects for use
// by MPI process's running on the same node.


// Stores names of all shared blocks that we have allocated
static std::unordered_set<std::string> shared_segments;


void *AllocSharedMemorySegment(char *name, int size)
{
    using namespace boost::interprocess;

    if(pct.is_local_master) {
        shared_memory_object segment(open_or_create, name, read_write);
        segment.truncate(size);
    }

    MPI_Barrier(pct.local_comm);
    shared_memory_object segment(open_only, name, read_write);
    mapped_region region(segment, read_write);
    void *rptr = static_cast<char*>(region.get_address());
    shared_segments.emplace(name);
    return rptr;
}

void FreeSharedMemory(char *name)
{
    using namespace boost::interprocess;
    shared_memory_object::remove(name);
    shared_segments.erase(name);
}

void FreeAllSharedMemory(void)
{
    using namespace boost::interprocess;
    for(auto it = shared_segments.cbegin();it != shared_segments.cend();++it) {
        shared_memory_object::remove(it->c_str());
    }
    shared_segments.erase(shared_segments.begin());
}
