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
#include <unordered_map>
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

using namespace boost::interprocess;

// This collection of functions is used to manage shared memory objects for use
// by MPI process's running on the same node.


// Stores names of all shared blocks that we have allocated
static std::unordered_map<std::string, mapped_region *> shared_segments;


// Called at initialization to cleanup any old shared memory segments
extern "C" void InitSharedMemory(void)
{

}

extern "C" void *AllocSharedMemorySegment(char *name, int size)
{
    mapped_region *region;
    void *rptr;
    // Maybe not the best way to do this but should help avoid a name collision
    std::string newname(name);
    newname = newname + "435lkj91kmSpv10rwiq";

    if(pct.is_local_master) {
        shared_memory_object::remove(newname.c_str());
        shared_memory_object segment(open_or_create, newname.c_str(), read_write);
        segment.truncate(size);
    }

    MPI_Barrier(pct.local_comm);
    try {
        shared_memory_object segment(open_only, newname.c_str(), read_write);
        region = new mapped_region(segment, read_write);
        rptr = static_cast<char*>(region->get_address());
    }
    catch(interprocess_exception &ex) {
        return NULL;
    }
    if(size != (int)region->get_size()) {
        delete region;
        return NULL;
    }
    shared_segments.emplace(newname, region);
    return rptr;
}

extern "C" void FreeSharedMemory(char *name)
{
    std::string newname(name);
    newname = newname + "435lkj91kmSpv10rwiq";
    shared_memory_object::remove(newname.c_str());
    shared_segments.erase(newname);
}

extern "C" void FreeAllSharedMemory(void)
{

    if(shared_segments.size() == 0) return;

    for(auto it = shared_segments.cbegin();it != shared_segments.cend();++it) {
        std::string name = it->first;
        mapped_region *region = it->second;
        delete region;
        if(pct.is_local_master) shared_memory_object::remove(name.c_str());
    }
    shared_segments.erase(shared_segments.begin());
}
