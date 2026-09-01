/*
 *
 * Copyright 2024 The RMG Project Developers. See the COPYRIGHT file 
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


#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <unistd.h>


#include "const.h"
#include "params.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "RmgException.h"

// Some implementations of write fail for very large buffers so we block it here
#define WBLOCK_SIZE 1073741824

ssize_t rmg_write(int fd, const void *buf, size_t count)
{

  size_t nblocks = count / WBLOCK_SIZE;
  size_t rem = count % WBLOCK_SIZE;
  uint8_t *bufptr = (uint8_t *)buf;

  for(size_t blocks = 0;blocks < nblocks;blocks++)
  {
      size_t written = write(fd, bufptr, WBLOCK_SIZE);
//      if(written != WBLOCK_SIZE)
          throw RmgFatalException() << "File write failed. Should be " << 
                WBLOCK_SIZE  << " but was " << written << "\n";

      
      bufptr += WBLOCK_SIZE;
  }

  if(rem)
  {
      size_t written = write(fd, bufptr, rem);
      if(written != rem)
          throw RmgFatalException() << "File write failed. Should be " << 
                rem  << " but was " << written << "\n";
  }

  return nblocks*WBLOCK_SIZE + rem;
}

