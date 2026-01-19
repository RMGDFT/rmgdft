/*  
 *  
 * Copyright 2026 The RMG Project Developers. See the COPYRIGHT file 
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
        

#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include "rmg_error.h"
#include "rmg_mangling.h"
#include "params.h"
#include "pe_control.h"
#include "transition.h"

#define         errore          RMG_FC_GLOBAL(errore, ERRORE)


// Here for fortran routines.
extern "C" void errore(char *where, char *message, int ierr, int where_len, int message_len)
{
  char tbuf[1000];
  memset(tbuf, 0, sizeof(tbuf));

  if(((where_len + message_len) > (int)sizeof(tbuf)) || (where_len < 0) || (message_len < 0)) {
     std::cout << "Unknown issue printing error message from fortran routines" << std::endl;
     printf("Unknown issue printing error message from fortran routines\n");fflush(NULL);
     rmg::printlog("Unknown issue printing error message from fortran routines\n");fflush(NULL);
     raise(SIGTERM);
  }

  strncpy(tbuf, message, sizeof(tbuf)-1);
  int iz = message_len;
  if(iz >= (int)sizeof(tbuf)) iz = sizeof(tbuf) - 1;
  tbuf[iz] = 0;
  printf("%s in ", tbuf);
  strncpy(tbuf, where, sizeof(tbuf)-1);
  iz = where_len;
  if(iz >= (int)sizeof(tbuf)) iz = sizeof(tbuf) - 1;
  tbuf[iz] = 0;
  printf("%s\n", tbuf);
  rmg::printlog("%s\n", tbuf);


  fflush (NULL);
#if (defined(_WIN32) || defined(_WIN64))
    Sleep(2000);
#else
    sleep (2);
#endif
    MPI_Abort( MPI_COMM_WORLD, 0 );

}
