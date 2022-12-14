/*
 *
 * Copyright 2017 The RMG Project Developers. See the COPYRIGHT file 
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


#ifndef RMG_Elpa_H
#define RMG_Elpa_H 1

#include <mpi.h>
#include "Scalapack.h"
#if USE_ELPA
//#include "elpa/elpa.h"
#endif

class Elpa: public Scalapack {

public:

  Elpa(int ngroups, int thisimg, int images_per_node, int N, int NB, int last, MPI_Comm rootcomm):Scalapack(ngroups, thisimg, images_per_node, N, NB, last, rootcomm){}
  void Init(void);
  int GetElpaCommRows(void);
  int GetElpaCommCols(void);
  void generalized_eigenvectors(double *a, double *b, double *ev, double *q, int is_already_decomposed, int *error);
  void generalized_eigenvectors(std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, int is_already_decomposed, int *error);

  ~Elpa(void);

private:
   int elpa_comm_rows, elpa_comm_cols;
   void *handle;
//   elpa_autotune_t autotune_handle;  

};

#endif
