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


#ifndef RMG_Dos_H
#define RMG_Dos_H 1

#include <string>
#include "Lattice.h"
#include <boost/multi_array.hpp>

// Ref: PRB89, 094515 (2014) and QE code
class Dos {

private:
    // Lattice object
    Lattice &L;

    bool tetraflag;
    int num_tetra;
    int nk_per_tetra;
    double gaus_broad;
    boost::multi_array<int, 2> tetra;

    void tetra_init(int kmesh[3], int kshift[3]);
    double Pij[4][20];
 
    double tot_dos_tetra(int nk, int nband, std::vector<double> &eigs, double Ef, std::vector<double> &e_list, std::vector<double> &dos_t);
    double tot_dos_gauss(int nk, int nband, std::vector<double> &eigs, double energy);
    void set_Pij();

public:
    Dos (int kpoint_mesh[3], int kpoint_is_shift[3], Lattice &L, double gaus_broad);
    void tot_dos(int nk, int nband, std::vector<double> eigs, double Ef);

    ~Dos(void);

};

#endif
