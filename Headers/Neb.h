/*
 *
 * Copyright 2019 The RMG Project Developers. See the COPYRIGHT file 
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


#ifndef RMG_Neb_H
#define RMG_Neb_H 1


#include <string>
#include <vector>
#include <set>
#include <complex>
#include <mutex>
#include "BaseGrid.h"
#include "ION.h"
#include "Kpoint.h"

template <typename T> class Neb {

private:
    // BaseGrid class (distributed)
    BaseGrid &BG;
    std::vector<ION> Atoms_initial;
    std::vector<ION> Atoms_final;
    int num_images;
    int max_steps;
    double totale_initial;
    double totale_final;
    double *L_coor;
    double *S_coor;
    double *R_coor;
    double *totale;
    double *all_frc;
    double *path_length;


public:
    Neb( BaseGrid &BG, int num_images, int max_steps, std::string input_initial, 
            std::string input_final, double totale_initial, double total_final); 

    ~Neb(void);

    void relax (double * vxc, double * vh, double * vnuc, double * rho, 
            double * rho_oppo, double * rhocore, double * rhoc, Kpoint<T> **Kptr);
    void climb_image();

};

#endif
