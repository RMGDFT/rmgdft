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


#ifndef RMG_Stress_H
#define RMG_Stress_H 1


//#include "Kpoint.h"
//#include "Lattice.h"
//#include "BaseGrid.h"
//#include "species.h"
//#include "Pw.h"



template <typename T> class Stress 
{
private:
 //   BaseGrid BG;
    
public:
    double stress_tensor[9];
    Stress(Kpoint<T> **Kpin, Lattice &L, BaseGrid &BG, Pw &pwaves,
            std::vector<ION> &atoms, std::vector<SPECIES> &species, double Exc, double *vxc, double *rho, double *rhocore, double *vtot);

    void Ewald_term(std::vector<ION> &atoms, std::vector<SPECIES> &species, Lattice &L, Pw &pwaves);
    void Kinetic_term(Kpoint<T> **Kpin, BaseGrid &BG, Lattice &L);
    void Local_term(std::vector<ION> &atoms, std::vector<SPECIES> &species, double *rho, Pw &pwaves);
    void Local_term1(double *rho, double *vnuc);
    void Hartree_term(double *rho, Pw &pwaves);
    void NonLocal_term(Kpoint<T> **Kpin, std::vector<ION> &atoms, std::vector<SPECIES> &species);
    void NonLocalQfunc_term(Kpoint<T> **Kpin, std::vector<ION> &atoms, std::vector<SPECIES> &species, double *vtot);
    void Exc_term(double Exc, double *vxc, double *rho);
    void Exc_gradcorr(double Exc, double *vxc, double *rho, double *rhocore);
    void Exc_Nlcc(double *vxc, double *rhocore);
    ~Stress(void);

};

#endif
