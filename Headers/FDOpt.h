/*
 *
 * Copyright 2022 The RMG Project Developers. See the COPYRIGHT file 
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

#ifndef RMG_FDOpt_H
#define RMG_FDOpt_H 1

#include <vector>

class FDOpt
{
private:
    static int pbasis;
    static int sbasis;
    static int order;
    static int num_orb;
    static int num_coeff;
    static std::vector<double> ke_fft;
    static std::vector<double> ke_fd;
    static std::vector<double> occ_weight;
    static std::vector<double> coeff;
    static std::vector<double> coeff_grad;
    static double *orbitals;
    static double *orbitals_b;
    static double *psi_psin;
    static double *work;
    static double kvec[3];

public:
    FDOpt(void);
    ~FDOpt(void);
    double Optimize(void);
    static double evaluate(void *,
             const double *x,
             double *g,
             const int n);
    static double stepbound(void *,
             const double *x,
             const double *g,
             const int n);
    static int progress(void *instance,
                         const double *x,
                         const double *g,
                         const double fx,
                         const double xnorm,
                         const double gnorm,
                         const double step,
                         int n,
                         int k,
                         int ls);

};

#endif
