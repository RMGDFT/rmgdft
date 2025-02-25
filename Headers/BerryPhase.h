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


#ifndef RMG_BerryPhase_H
#define RMG_BerryPhase_H 1


class BerryPhase
{
private:
    
public:
    int num_kort, num_kort_pe, num_kpp;
    int BerryPhase_dir;
    std::vector<double> kweight_string;
    BerryPhase();

    ~BerryPhase(void);

    void CalcBP(Kpoint<std::complex<double>> **Kptr);
    void CalcBP(Kpoint<double> **Kptr);
};

#endif
