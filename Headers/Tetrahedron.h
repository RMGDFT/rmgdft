/*
 *
 * Copyright 2023 The RMG Project Developers. See the COPYRIGHT file 
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


#ifndef RMG_Tetrahedron_H
#define RMG_Tetrahedron_H 1

#include <vector>
#include <cstdint>
#include "rmgtypedefs.h"
#include "BaseGrid.h"
#include "Lattice.h"
#include "Kpoint.h"
#include <boost/multi_array.hpp>

class Tetrahedron
{

    private:

    public:
      Tetrahedron (void);
      ~Tetrahedron (void);
      template <typename KpointType>
      double FillTetra (Kpoint<KpointType> **Kptr);
      double ef, ef_up, ef_down;
};
#endif
