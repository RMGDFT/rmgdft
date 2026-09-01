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

#include <complex>
#include <cfloat>
#include <cmath>
#include <cfenv>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "Subdiag.h"
#include "ErrorFuncs.h"
#include "transition.h"


// Sets up the parallel decomposition across processing elements
// A processing element in this case is not necessarily a single
// MPI process but could be a group of MPI processes (e.g. scalapack)
// 
// n = size of square matrix to be diagonalized
// FS_NPES = number of processing elements
// THISPE = index of this PE starting from 0
// fs_eigstart, fs_eigstop, fs_eigcounts = arrays used in MPI calls
// eigstart = index of starting eigenvector this PE is responsible for starting from 0
// eigstop = index of last eigenvector this PE is responsible for - 1
// n_start = start of window
// n_win = width of window
//
//
//  ASCII diagram below
//
//     |      |           |             |            |                   |
//     0      n_start     eig_start     eig_stop     n_start+n_win       n               
//
void FoldedSpectrumSetup(int n, int FS_NPES, int THISPE, 
                         int& eig_start, int& eig_stop, int& eig_step,
                         int& n_start, int& n_win,
                         int *fs_eigstart, int *fs_eigstop, int *fs_eigcounts, int blocksize)
{

    for(int idx = 0;idx < FS_NPES;idx++) {
        double t1 = (double)n;
        t1 = t1 / ((double)FS_NPES);
        double t2 = t1 * (double)idx;
        fs_eigstart[idx] = (int)rint(t2);
        fs_eigstop[idx] = (int)rint(t1 + t2);
        fs_eigcounts[idx] = fs_eigstop[idx] - fs_eigstart[idx];
        fs_eigstart[idx] *= n;
        fs_eigstop[idx] *= n;
        fs_eigcounts[idx] *= n;

    }


    // Folded spectrum method is parallelized over PE's. Each PE gets assigned
    // a subset of the eigenvectors.
    double t1 = (double)n;
    t1 = t1 / ((double)FS_NPES);
    double t2 = t1 * (double)THISPE;
    eig_start = (int)rint(t2);
    eig_stop = (int)rint(t1 + t2);
    eig_step = eig_stop - eig_start;
    if(THISPE == (FS_NPES - 1)) eig_stop = n;

    // Set width of window in terms of a percentage of n. Larger values will be slower but
    // exhibit behavior closer to full diagonalization.
    if((ct.folded_spectrum_width < 0.15) || (ct.folded_spectrum_width > 1.0)) {
        rmg_printf("Folded spectrum width of %8.4f is outside valid range (0.15,1.0). Resetting to default of 0.3.\n", ct.folded_spectrum_width);
        ct.folded_spectrum_width = 0.3;
    }

    double r_width = ct.folded_spectrum_width;
    t1 = (double)n;
    n_win = (int)(r_width * t1);

    // Find start of interval
    int ix = n_win - eig_step;
    if(ix < 4)
        rmg_error_handler(__FILE__, __LINE__, "Too few PE's to use folded spectrum method for this problem");
    if(ix % 2) {
        n_win = n_win + 1;
        ix = n_win - eig_step;
    }

    n_start = eig_start - ix/2;
    // force n_start to begin on the closest lower block boundary
    n_start = blocksize * (n_start / blocksize);
//std::cout << "BBBBBBBB " << n_start << "  " << blocksize << std::endl;
    if(n_start < 0) n_start = 0;
    if((n_start + n_win) > n) {
        if(blocksize == 1) {
            n_start = n - n_win;
        }
        else {
            n_start = blocksize * (n / blocksize);
            n_win = n - n_start;
        }
    }

}

