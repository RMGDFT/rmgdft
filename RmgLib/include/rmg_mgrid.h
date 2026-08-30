/*
 *
 * Copyright (c) 2014, Emil Briggs
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
*/


#ifndef RMG_mgrid_H
#define RMG_mgrid_H 1

/* Maximum number of multigrid levels */
#define         MAX_MG_LEVELS   8

#ifdef __cplusplus
#include "Lattice.h"
#include "TradeImages.h"
#include "rmg_error.h"
#include "boundary_conditions.h"

namespace rmg
{

    class mgrid {

    private:
        Lattice *L;
        TradeImages *T;
        rmg::grid *G;
        int density;
        double zmax;
        int ibrav;
        int level_flag;
        bool central_trade;
        static int level_warning;
        std::vector<double> kvec = {0.0, 0.0, 0.0};
        double kmag=0.0;
        int boundary_flag = PERIODIC;   // only thing supported for now

        std::array<double, MAX_MG_LEVELS> hx;
        std::array<double, MAX_MG_LEVELS> hy;
        std::array<double, MAX_MG_LEVELS> hz;


        // Timer mode 0=off (default) 1=on
        bool timer_mode;

    public:
        mgrid(Lattice *lptr, TradeImages *tptr, rmg::grid *G, int density_in, double zmax_in);
       ~mgrid(void);
       

        // Public in case the caller wants to adjust these for some reasons
        std::array<int, MAX_MG_LEVELS> pre_cyc = {3, 3, 3, 3, 3, 3, 3, 3};
        std::array<int, MAX_MG_LEVELS> post_cyc = {3, 3, 3, 3, 3, 3, 3, 3};
        std::array<int, MAX_MG_LEVELS> mu_cyc = {1, 1, 1, 1, 1, 1, 1, 1};
        std::array<std::vector<double>, MAX_MG_LEVELS> rms_residuals;


        // Level 0 grid offsets and dimensions
        int gxsize;
        int gysize;
        int gzsize;
        int gxoffset;
        int gyoffset;
        int gzoffset;

        // This vector holds the maximum offset usable for any given multigrid level.
        static std::vector<int> toffsets;

        void set_timer_mode(bool verbose);
        void set_kpoints(double *kvec_in, double kmag_in)
        {
            kvec[0] = kvec_in[0];
            kvec[1] = kvec_in[1];
            kvec[2] = kvec_in[2];
            kmag = kmag_in;
        }

        double rel_sradius(int level);

        template <typename RmgType> void anchor_residual(int id, int level, int n, RmgType *r);
        template <typename RmgType> void anchor_residual(int level, int dimx, int dimy, int dimz, RmgType *r);

        template <typename RmgType> void mg_restrict (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

        template <typename RmgType> void mg_prolong (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

        template <typename RmgType> void mg_prolong_cubic (RmgType * full, RmgType * half, int dimx, int dimy, int dimz, int dx2, int dy2, int dz2, int xoffset, int yoffset, int zoffset);

        template <typename RmgType> void eval_residual (RmgType * mat, RmgType * f_mat, RmgType *work, int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, RmgType * res, double *pot);

        template <typename RmgType> void solv_pois (RmgType * vmat, RmgType * fmat, RmgType * work,
                    int dimx, int dimy, int dimz, double gridhx, double gridhy, double gridhz, double step, double k, double *pot);


        int MG_SIZE (int curdim, int curlevel, int global_dim, int global_offset, int global_pdim, int *roffset, int bctype);

        template <typename RmgType> void mgrid_solv (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                     int dimx, int dimy, int dimz,
                     int level, int max_levels, double k, double *pot, int pxdim, int pydim, int pzdim);

        template <typename RmgType> void mgrid_solv_pois (RmgType * v_mat, RmgType * f_mat, RmgType * work,
                     int dimx, int dimy, int dimz,
                     int level, int max_levels, int pxdim, int pydim, int pzdim);

    };

}

#endif
#endif

