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

// Batched multigrid preconditioner for the Davidson residuals of k-point calculations.
//
// It applies the same V-cycle as DavPreconditionerOne (2nd order combined operator, Chebyshev
// smoothing, full weighting restriction, linear prolongation, maxpts sweeps on the coarsest level)
// to a batch of states at once, so that every halo exchange carries all the states of the batch.
// With SYCL the grids live on the device and every step is a kernel; otherwise the same loops
// run on the OpenMP threads, which is how the batched version is checked against the original.
// Only cubic, orthorhombic and tetragonal primitive cells, collinear spin and coalesce_factor 1.

#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <map>
#include "FiniteDiff.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "rmg_reduce.h"
#include "Kpoint.h"
#include "rmg_mgrid.h"
#include "Solvers.h"
#include "TradeImages.h"
#include "transition.h"

namespace {

struct cf { float re, im; };

struct MgLevel {
    int dx, dy, dz;          // points on this rank
    int xo, yo, zo;          // offsets of the restriction to the next level (MG_SIZE)
    int nx, ny, nz;          // global points
    double hx, hy, hz;
    size_t sz;               // (dx+2)*(dy+2)*(dz+2)
    float th2;               // diagonal of the 2nd order combined operator
    cf cpx, cmx, cpy, cmy, cpz, cmz;
    int presweeps;
    std::vector<double> pcoefs;
    bool bottom;
};

class BatchMg {
public:
    BatchMg(Kpoint<std::complex<double>> *kptr, int levels, double k, int maxb);
    ~BatchMg(void);
    bool ok = false;
    void set_k(double k_in) { k = k_in; }
    void run(std::complex<double> *res, double *eigs, int nst);

private:
    Kpoint<std::complex<double>> *kptr;
    std::vector<MgLevel> L;
    double k;
    int maxb;
    bool dev = false;
    MPI_Comm comm;
    int nb[6];               // W E S N D U neighbours
    std::vector<cf *> v, f, r;
    cf *sbuf = nullptr, *rbuf = nullptr;
    size_t hbuf_n = 0;
    cf *host_stage = nullptr;
#if SYCL_ENABLED
    sycl::queue *q = nullptr;
#endif

    template <typename F> void pfor(size_t n, F fn);
    void sync(void);
    cf *alloc(size_t n);
    void release(cf *p);
    void trade(cf *a, int l, int nst);
    void solve(int l, int nst);
};

template <typename F> void BatchMg::pfor(size_t n, F fn)
{
    if(n == 0) return;
#if SYCL_ENABLED
    if(dev)
    {
        q->parallel_for(sycl::range<1>(n), [=](sycl::id<1> i) { fn((size_t)i[0]); });
        return;
    }
#endif
#pragma omp parallel for schedule(static)
    for(size_t i = 0;i < n;i++) fn(i);
}

void BatchMg::sync(void)
{
#if SYCL_ENABLED
    if(dev) q->wait();
#endif
}

cf *BatchMg::alloc(size_t n)
{
#if SYCL_ENABLED
    if(dev) return sycl::malloc_device<cf>(std::max(n, (size_t)1), *q);
#endif
    return new cf[std::max(n, (size_t)1)];
}

void BatchMg::release(cf *p)
{
    if(!p) return;
#if SYCL_ENABLED
    if(dev) { sycl::free(p, *q); return; }
#endif
    delete [] p;
}

BatchMg::BatchMg(Kpoint<std::complex<double>> *kptr_in, int levels, double k_in, int maxb_in) :
    kptr(kptr_in), k(k_in), maxb(maxb_in)
{
    rmg::grid *G = kptr->G;
    Lattice *Lat = kptr->L;
    int ibrav = Lat->get_ibrav_type();
    if(ibrav != CUBIC_PRIMITIVE && ibrav != ORTHORHOMBIC_PRIMITIVE && ibrav != TETRAGONAL_PRIMITIVE) return;

    rmg::mgrid MG(Lat, kptr->T, G, 1, ct.max_zvalence);
    MG.set_kpoints(kptr->kp.kvec, kptr->kp.kmag);
    FiniteDiff FD(Lat);

    int pxdim = G->get_PX0_GRID(1), pydim = G->get_PY0_GRID(1), pzdim = G->get_PZ0_GRID(1);
    int gx = G->get_NX_GRID(1), gy = G->get_NY_GRID(1), gz = G->get_NZ_GRID(1);
    int ox = G->get_PX_OFFSET(1), oy = G->get_PY_OFFSET(1), oz = G->get_PZ_OFFSET(1);
    double hx0 = G->get_hxgrid(1), hy0 = G->get_hygrid(1), hz0 = G->get_hzgrid(1);

    int dx = pxdim, dy = pydim, dz = pzdim;
    for(int level = 0;;level++)
    {
        MgLevel lv;
        lv.dx = dx; lv.dy = dy; lv.dz = dz;
        int fac = 1 << level;
        lv.nx = gx / fac; lv.ny = gy / fac; lv.nz = gz / fac;
        lv.hx = hx0 * fac; lv.hy = hy0 * fac; lv.hz = hz0 * fac;
        lv.sz = (size_t)(dx + 2) * (size_t)(dy + 2) * (size_t)(dz + 2);
        int dx2 = MG.MG_SIZE(dx, level, gx, ox, pxdim, &lv.xo, PERIODIC);
        int dy2 = MG.MG_SIZE(dy, level, gy, oy, pydim, &lv.yo, PERIODIC);
        int dz2 = MG.MG_SIZE(dz, level, gz, oz, pzdim, &lv.zo, PERIODIC);
        lv.bottom = (level >= levels) || (dx2 < 0) || (dy2 < 0) || (dz2 < 0);

        std::complex<float> cmx[12], cpx[12], cmy[12], cpy[12], cmz[12], cpz[12];
        FD.fd_combined_coeffs(2, lv.hx, 0, cmx, cpx, kptr->kp.kvec);
        FD.fd_combined_coeffs(2, lv.hx, 1, cmy, cpy, kptr->kp.kvec);
        FD.fd_combined_coeffs(2, lv.hx, 2, cmz, cpz, kptr->kp.kvec);
        lv.th2 = (float)FD.fd_coeff0(2, lv.hx);
        lv.cpx = {cpx[0].real(), cpx[0].imag()}; lv.cmx = {cmx[0].real(), cmx[0].imag()};
        lv.cpy = {cpy[0].real(), cpy[0].imag()}; lv.cmy = {cmy[0].real(), cmy[0].imag()};
        lv.cpz = {cpz[0].real(), cpz[0].imag()}; lv.cmz = {cmz[0].real(), cmz[0].imag()};

        if(lv.bottom)
        {
            lv.presweeps = std::max(std::max(lv.nx, lv.ny), lv.nz);
            lv.pcoefs = MG.pois_chebyshev_coeffs(lv.nx, lv.ny, lv.nz, lv.hx, lv.hy, lv.hz, 0.0, lv.presweeps);
        }
        else
        {
            lv.presweeps = MG.pre_cyc[level];
            lv.pcoefs = MG.pois_chebyshev_coeffs(lv.nx, lv.ny, lv.nz, lv.hx, lv.hy, lv.hz, 1.0, lv.presweeps);
            if(MG.post_cyc[level] != lv.presweeps || MG.mu_cyc[level] != 1) return;
        }
        L.push_back(lv);
        if(lv.bottom) break;
        dx = dx2; dy = dy2; dz = dz2;
    }

    comm = pct.grid_comm;
    int ii, jj, kk;
    G->pe2xyz(pct.gridpe, &ii, &jj, &kk);
    int PX = G->get_PE_X(), PY = G->get_PE_Y(), PZ = G->get_PE_Z();
    nb[0] = G->xyz2pe((ii - 1 + PX) % PX, jj, kk);
    nb[1] = G->xyz2pe((ii + 1) % PX, jj, kk);
    nb[2] = G->xyz2pe(ii, (jj - 1 + PY) % PY, kk);
    nb[3] = G->xyz2pe(ii, (jj + 1) % PY, kk);
    nb[4] = G->xyz2pe(ii, jj, (kk - 1 + PZ) % PZ);
    nb[5] = G->xyz2pe(ii, jj, (kk + 1) % PZ);

#if SYCL_ENABLED
    dev = true;
    q = new sycl::queue(ct.sycl_Q.get_context(), ct.sycl_Q.get_device(), sycl::property::queue::in_order());
#endif

    for(size_t l = 0;l < L.size();l++)
    {
        v.push_back(alloc((size_t)maxb * L[l].sz));
        f.push_back(alloc((size_t)maxb * L[l].sz));
        r.push_back(alloc((size_t)maxb * L[l].sz));
    }
    // Largest face, including the image points of the other two directions
    size_t face = 0;
    face = std::max(face, (size_t)(L[0].dy + 2) * (L[0].dz + 2));
    face = std::max(face, (size_t)(L[0].dx + 2) * (L[0].dz + 2));
    face = std::max(face, (size_t)(L[0].dx + 2) * (L[0].dy + 2));
    hbuf_n = 2 * face * maxb;
    sbuf = alloc(hbuf_n);
    rbuf = alloc(hbuf_n);
    host_stage = new cf[(size_t)maxb * L[0].sz];
    ok = true;
}

BatchMg::~BatchMg(void)
{
    for(auto p : v) release(p);
    for(auto p : f) release(p);
    for(auto p : r) release(p);
    release(sbuf);
    release(rbuf);
    delete [] host_stage;
#if SYCL_ENABLED
    if(q) delete q;
#endif
}

// Periodic halo exchange of all nst grids of level l: x, then y including the x images,
// then z including the x and y images, like a full trade_images.
void BatchMg::trade(cf *a, int l, int nst)
{
    const int dx = L[l].dx, dy = L[l].dy, dz = L[l].dz;
    const size_t sz = L[l].sz;
    const int incx = (dy + 2) * (dz + 2), incy = dz + 2;
    for(int dir = 0;dir < 3;dir++)
    {
        int n1, n2;        // face extents (slow, fast) and the first/last interior plane
        if(dir == 0) { n1 = dy + 2; n2 = dz + 2; }
        else if(dir == 1) { n1 = dx + 2; n2 = dz + 2; }
        else { n1 = dx + 2; n2 = dy + 2; }
        const int dlen = (dir == 0) ? dx : (dir == 1) ? dy : dz;
        const size_t fsz = (size_t)n1 * (size_t)n2;
        const size_t nf = fsz * (size_t)nst;
        cf *s = sbuf, *rr = rbuf;
        // face index -> grid index for plane p of direction dir
        auto gidx = [=](size_t st, size_t i1, size_t i2, int p) -> size_t {
            if(dir == 0) return st*sz + (size_t)p*incx + i1*incy + i2;
            if(dir == 1) return st*sz + i1*incx + (size_t)p*incy + i2;
            return st*sz + i1*incx + i2*incy + (size_t)p;
        };
        // Pack the low interior plane (1) into s[0..nf) and the high interior plane (dlen) into s[nf..2nf)
        pfor(nf, [=](size_t i) {
            size_t st = i / fsz, rem = i % fsz, i1 = rem / n2, i2 = rem % n2;
            s[i] = a[gidx(st, i1, i2, 1)];
            s[nf + i] = a[gidx(st, i1, i2, dlen)];
        });
        sync();
        int lo = nb[2*dir], hi = nb[2*dir+1];
        if(lo == pct.gridpe && hi == pct.gridpe)
        {
            pfor(nf, [=](size_t i) { rr[i] = s[nf + i]; rr[nf + i] = s[i]; });
        }
        else
        {
            MPI_Request req[4];
            int cnt = (int)(2 * nf);   // floats
            // rr[0..nf) <- high plane of the low neighbour (my image 0); rr[nf..) <- low plane of the high neighbour
            MPI_Irecv((float *)rr, cnt, MPI_FLOAT, lo, 101 + 2*dir, comm, &req[0]);
            MPI_Irecv((float *)(rr + nf), cnt, MPI_FLOAT, hi, 100 + 2*dir, comm, &req[1]);
            MPI_Isend((float *)s, cnt, MPI_FLOAT, lo, 100 + 2*dir, comm, &req[2]);
            MPI_Isend((float *)(s + nf), cnt, MPI_FLOAT, hi, 101 + 2*dir, comm, &req[3]);
            MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
        }
        pfor(nf, [=](size_t i) {
            size_t st = i / fsz, rem = i % fsz, i1 = rem / n2, i2 = rem % n2;
            a[gidx(st, i1, i2, 0)] = rr[i];
            a[gidx(st, i1, i2, dlen + 1)] = rr[nf + i];
        });
        sync();
    }
}

void BatchMg::solve(int l, int nst)
{
    MgLevel &lv = L[l];
    const int dx = lv.dx, dy = lv.dy, dz = lv.dz;
    const size_t sz = lv.sz;
    const int incx = (dy + 2) * (dz + 2), incy = dz + 2;
    const size_t nall = sz * nst;
    const size_t nin = (size_t)dx * dy * dz * nst;
    cf *V = v[l], *F = f[l], *R = r[l];
    const float th2 = lv.th2;
    const cf cpx = lv.cpx, cmx = lv.cmx, cpy = lv.cpy, cmy = lv.cmy, cpz = lv.cpz, cmz = lv.cmz;
    const float kk = (float)k;
    const double pscale = std::pow(0.5, l);
    const double Zfac = ct.max_zvalence * lv.hx * lv.hx;
    const double diag = -(double)lv.th2;

    // interior index -> grid index
    auto in2g = [=](size_t i) -> size_t {
        size_t st = i / ((size_t)dx*dy*dz), rem = i % ((size_t)dx*dy*dz);
        size_t ix = rem / ((size_t)dy*dz) + 1, iy = (rem / dz) % dy + 1, iz = rem % dz + 1;
        return st*sz + ix*incx + iy*incy + iz;
    };
    // A v at grid index g (2nd order combined operator), complex float
    auto apply = [=](const cf *a, size_t g) -> cf {
        cf c = a[g], xp = a[g+incx], xm = a[g-incx], yp = a[g+incy], ym = a[g-incy], zp = a[g+1], zm = a[g-1];
        float re = th2*c.re
                 + cpz.re*zp.re - cpz.im*zp.im + cmz.re*zm.re - cmz.im*zm.im
                 + cpy.re*yp.re - cpy.im*yp.im + cmy.re*ym.re - cmy.im*ym.im
                 + cpx.re*xp.re - cpx.im*xp.im + cmx.re*xm.re - cmx.im*xm.im;
        float im = th2*c.im
                 + cpz.re*zp.im + cpz.im*zp.re + cmz.re*zm.im + cmz.im*zm.re
                 + cpy.re*yp.im + cpy.im*yp.re + cmy.re*ym.im + cmy.im*ym.re
                 + cpx.re*xp.im + cpx.im*xp.re + cmx.re*xm.im + cmx.im*xm.re;
        return cf{re, im};
    };
    // One Jacobi/Chebyshev sweep: v += scale*(A v + k v - f). A v needs the old v of the
    // neighbours, so it goes to R first.
    auto sweep = [&](double step) {
        float scale = (float)(step / (diag + Zfac));
        pfor(nin, [=](size_t i) { size_t g = in2g(i); R[g] = apply(V, g); });
        pfor(nin, [=](size_t i) {
            size_t g = in2g(i);
            cf av = R[g], vv = V[g], ff = F[g];
            V[g].re = vv.re + scale*(av.re + kk*vv.re - ff.re);
            V[g].im = vv.im + scale*(av.im + kk*vv.im - ff.im);
        });
    };

    trade(F, l, nst);
    pfor(nall, [=](size_t i) { V[i] = cf{0.0f, 0.0f}; });
    for(int c = 0;c < lv.presweeps;c++)
    {
        sweep(pscale * lv.pcoefs[c]);
        trade(V, l, nst);
    }
    if(lv.bottom) return;

    // residual R = f - (A v + k v)
    pfor(nin, [=](size_t i) {
        size_t g = in2g(i);
        cf av = apply(V, g), vv = V[g], ff = F[g];
        R[g] = cf{ff.re - (av.re + kk*vv.re), ff.im - (av.im + kk*vv.im)};
    });
    trade(R, l, nst);

    // full weighting restriction into the interior of f[l+1]
    {
        MgLevel &nl = L[l+1];
        const int dx2 = nl.dx, dy2 = nl.dy, dz2 = nl.dz;
        const size_t sz2 = nl.sz;
        const int incx2 = (dy2 + 2) * (dz2 + 2), incy2 = dz2 + 2;
        const int xo = lv.xo, yo = lv.yo, zo = lv.zo;
        cf *F2 = f[l+1];
        const size_t n2 = (size_t)dx2 * dy2 * dz2 * nst;
        pfor(n2, [=](size_t i) {
            size_t st = i / ((size_t)dx2*dy2*dz2), rem = i % ((size_t)dx2*dy2*dz2);
            int ix = (int)(rem / ((size_t)dy2*dz2)) + 1, iy = (int)((rem / dz2) % dy2) + 1, iz = (int)(rem % dz2) + 1;
            int x0 = 2*ix - 1 + xo, y0 = 2*iy - 1 + yo, z0 = 2*iz - 1 + zo;
            const cf *A = R + st*sz;
            float sr[4] = {0,0,0,0}, si[4] = {0,0,0,0};   // center, face, edge, corner
            for(int ddx = -1;ddx <= 1;ddx++)
            for(int ddy = -1;ddy <= 1;ddy++)
            for(int ddz = -1;ddz <= 1;ddz++)
            {
                int nz = (ddx != 0) + (ddy != 0) + (ddz != 0);
                cf val = A[(size_t)(x0+ddx)*incx + (size_t)(y0+ddy)*incy + (size_t)(z0+ddz)];
                sr[nz] += val.re; si[nz] += val.im;
            }
            const float scale = 1.0f / 64.0f;
            F2[st*sz2 + (size_t)ix*incx2 + (size_t)iy*incy2 + iz] =
                cf{scale*(8.0f*sr[0] + 4.0f*sr[1] + 2.0f*sr[2] + sr[3]),
                   scale*(8.0f*si[0] + 4.0f*si[1] + 2.0f*si[2] + si[3])};
        });
    }

    solve(l + 1, nst);

    // linear prolongation of v[l+1] into R, then v += R
    {
        MgLevel &nl = L[l+1];
        const int dy2 = nl.dy, dz2 = nl.dz;
        const size_t sz2 = nl.sz;
        const int incxr = (dy2 + 2) * (dz2 + 2), incyr = dz2 + 2;
        const int xo = lv.xo, yo = lv.yo, zo = lv.zo;
        const cf *H = v[l+1];
        // coarse points: ix in [1-xo, dx/2+1] etc.
        const int ax0 = 1 - xo, ay0 = 1 - yo, az0 = 1 - zo;
        const int ax1 = dx/2 + 1, ay1 = dy/2 + 1, az1 = dz/2 + 1;
        const size_t cnx = std::max(0, ax1 - ax0 + 1), cny = std::max(0, ay1 - ay0 + 1), cnz = std::max(0, az1 - az0 + 1);
        const size_t nc = cnx * cny * cnz * nst;
        pfor(nc, [=](size_t i) {
            size_t st = i / (cnx*cny*cnz), rem = i % (cnx*cny*cnz);
            int ix = ax0 + (int)(rem / (cny*cnz)), iy = ay0 + (int)((rem / cnz) % cny), iz = az0 + (int)(rem % cnz);
            R[st*sz + (size_t)(2*ix - 1 + xo)*incx + (size_t)(2*iy - 1 + yo)*incy + (size_t)(2*iz - 1 + zo)] =
                H[st*sz2 + (size_t)ix*incxr + (size_t)iy*incyr + iz];
        });
        // The other fine points only read coarse ones. Parity of (ix-(1-xo)), (iy-(1-yo)), (iz-(1-zo))
        // selects the stencil: odd in a direction means average the two coarse neighbours in it.
        const int bx = 1 - xo, by = 1 - yo, bz = 1 - zo;
        const int fnx = dx + 1 - bx, fny = dy + 1 - by, fnz = dz + 1 - bz;   // ix in [bx, dx]
        const size_t nfp = (size_t)std::max(fnx,0) * std::max(fny,0) * std::max(fnz,0) * nst;
        pfor(nfp, [=](size_t i) {
            size_t st = i / ((size_t)fnx*fny*fnz), rem = i % ((size_t)fnx*fny*fnz);
            int ix = bx + (int)(rem / ((size_t)fny*fnz)), iy = by + (int)((rem / fnz) % fny), iz = bz + (int)(rem % fnz);
            int px = (ix - bx) & 1, py = (iy - by) & 1, pz = (iz - bz) & 1;
            if(!px && !py && !pz) return;     // coarse point
            cf *A = R + st*sz;
            float w = 1.0f / (float)(1 << (px + py + pz));
            float sr = 0.0f, si = 0.0f;
            for(int ddx = -px;ddx <= px;ddx += 2)
            for(int ddy = -py;ddy <= py;ddy += 2)
            for(int ddz = -pz;ddz <= pz;ddz += 2)
            {
                cf val = A[(size_t)(ix+ddx)*incx + (size_t)(iy+ddy)*incy + (size_t)(iz+ddz)];
                sr += val.re; si += val.im;
            }
            A[(size_t)ix*incx + (size_t)iy*incy + iz] = cf{w*sr, w*si};
        });
        pfor(nin, [=](size_t i) {
            size_t g = in2g(i);
            V[g].re += R[g].re; V[g].im += R[g].im;
        });
    }
    trade(V, l, nst);
    for(int c = 0;c < lv.presweeps;c++)
    {
        sweep(pscale * lv.pcoefs[c]);
        if(c < lv.presweeps - 1) trade(V, l, nst);
    }
    trade(V, l, nst);
}

void BatchMg::run(std::complex<double> *res, double *eigs, int nstates)
{
    rmg::grid *G = kptr->G;
    const int dx = L[0].dx, dy = L[0].dy, dz = L[0].dz;
    const size_t pbasis = (size_t)dx * dy * dz;
    const size_t sz = L[0].sz;
    const int incx = (dy + 2) * (dz + 2), incy = dz + 2;
    const double gbasis = (double)G->get_NX_GRID(1) * (double)G->get_NY_GRID(1) * (double)G->get_NZ_GRID(1);

    // Cell neutralization: mean of the real part of each residual
    std::vector<double> t1(nstates, 0.0);
#pragma omp parallel for
    for(int st = 0;st < nstates;st++)
    {
        double s = 0.0;
        for(size_t idx = 0;idx < pbasis;idx++) s += std::real(res[st*pbasis + idx]);
        t1[st] = s;
    }
    MPI_Allreduce(MPI_IN_PLACE, t1.data(), nstates, MPI_DOUBLE, MPI_SUM, comm);
    for(int st = 0;st < nstates;st++) t1[st] /= gbasis;

    for(int st0 = 0;st0 < nstates;st0 += maxb)
    {
        int nst = std::min(maxb, nstates - st0);
        // f[0] = res - t1 in single precision, image points zero
#pragma omp parallel for
        for(int s = 0;s < nst;s++)
        {
            cf *h = host_stage + (size_t)s*sz;
            std::fill(h, h + sz, cf{0.0f, 0.0f});
            const std::complex<double> *src = res + (size_t)(st0 + s)*pbasis;
            double t = t1[st0 + s];
            for(int ix = 0;ix < dx;ix++)
                for(int iy = 0;iy < dy;iy++)
                    for(int iz = 0;iz < dz;iz++)
                    {
                        std::complex<double> c = src[((size_t)ix*dy + iy)*dz + iz] - t;
                        h[(size_t)(ix+1)*incx + (size_t)(iy+1)*incy + iz + 1] = cf{(float)c.real(), (float)c.imag()};
                    }
        }
#if SYCL_ENABLED
        if(dev) q->memcpy(f[0], host_stage, (size_t)nst*sz*sizeof(cf)).wait();
        else
#endif
        std::copy(host_stage, host_stage + (size_t)nst*sz, f[0]);

        solve(0, nst);

#if SYCL_ENABLED
        if(dev) q->memcpy(host_stage, v[0], (size_t)nst*sz*sizeof(cf)).wait();
        else
#endif
        std::copy(v[0], v[0] + (size_t)nst*sz, host_stage);

#pragma omp parallel for
        for(int s = 0;s < nst;s++)
        {
            const cf *h = host_stage + (size_t)s*sz;
            std::complex<double> *dst = res + (size_t)(st0 + s)*pbasis;
            double shift = eigs[st0 + s] * t1[st0 + s];
            for(int ix = 0;ix < dx;ix++)
                for(int iy = 0;iy < dy;iy++)
                    for(int iz = 0;iz < dz;iz++)
                    {
                        cf c = h[(size_t)(ix+1)*incx + (size_t)(iy+1)*incy + iz + 1];
                        dst[((size_t)ix*dy + iy)*dz + iz] = std::complex<double>((double)c.re + shift, (double)c.im);
                    }
        }
    }
}

}  // namespace

// Returns false when the batched preconditioner does not apply; the caller then uses the original one.
// Input key batch_preconditioner.
bool DavPreconditionerBatch(Kpoint<std::complex<double>> *kptr, std::complex<double> *res, double *eigs,
                            int notconv, double avg_potential)
{
    if(!ct.batch_preconditioner) return false;
    if(ct.noncoll_factor != 1 || pct.coalesce_factor != 1 || notconv < 1) return false;

    // One instance per k-point object (the operator coefficients depend on the k-vector)
    static std::map<void *, BatchMg *> instances;
    BatchMg *&B = instances[(void *)kptr];
    if(!B) B = new BatchMg(kptr, ct.eig_parm.levels, -avg_potential, 64);
    if(!B->ok) return false;
    B->set_k(-avg_potential);

    B->run(res, eigs, notconv);
    return true;
}
