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

// Batched charge density for k-point calculations.
//
// Same result as the GetNewRhoOne loop of GetNewRhoPre: every occupied orbital is interpolated
// to the fine grid with the ratio 2 Prolong stencil in single precision, normalised on the fine
// grid and |psi|^2 is accumulated with its weight. Here a batch of states shares each halo
// exchange and every step is one kernel. States with zero weight contribute exactly zero and are
// skipped. With SYCL the grids live on the device; otherwise the same loops run on the OpenMP
// threads, which is how the batched version is checked against the original.
// Only cubic, orthorhombic and tetragonal primitive cells, collinear spin, coalesce_factor 1,
// no meta-GGA.

#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <array>
#include <omp.h>
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "transition.h"
#include "RmgTimer.h"
#include "Kpoint.h"
#include "Prolong.h"
#include "Solvers.h"
#include "common_prototypes.h"

namespace {

struct cf { float re, im; };

class BatchRho {
public:
    BatchRho(Prolong &P, int maxb);
    ~BatchRho(void);
    bool ok = false;
    template <typename KpointType> void run(Kpoint<KpointType> **Kpts, double *work);

private:
    int maxb, ord, w;
    int hx, hy, hz;          // coarse points on this rank
    size_t csz;              // padded coarse grid (hx+ord)*(hy+ord)*(hz+ord)
    size_t asz, bsz, fsz;    // after the x pass, after the y pass, fine grid
    float a0, a1[MAX_PROLONG_ORDER];
    bool dev = false;
    MPI_Comm comm;
    int nb[6];               // W E S N D U neighbours
    cf *c = nullptr, *A = nullptr, *B = nullptr;
    double *dens = nullptr, *acc = nullptr, *sums = nullptr;
    cf *sbuf = nullptr, *rbuf = nullptr;
    cf *host_stage = nullptr;
#if SYCL_ENABLED
    sycl::queue *q = nullptr;
#endif

    template <typename F> void pfor(size_t n, F fn);
    void sync(void);
    template <typename T> T *alloc(size_t n);
    template <typename T> void release(T *p);
    void trade(int nst);
    void batch(int nst, std::array<double, 64> &wgt, bool normalise);
};

template <typename F> void BatchRho::pfor(size_t n, F fn)
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

void BatchRho::sync(void)
{
#if SYCL_ENABLED
    if(dev) q->wait();
#endif
}

template <typename T> T *BatchRho::alloc(size_t n)
{
#if SYCL_ENABLED
    if(dev) return sycl::malloc_device<T>(std::max(n, (size_t)1), *q);
#endif
    return new T[std::max(n, (size_t)1)];
}

template <typename T> void BatchRho::release(T *p)
{
    if(!p) return;
#if SYCL_ENABLED
    if(dev) { sycl::free(p, *q); return; }
#endif
    delete [] p;
}

BatchRho::BatchRho(Prolong &P, int maxb_in) : maxb(maxb_in)
{
    int ibrav = Rmg_L.get_ibrav_type();
    if(ibrav != CUBIC_PRIMITIVE && ibrav != ORTHORHOMBIC_PRIMITIVE && ibrav != TETRAGONAL_PRIMITIVE) return;
    if(P.ratio != 2 || P.order < 2 || P.order > MAX_PROLONG_ORDER || (P.order % 2)) return;
    ord = P.order;
    w = ord / 2;
    a0 = P.af[0][w - 1];
    for(int k = 0;k < ord;k++) a1[k] = P.af[1][k];

    hx = Rmg_G->get_PX0_GRID(1); hy = Rmg_G->get_PY0_GRID(1); hz = Rmg_G->get_PZ0_GRID(1);
    // The halo of width ord/2 comes from the nearest neighbour only
    if(hx < w || hy < w || hz < w) return;
    csz = (size_t)(hx + ord) * (hy + ord) * (hz + ord);
    asz = (size_t)(2*hx) * (hy + ord) * (hz + ord);
    bsz = (size_t)(2*hx) * (2*hy) * (hz + ord);
    fsz = (size_t)(2*hx) * (2*hy) * (2*hz);

    comm = pct.grid_comm;
    int ii, jj, kk;
    Rmg_G->pe2xyz(pct.gridpe, &ii, &jj, &kk);
    int PX = Rmg_G->get_PE_X(), PY = Rmg_G->get_PE_Y(), PZ = Rmg_G->get_PE_Z();
    nb[0] = Rmg_G->xyz2pe((ii - 1 + PX) % PX, jj, kk);
    nb[1] = Rmg_G->xyz2pe((ii + 1) % PX, jj, kk);
    nb[2] = Rmg_G->xyz2pe(ii, (jj - 1 + PY) % PY, kk);
    nb[3] = Rmg_G->xyz2pe(ii, (jj + 1) % PY, kk);
    nb[4] = Rmg_G->xyz2pe(ii, jj, (kk - 1 + PZ) % PZ);
    nb[5] = Rmg_G->xyz2pe(ii, jj, (kk + 1) % PZ);

#if SYCL_ENABLED
    dev = true;
    q = new sycl::queue(ct.sycl_Q.get_context(), ct.sycl_Q.get_device(), sycl::property::queue::in_order());
#endif

    c = alloc<cf>((size_t)maxb * csz);
    A = alloc<cf>((size_t)maxb * asz);
    B = alloc<cf>((size_t)maxb * bsz);
    dens = alloc<double>((size_t)maxb * fsz);
    acc = alloc<double>(fsz);
    sums = alloc<double>(maxb);
    size_t face = 0;
    face = std::max(face, (size_t)(hy + ord) * (hz + ord));
    face = std::max(face, (size_t)(hx + ord) * (hz + ord));
    face = std::max(face, (size_t)(hx + ord) * (hy + ord));
    sbuf = alloc<cf>(2 * face * w * maxb);
    rbuf = alloc<cf>(2 * face * w * maxb);
    host_stage = new cf[(size_t)maxb * csz];
    ok = true;
}

BatchRho::~BatchRho(void)
{
    release(c); release(A); release(B); release(dens); release(acc); release(sums);
    release(sbuf); release(rbuf);
    delete [] host_stage;
#if SYCL_ENABLED
    if(q) delete q;
#endif
}

// Periodic halo exchange of width w for the nst coarse grids in c: x, then y including the x
// images, then z including the x and y images, like trade_imagesx with FULL_TRADE.
void BatchRho::trade(int nst)
{
    const int px = hx + ord, py = hy + ord, pz = hz + ord;
    const size_t sz = csz;
    const int incx = py * pz, incy = pz;
    const int ww = w;
    cf *a = c;
    for(int dir = 0;dir < 3;dir++)
    {
        int n1, n2;
        if(dir == 0) { n1 = py; n2 = pz; }
        else if(dir == 1) { n1 = px; n2 = pz; }
        else { n1 = px; n2 = py; }
        const int dlen = (dir == 0) ? hx : (dir == 1) ? hy : hz;
        const size_t fs = (size_t)n1 * n2;
        const size_t pl = fs * ww;               // w planes of one state
        const size_t nf = pl * (size_t)nst;
        cf *s = sbuf, *rr = rbuf;
        auto gidx = [=](size_t st, size_t i1, size_t i2, int p) -> size_t {
            if(dir == 0) return st*sz + (size_t)p*incx + i1*incy + i2;
            if(dir == 1) return st*sz + i1*incx + (size_t)p*incy + i2;
            return st*sz + i1*incx + i2*incy + (size_t)p;
        };
        // low interior planes w..2w-1 into s[0..nf), high interior planes dlen..dlen+w-1 into s[nf..2nf)
        pfor(nf, [=](size_t i) {
            size_t st = i / pl, rem = i % pl, p = rem / fs, r2 = rem % fs, i1 = r2 / n2, i2 = r2 % n2;
            s[i] = a[gidx(st, i1, i2, ww + (int)p)];
            s[nf + i] = a[gidx(st, i1, i2, dlen + (int)p)];
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
            // rr[0..nf) <- high planes of the low neighbour; rr[nf..) <- low planes of the high neighbour
            MPI_Irecv((float *)rr, cnt, MPI_FLOAT, lo, 111 + 2*dir, comm, &req[0]);
            MPI_Irecv((float *)(rr + nf), cnt, MPI_FLOAT, hi, 110 + 2*dir, comm, &req[1]);
            MPI_Isend((float *)s, cnt, MPI_FLOAT, lo, 110 + 2*dir, comm, &req[2]);
            MPI_Isend((float *)(s + nf), cnt, MPI_FLOAT, hi, 111 + 2*dir, comm, &req[3]);
            MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
        }
        pfor(nf, [=](size_t i) {
            size_t st = i / pl, rem = i % pl, p = rem / fs, r2 = rem % fs, i1 = r2 / n2, i2 = r2 % n2;
            a[gidx(st, i1, i2, (int)p)] = rr[i];
            a[gidx(st, i1, i2, dlen + ww + (int)p)] = rr[nf + i];
        });
        sync();
    }
}

// Interpolates the nst states in c (interior filled) and adds wgt[st] * |psi_f|^2, divided by the
// fine grid norm when normalise is set, to acc.
void BatchRho::batch(int nst, std::array<double, 64> &wgt, bool normalise)
{
    trade(nst);

    const int o = ord, ic = w - 1;
    const float c0 = a0;
    std::array<float, MAX_PROLONG_ORDER> c1;
    for(int k = 0;k < MAX_PROLONG_ORDER;k++) c1[k] = (k < ord) ? a1[k] : 0.0f;
    const int dx = 2*hx, dy = 2*hy, dz = 2*hz;
    const int py = hy + ord, pz = hz + ord;
    const size_t cs = csz, as = asz, bs = bsz, fs = fsz;
    const cf *C = c;
    cf *AA = A, *BB = B;
    double *D = dens;

    // x pass: A[ix][iy][iz], ix fine, iy and iz padded coarse
    pfor((size_t)nst * as, [=](size_t i) {
        size_t st = i / as, rem = i % as;
        int ix = (int)(rem / ((size_t)py*pz)), iy = (int)((rem / pz) % py), iz = (int)(rem % pz);
        const cf *h = C + st*cs + (size_t)(ix/2 + 1)*py*pz + (size_t)iy*pz + iz;
        const size_t str = (size_t)py*pz;
        cf sum{0.0f, 0.0f};
        if(ix & 1)
            for(int k = 0;k < o;k++) { sum.re = sum.re + c1[k]*h[k*str].re; sum.im = sum.im + c1[k]*h[k*str].im; }
        else
            sum = cf{c0*h[ic*str].re, c0*h[ic*str].im};
        AA[i] = sum;
    });
    // y pass: B[ix][iy][iz], ix and iy fine, iz padded coarse
    pfor((size_t)nst * bs, [=](size_t i) {
        size_t st = i / bs, rem = i % bs;
        int ix = (int)(rem / ((size_t)dy*pz)), iy = (int)((rem / pz) % dy), iz = (int)(rem % pz);
        const cf *h = AA + st*as + (size_t)ix*py*pz + (size_t)(iy/2 + 1)*pz + iz;
        cf sum{0.0f, 0.0f};
        if(iy & 1)
            for(int k = 0;k < o;k++) { sum.re = sum.re + c1[k]*h[k*pz].re; sum.im = sum.im + c1[k]*h[k*pz].im; }
        else
            sum = cf{c0*h[ic*pz].re, c0*h[ic*pz].im};
        BB[i] = sum;
    });
    // z pass and |psi_f|^2 in double
    pfor((size_t)nst * fs, [=](size_t i) {
        size_t st = i / fs, rem = i % fs;
        int ix = (int)(rem / ((size_t)dy*dz)), iy = (int)((rem / dz) % dy), iz = (int)(rem % dz);
        const cf *h = BB + st*bs + (size_t)ix*dy*pz + (size_t)iy*pz + (iz/2 + 1);
        cf sum{0.0f, 0.0f};
        if(iz & 1)
            for(int k = 0;k < o;k++) { sum.re = sum.re + c1[k]*h[k].re; sum.im = sum.im + c1[k]*h[k].im; }
        else
            sum = cf{c0*h[ic].re, c0*h[ic].im};
        double re = (double)sum.re, im = (double)sum.im;
        D[i] = re*re + im*im;
    });

    std::vector<double> norm(nst, 1.0);
    if(normalise)
    {
#if SYCL_ENABLED
        if(dev)
        {
            q->memset(sums, 0, nst*sizeof(double));
            for(int st = 0;st < nst;st++)
            {
                const double *Ds = D + (size_t)st*fs;
                q->parallel_for(sycl::range<1>(fs), sycl::reduction(sums + st, sycl::plus<double>()),
                                [=](sycl::id<1> j, auto &s) { s += Ds[j[0]]; });
            }
            q->memcpy(norm.data(), sums, nst*sizeof(double)).wait();
        }
        else
#endif
        for(int st = 0;st < nst;st++)
        {
            double s = 0.0;
            const double *Ds = D + (size_t)st*fs;
#pragma omp parallel for reduction(+:s)
            for(size_t j = 0;j < fs;j++) s += Ds[j];
            norm[st] = s;
        }
        MPI_Allreduce(MPI_IN_PLACE, norm.data(), nst, MPI_DOUBLE, MPI_SUM, comm);
        for(int st = 0;st < nst;st++) norm[st] = 1.0 / norm[st] / get_vel_f();
    }
    std::array<double, 64> f;
    for(int st = 0;st < 64;st++) f[st] = (st < nst) ? wgt[st] * norm[st] : 0.0;
    double *ACC = acc;
    pfor(fs, [=](size_t j) {
        double s = 0.0;
        for(int st = 0;st < nst;st++) s += f[st] * D[(size_t)st*fs + j];
        ACC[j] += s;
    });
    sync();
}

template <typename KpointType> void BatchRho::run(Kpoint<KpointType> **Kpts, double *work)
{
    const int py = hy + ord, pz = hz + ord;
    pfor(fsz, [=, a = acc](size_t j) { a[j] = 0.0; });
    sync();

    for(int kpt = 0;kpt < ct.num_kpts_pe;kpt++)
    {
        double kw = Kpts[kpt]->kp.kweight;
        std::vector<int> sts;
        for(int st = 0;st < Kpts[kpt]->nstates;st++)
            if(kw * Kpts[kpt]->Kstates[st].occupation[0] != 0.0) sts.push_back(st);

        for(size_t s0 = 0;s0 < sts.size();s0 += maxb)
        {
            int nst = (int)std::min((size_t)maxb, sts.size() - s0);
            std::array<double, 64> wgt;
#pragma omp parallel for
            for(int s = 0;s < nst;s++)
            {
                State<KpointType> &S = Kpts[kpt]->Kstates[sts[s0 + s]];
                cf *h = host_stage + (size_t)s*csz;
                const KpointType *src = S.psi;
                for(int ix = 0;ix < hx;ix++)
                    for(int iy = 0;iy < hy;iy++)
                        for(int iz = 0;iz < hz;iz++)
                        {
                            std::complex<double> v = src[((size_t)ix*hy + iy)*hz + iz];
                            h[(size_t)(ix+w)*py*pz + (size_t)(iy+w)*pz + iz + w] = cf{(float)std::real(v), (float)std::imag(v)};
                        }
            }
            for(int s = 0;s < nst;s++) wgt[s] = kw * Kpts[kpt]->Kstates[sts[s0 + s]].occupation[0];
#if SYCL_ENABLED
            if(dev) q->memcpy(c, host_stage, (size_t)nst*csz*sizeof(cf)).wait();
            else
#endif
            std::copy(host_stage, host_stage + (size_t)nst*csz, c);
            batch(nst, wgt, ct.norm_conserving_pp);
        }
    }

#if SYCL_ENABLED
    if(dev) q->memcpy(work, acc, fsz*sizeof(double)).wait();
    else
#endif
    std::copy(acc, acc + fsz, work);
}

}  // namespace

// Fills work (fine grid, this rank) with the density of all k-points of this rank, as the
// GetNewRhoOne loop does. Returns false when the batched path does not apply.
template <typename KpointType> bool GetNewRhoBatch(Kpoint<KpointType> **Kpts, Prolong &P, double *work)
{
    if(!ct.batch_density) return false;
    if(!std::is_same_v<KpointType, std::complex<double>>) return false;
    if(ct.noncoll_factor != 1 || ct.xc_is_meta || ct.prolong_order == 0) return false;
    if(Rmg_G->default_FG_RATIO != 2) return false;
    int cfac = 1;
    if(ct.coalesce_states) cfac = pct.coalesce_factor;
    if(Rmg_G->get_PX0_GRID(1) >= 5) cfac = 1;
    if(cfac != 1) return false;

    static BatchRho *B = nullptr;
    if(!B) B = new BatchRho(P, 16);
    if(!B->ok) return false;
    B->run(Kpts, work);
    return true;
}

template bool GetNewRhoBatch<double>(Kpoint<double> **, Prolong &, double *);
template bool GetNewRhoBatch<std::complex<double>>(Kpoint<std::complex<double>> **, Prolong &, double *);
