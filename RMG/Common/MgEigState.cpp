

#include "TradeImages.h"
#include "FiniteDiff.h"
#include "Mgrid.h"
#include "BlasWrappers.h"
#include "auxiliary.h"
#include "const.h"
#include "rmgtypedefs.h"
#include "typedefs.h"
#include "common_prototypes.h"
#include "common_prototypes1.h"
#include "rmg_alloc.h"
#include "rmg_error.h"
#include "RmgTimer.h"


extern STATE *states;

using namespace std;

static std::mutex vtot_sync_mutex;

template <typename RmgType>
void MgEigState (STATE * sp, int tid, rmg_double_t * vtot_psi)
{

    RmgTimer RT("Mg_eig");

    int idx, cycles, P0_BASIS;
    int nits, pbasis, sbasis;
    rmg_double_t eig, diag, t1, t2, t3, t4;
    rmg_double_t *work1, *nv, *ns, *res2;
    rmg_double_t *tmp_psi, *res, *saved_psi;
    rmg_double_t *nvtot_psi;
    RmgType *tmp_psi_t, *work2_t, *res_t, *res2_t, *sg_psi_t;
    int eig_pre[6] = { 0, 3, 6, 2, 2, 2 };
    int eig_post[6] = { 0, 3, 6, 2, 2, 2 };
    int ione = 1;
    int dimx, dimy, dimz, levels, potential_acceleration;
    rmg_double_t hxgrid, hygrid, hzgrid, sb_step;
    rmg_double_t tarr[8];
    RmgType *sg_twovpsi_t, *work1_t;
    BaseGrid G;
    TradeImages T;
    Mgrid MG;


    P0_BASIS = G.get_P0_BASIS();

    nits = ct.eig_parm.gl_pre + ct.eig_parm.gl_pst;
    dimx = sp->dimx;
    dimy = sp->dimy;
    dimz = sp->dimz;
    hxgrid = sp->hxgrid;
    hygrid = sp->hygrid;
    hzgrid = sp->hzgrid;
    levels = ct.eig_parm.levels;
    if ((ct.runflag == 0) && (ct.scf_steps < 2)) {
        levels = 0;
    }

    sb_step = 1.0;
    pbasis = sp->pbasis;
    sbasis = sp->sbasis;

    /* Grab some memory */
    res2_t = new RmgType[sbasis];
    work2_t = new RmgType[4 * sbasis];
    work1_t = new RmgType[4 * sbasis];
    work1 = new rmg_double_t[4 * sbasis];

    sg_psi_t = new RmgType[sbasis];
    res = new rmg_double_t[sbasis];
    sg_twovpsi_t = new RmgType[sbasis];

#if !BATCH_NLS
    ns = new rmg_double_t[sbasis];
    nv = new rmg_double_t[sbasis];
#endif
    res2 = new rmg_double_t[sbasis];
    saved_psi = new rmg_double_t[sbasis];
    nvtot_psi = new rmg_double_t[sbasis];
    tmp_psi_t = new RmgType[sbasis];
    res_t = new RmgType[sbasis];

    tmp_psi = sp->psiR;


    if(ct.eig_parm.mucycles > 1)
        mix_betaxpsi1(sp);

    /* Get the non-local operator and S acting on psi (nv and ns, respeget_vel()y) */
#if !BATCH_NLS
    app_nls (tmp_psi, NULL, nv, NULL, ns, NULL, pct.oldsintR_local, NULL, sp->istate, sp->kidx);

#else
    nv = &pct.nv[sp->istate * P0_BASIS]; 
    ns = &pct.ns[sp->istate * P0_BASIS]; 
#endif


    // Copy double precision ns into temp single precision array */
    for(idx = 0;idx < pbasis;idx++) {
        work1_t[idx] = (RmgType)ns[idx];
    }
  

    /*Apply double precision Mehrstellen right hand operator to ns and save in res2 */
    CPP_app_cir_driver<RmgType> (work1_t, res2_t, dimx, dimy, dimz, ct.kohn_sham_fd_order);

    // Copy double precision psi into single precison array
    for(idx = 0;idx < pbasis;idx++) {
        tmp_psi_t[idx] = (RmgType)tmp_psi[idx];
    }

    // Setup some potential acceleration stuff
    potential_acceleration = ((ct.potential_acceleration_constant_step > 0.0) || (ct.potential_acceleration_poisson_step > 0.0));
    if(potential_acceleration) {
        for(idx = 0;idx <P0_BASIS;idx++) {
            nvtot_psi[idx] = vtot_psi[idx];
            saved_psi[idx] = tmp_psi[idx];
        }
    }




    /* Smoothing cycles */
    for (cycles = 0; cycles <= nits; cycles++)
    {


        /* Apply Mehrstellen left hand operators */
        diag = CPP_app_cil_driver<RmgType> (tmp_psi_t, work2_t, dimx, dimy, dimz, hxgrid, hygrid, hzgrid, ct.kohn_sham_fd_order);
        diag = -1.0 / diag;


        // Copy saved application to ns to res
        QMD_copy(pbasis, res2_t, 1, res_t, 1);


        if(potential_acceleration) {
            /* Generate 2 * V * psi */
            CPP_genvpsi (tmp_psi_t, sg_twovpsi_t, nvtot_psi, nv, NULL, 0.0, dimx, dimy, dimz);
        }
        else { 
            CPP_genvpsi (tmp_psi_t, sg_twovpsi_t, vtot_psi, nv, NULL, 0.0, dimx, dimy, dimz);
        }



        /* B operating on 2*V*psi stored in work1 */
        CPP_app_cir_driver<RmgType> (sg_twovpsi_t, work1_t, dimx, dimy, dimz, ct.kohn_sham_fd_order);
        for(idx = 0; idx < dimx * dimy * dimz; idx++) work1_t[idx] += TWO * nv[idx];


        RmgType tscal = -ONE;
        QMD_axpy (pbasis, tscal, work2_t, ione, work1_t, ione);

        /* If this is the first time through compute the eigenvalue */
        if ((cycles == 0) || (potential_acceleration != 0)) 
        {

            eig = 0.0;
            t2 = 0.0;

            for (idx = 0; idx < pbasis; idx++)
            {

                t2 += (rmg_double_t)tmp_psi_t[idx] * (rmg_double_t)res_t[idx];
                eig += (rmg_double_t)tmp_psi_t[idx] * (rmg_double_t)work1_t[idx];

            }

            tarr[0] = t2;
            tarr[1] = eig;
            idx = 2;


            /*If diagonalization is done every step, do not calculate eigenvalues, use those
             * from diagonalization, except for the first step, since at that time eigenvalues 
	     * are not defined yet*/
            if ((ct.diag == 1) && (potential_acceleration == 0) && (ct.scf_steps < ct.end_diag))
            {
                if (ct.scf_steps == 0)
                {
                    global_sums (tarr, &idx, pct.grid_comm);
                    eig = tarr[1] / (TWO * tarr[0]);
                    sp->eig[0] = eig;
                    sp->oldeig[0] = eig;
                }
		else
                    eig = sp->eig[0];
	    }
            else
            {
                global_sums (tarr, &idx, pct.grid_comm);
                eig = tarr[1] / (TWO * tarr[0]);
                sp->eig[0] = eig;
                if(ct.scf_steps == 0) {
                    sp->oldeig[0] = eig;
                }
            }
            
            if(potential_acceleration) {
                t1 = eig;
                eig = 0.3 * eig + 0.7 * sp->oldeig[0];
                sp->oldeig[0] = t1;
            }

        }

        /* Now either smooth the wavefunction or do a multigrid cycle */
        if (cycles == ct.eig_parm.gl_pre)
        {


            t1 = TWO * eig;
            for (idx = 0; idx <P0_BASIS; idx++)
            {

                res_t[idx] = t1 * res_t[idx] - work1_t[idx];

            }

            /* Pack the residual data into multigrid array */
            CPP_pack_ptos<RmgType> (sg_psi_t, res_t, dimx, dimy, dimz);

            T.trade_images<RmgType> (sg_psi_t, dimx, dimy, dimz, FULL_TRADE);


            /* Smooth it once and store the smoothed residual in work1 */
            t1 = 145.0;

            CPP_app_smooth<RmgType> (sg_psi_t, work1_t, G.get_PX0_GRID(), G.get_PY0_GRID(), G.get_PZ0_GRID());


            if(potential_acceleration) {
                t1 = eig - states[0].eig[0];
                t1 = -t1*t1 / 10.0;
            }
            else {
                t1 = 0.0;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            MG.mgrid_solv (sg_twovpsi_t, work1_t, work2_t,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, get_neighbors(), levels, eig_pre, eig_post, 1, sb_step, t1,
                        G.get_NX_GRID(), G.get_NY_GRID(), G.get_NZ_GRID(),
                        G.get_PX_OFFSET(), G.get_PY_OFFSET(), G.get_PZ_OFFSET(),
                        G.get_PX0_GRID(), G.get_PY0_GRID(), G.get_PZ0_GRID(), ct.boundaryflag);



            /* The correction is in a smoothing grid so we use this
             * routine to update the orbital which is stored in a physical grid.
             */

            t1 = -1.0;
            CPP_pack_stop_axpy<RmgType> (sg_twovpsi_t, tmp_psi_t, t1, dimx, dimy, dimz);


        }
        else
        {


            t1 = TWO * eig;
            t2 = ZERO;
            t4 = ct.eig_parm.gl_step * diag;

            for (idx = 0; idx <P0_BASIS; idx++)
            {

                t3 = t1 * (rmg_double_t)res_t[idx] - (rmg_double_t)work1_t[idx];
                t2 += t3 * t3;
                tmp_psi_t[idx] += t4 * t3;

            }

            if (cycles == 0)
            {

                t2 = real_sum_all (t2, pct.grid_comm);
                t1 = (rmg_double_t) (ct.psi_nbasis);
                sp->res = ct.hmaxgrid * ct.hmaxgrid * sqrt (t2 / t1) * 0.25;

            }

        }

    }                           /* end for */


    if(potential_acceleration) {

        // Save potential used for this orbital and update potential for future orbitals
        for(idx = 0;idx <P0_BASIS;idx++) {
            sp->dvhxc[idx] = nvtot_psi[idx];
        }


        if(ct.potential_acceleration_constant_step > 0.0) {

            t1 = 1.8 * ct.potential_acceleration_constant_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;

            vtot_sync_mutex.lock();
            for(idx = 0;idx <P0_BASIS;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * PI * sp->occupation[0] * tmp_psi_t[idx] * (tmp_psi_t[idx] - saved_psi[idx]);
            }
            vtot_sync_mutex.unlock();

        }

        if(ct.potential_acceleration_poisson_step > 0.0) {

            // construct delta_rho
            for(idx = 0;idx <P0_BASIS;idx++) {
                res_t[idx] = -4.0 * PI * sp->occupation[0] *
                           (tmp_psi_t[idx] - saved_psi[idx]) * (2.0*saved_psi[idx] + (tmp_psi_t[idx] - saved_psi[idx]));
            }

            // zero out solution vector
            for(idx = 0;idx <P0_BASIS;idx++) {
                sg_twovpsi_t[idx] = 0.0;
            }

            /* Pack delta_rho into multigrid array */
            CPP_pack_ptos<RmgType> (sg_psi_t, res_t, dimx, dimy, dimz);
            T.trade_images<RmgType> (sg_psi_t, dimx, dimy, dimz, FULL_TRADE);
            /* Smooth it once and store the smoothed charge in res */
            CPP_app_smooth1<RmgType> (sg_psi_t, res_t, G.get_PX0_GRID(), G.get_PY0_GRID(), G.get_PZ0_GRID());

            // neutralize cell with a constant background charge
            t2 = 0.0;
            for(idx = 0;idx <P0_BASIS;idx++) {
                t2 += res_t[idx];
            }
            t2 = real_sum_all(t2, pct.grid_comm) / (G.get_NX_GRID() * G.get_NY_GRID() * G.get_NZ_GRID());
            for(idx = 0;idx <P0_BASIS;idx++) {
                res_t[idx] -= t2;
            }

            /* Do multigrid step with solution returned in sg_twovpsi */
            eig_pre[0] = 2;
            eig_post[0] = 2;
            levels=1;
            MG.mgrid_solv (sg_twovpsi_t, res_t, work2_t,
                        dimx, dimy, dimz, hxgrid,
                        hygrid, hzgrid, 0, G.get_neighbors(), levels, eig_pre, eig_post, 1, 1.0, 0.0,
                        G.get_NX_GRID(), G.get_NY_GRID(), G.get_NZ_GRID(),
                        G.get_PX_OFFSET(), G.get_PY_OFFSET(), G.get_PZ_OFFSET(),
                        G.get_PX0_GRID(), G.get_PY0_GRID(), G.get_PZ0_GRID(), ct.boundaryflag);

            for(idx = 0;idx <P0_BASIS;idx++) {
                res_t[idx] = 0.0;
            }
            CPP_pack_stop_axpy<RmgType> (sg_twovpsi_t, res_t, 1.0, dimx, dimy, dimz);
            t1 = ct.potential_acceleration_poisson_step;
            if(sp->occupation[0] < 0.5) t1 = 0.0;
            vtot_sync_mutex.lock();
            for(idx = 0;idx <P0_BASIS;idx++) {
               vtot_psi[idx] = vtot_psi[idx] + t1 * res_t[idx];
            }
            vtot_sync_mutex.unlock();
        }


    } // end if

    // Copy single precision orbital back to double precision
    for(idx = 0;idx < pbasis;idx++) {

        tmp_psi[idx] = (rmg_double_t)tmp_psi_t[idx];

    }

    /* Release our memory */
    delete [] res_t;
    delete [] tmp_psi_t;
    delete [] nvtot_psi;
    delete [] saved_psi;
    delete [] res2;

#if !BATCH_NLS
    delete [] nv;
    delete [] ns;
#endif

    delete [] sg_twovpsi_t;
    delete [] res;
    delete [] sg_psi_t;

    delete [] work1;
    delete [] work1_t;
    delete [] work2_t;
    delete [] res2_t;




} // end MgEigState

extern "C" void mg_eig_state(STATE *sp, int tid, rmg_double_t *vtot_psi)
{
    MgEigState<rmg_double_t> (sp, tid, vtot_psi);
}
extern "C" void mg_eig_state_f(STATE *sp, int tid, rmg_double_t *vtot_psi)
{
    MgEigState<rmg_float_t> (sp, tid, vtot_psi);
}

extern "C" void mg_eig_state_driver (STATE * sp, int tid, rmg_double_t * vtot_psi, int precision)
{

#if !GAMMA_PT

    // Single precision code path not programmed for non gamma calculations
//    mg_eig_state (sp, tid, vtot_psi);

#else

    if(precision == sizeof(rmg_double_t)) {

        MgEigState<rmg_double_t> (sp, tid, vtot_psi);

    }
    else {

        MgEigState<rmg_float_t> (sp, tid, vtot_psi);

    }

#endif

}

