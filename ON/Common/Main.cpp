/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/****f* QMD-MGDFT/md.c *****
 * NAME
 *   Ab initio O(n) real space code with 
 *   localized orbitals and multigrid acceleration
 *   Version: 3.0.0
 * COPYRIGHT
 *   Copyright (C) 2001  Wenchang Lu,
 *                       Jerzy Bernholc
 * FUNCTION
 *   int main(int argc, char **argv)
 *   Main program
 *   Read-in all informations, structures, pseudopentials, etc. 
 *   Then enters the main driver loop. 
 * INPUTS
 *   when we run it, we need to give the input control 
 *   file name in the first argument
 *   for example, md in.diamond8
 * OUTPUT
 *   nothing
 * PARENTS
 *   This is grand-grand-....
 * CHILDREN
 *   run.c
 * SEE ALSO
 *   main.h for structure defination
 * SOURCE
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <unordered_map>
#include "const.h"
#include "RmgTimer.h"
#include "RmgException.h"
#include "rmgtypedefs.h"
#include "params.h"
#include "typedefs.h"
#include "rmg_error.h"
#include "InputKey.h"
#include "blas.h"
//#include "main.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"
#include "rmgthreads.h"
#include "rmg_gemm.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "LocalObject.h"
#include "LdaU_on.h"
#include "Exxbase.h"
#include "Exx_on.h"

#if QMCPACK_SUPPORT
    #include "WriteEshdf.h"
#endif


LocalObject<double> *LocalOrbital;
LocalObject<double> *H_LocalOrbital;
LocalObject<double> *LocalProj;
LocalObject<double> *LocalAtomicOrbital;
LdaU_on *ldaU_on;
Exx_on<double> *Exx_onscf;

double *Kbpsi_mat;
double *Kbpsi_mat_local;

/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;


std::vector<int> num_nonlocal_ion;
std::vector<int> ionidx_allproc;
int max_ion_nonlocal;
STATE *states;
double *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot, *vtot_c;
double *rho_oppo, *rho_tot;
int MXLLDA, MXLCOL;
double *sg_twovpsi, *sg_res;
double *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *Hij_00, *Bij_00;
double *work_matrix_row, *coefficient_matrix_row, *nlarray1;
double *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
int *state_begin;
int *state_end;
double *vxc_old, *vh_old, *vh_corr;


int mpi_nprocs;
int mpi_myrank;

//STATE *states, *states1;

std::vector<ION> Atoms;
std::vector<SPECIES> Species;


/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;

std::unordered_map<std::string, InputKey *> ControlMap;

using namespace std;

MPI_Comm COMM_PEX, COMM_PEY, COMM_PEZ, COMM_3D;
MPI_Comm COMM_EN, COMM_EN1, COMM_EN2;



int main(int argc, char **argv)
{

    // Set branch type
    ct.rmg_branch = RMG_ON;
    ct.save_args(argc, argv);

    FiniteDiff::allocation_limit = 0;

    // double a[4], b[4];
    // MPI_Init_thread(&argc, &argv, ct.mpi_threadlevel, &provided);
    //MPI_Init(&argc, &argv);

    //DiagElemental(4, a, b);    


    ct.mpi_threadlevel = MPI_THREAD_SERIALIZED;
    ct.mpi_threadlevel = 0;
    // mpi::Initialize(argc, argv);
    RmgTimer *RT = new RmgTimer("1-TOTAL");

    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;
    // tem();

    ct.images_per_node = 1;
    ct.proj_nophase = 1;
    try 
    {
        InitIo(argc, argv, ControlMap);
        double A[1], B[1],C[1];
        rmg::gemm("N", "N", 1, 1, 1, 0.5, A, 1, B, 1, 0.0, C, 1);


        ReadBranchON(ct.cfile, ct, ControlMap);
        WriteXyz(ct.cfile);
        states = init_states();
        ReadOrbitals (ct.cfile, states, Atoms, pct.img_comm);

        MPI_Barrier(pct.img_comm);

        RmgTimer *RTi = new RmgTimer("1-TOTAL: init");

        init_dimension(&MXLLDA, &MXLCOL);

        if (pct.gridpe == 0)
            rmg::printlog("\n  MXLLDA: %d ", MXLLDA);

        /* allocate memory for matrixs  */
        allocate_matrix();

        /* Perform some necessary initializations no matter localized or not  
         */
        vxc_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];
        vh_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];

        vh_corr = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)]();


        InitON(vh, rho, rho_oppo, rhocore, rhoc, states, vnuc, vxc, vh_old, vxc_old, ControlMap);

        MPI_Barrier(pct.img_comm);

        delete(RTi);
        /* Dispatch to the correct driver routine */

        RmgTimer *RTq = new RmgTimer("1-TOTAL: quench");
        switch (ct.forceflag)
        {
            case MD_QUENCH:            /* Quench the electrons */
            case BAND_STRUCTURE:

                quench(states, vxc, vh, vnuc, vh_old, vxc_old, rho, rho_oppo, rhoc, rhocore);
                break;
            case TDDFT:

                if(!ct.restart_tddft) 
                {
                    quench(states, vxc, vh, vnuc, vh_old, vxc_old, rho, rho_oppo, rhoc, rhocore);
                    write_restart(ct.outfile, vh, vxc, vh_old, vxc_old, rho, rho_oppo, &states[0]);

                }

                OnTddft (vxc, vh, vnuc, rho, rho_oppo, rhocore, rhoc, *LocalOrbital, 
                        *H_LocalOrbital, *LocalProj);
                break;
            case Exx_only:
                {    
                    int nstates_occ = 0;
                    std::vector<double> occs;
                    // calculate num of occupied states
                    for(int i=0;i < ct.num_states;i++) 
                    {
                        if(states[i].occupation[0] > 1.0e-6)
                        {
                            occs.push_back(states[i].occupation[0]);
                            nstates_occ++;
                        }
                        else
                        {
                            break;
                        }
                    }

                    double *psi = new double[Rmg_G->get_P0_BASIS(1) * ct.num_states];
                    Exxbase<double> Exx(*Rmg_G, *Rmg_halfgrid, Rmg_L, "tempwave", nstates_occ, occs.data(), psi, ct.exx_mode);
                    if(ct.exx_mode == EXX_DIST_FFT)
                        Exx.ReadWfsFromSingleFile();

                    Exx.Vexx_integrals(ct.exx_int_file);
                }

                break;


            default:
                rmg::printlog("\n undifined MD Method");
                exit(0);
        }

        delete(RTq);

        /* Save data to output file */
        RmgTimer *RTw = new RmgTimer("1-TOTAL: write");

        write_restart(ct.outfile, vh, vxc, vh_old, vxc_old, rho, rho_oppo, &states[0]); 

        RmgTimer *RTO = new RmgTimer("WriteOrbitals");
        LocalOrbital->WriteOrbitalsToSingleFiles(ct.outfile, *Rmg_G);
        if(ct.num_ions == 1)
        {
            int pbasis = Rmg_G->get_P0_BASIS(1);
            int num_orb = LocalOrbital->num_tot;
            if(num_orb != LocalOrbital->num_thispe)
            {
                rmg::printlog("Main.cpp:  num_tot %d != num_thispe %d", num_orb, LocalOrbital->num_thispe);
                exit(0);
            }
            double *Cij_glob = new double[num_orb * num_orb];
            mat_dist_to_global(zz_dis, pct.desca, Cij_glob);

            double one = 1.0, zero = 0.0;
            dgemm("N", "N", &pbasis, &num_orb, &num_orb , &one, LocalOrbital->storage_cpu, &pbasis, 
                    Cij_glob, &num_orb, &zero, H_LocalOrbital->storage_cpu, &pbasis);

            for(int idx = 0; idx < num_orb * pbasis; idx++) 
                LocalOrbital->storage_cpu[idx] = H_LocalOrbital->storage_cpu[idx];
            LocalOrbital->WriteOrbitalsToSingleFiles(ct.outfile, *Rmg_G);
        }
        delete RTO;

        // test conditions
        check_tests();

        if(ct.write_qmcpack_restart)
        {

            double *eig = new double[LocalOrbital->num_tot];
            double *occ = new double[LocalOrbital->num_tot];
            for(int is = 0; is < LocalOrbital->num_tot; is++)
            {
                eig[is] = states[is].eig[0];
                occ[is] = states[is].occupation[0];
            }
            std::string fname(ct.outfile);
            WriteWavefunctions(fname, *LocalOrbital, zz_dis, *Rmg_G, eig, occ);
            MPI_Barrier(MPI_COMM_WORLD);
            delete [] eig;
            delete [] occ;
#if QMCPACK_SUPPORT
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if(rank == 0)
            {
                WriteQmcpackRestart(fname);
            }
#else
            rmg::printlog ("Unable to write QMCPACK file since RMG was not built with HDF and QMCPACK support.\n");
#endif

        }

        /* Save state information to file */

        if(ct.write_qmcpack_restart_localized)
        {

#if QMCPACK_SUPPORT
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            std::string fname(ct.outfile);
            //  double *Cij_glob = new double[num_orb * num_orb];
            //  mat_dist_to_global(zz_dis, pct.desca, Cij_glob);

            WriteCij(fname, zz_dis);

            // if(pct.gridpe == 0)
            //{
            //   std::string cijname = fname + "_spin" + std::to_string(pct.spinpe) + "_Cij";
            //  int fhand = open(cijname.c_str(), O_CREAT | O_TRUNC | O_RDWR, S_IREAD | S_IWRITE);
            // if (fhand < 0) {
            //    rmg::printlog("Can't open restart file %s", cijname.c_str());
            //   rmg::error("Terminating.");
            //  }
            // 
            // size_t wsize = ct.num_states * ct.num_states * sizeof(double);
            // write (fhand, Cij_glob, wsize);
            // close(fhand);
            // fflush(NULL);
            //  }

            MPI_Barrier(MPI_COMM_WORLD);
            if(rank == 0)
            {
                WriteQmcpackRestartLocalized(fname);
            }
#else
            rmg::printlog ("Unable to write QMCPACK file since RMG was not built with HDF and QMCPACK support.\n");
#endif

        }

        /* Save state information to file */
        // write_states_info(ct.outfile, &states[0]);

        delete(RTw);
        MPI_Barrier(pct.img_comm);

    }
    // Catch exceptions issued by us.
    catch(RmgFatalException const &e) {
        std::cout << e.rwhat() << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // By std
    catch (std::exception &e) {
        std::cout << "Caught a std exception: " << e.what () << std::endl;
        MPI_Finalize();
        exit(0);
    } 

    // Catchall
    catch (...) {
        std::cout << "Caught an unknown exception of some type." << std::endl;
        MPI_Finalize();
        exit(0);
    } 

    delete(RT);

    check_tests();
    if(pct.imgpe == 0) fclose(ct.logfile);
    int override_rank = 0;
    if(pct.imgpe==0) MPI_Comm_rank (pct.img_comm, &override_rank);
    RmgPrintTimings(pct.img_comm, ct.logname, ct.scf_steps, pct.num_owned_ions * ct.num_kpts_pe, override_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    RmgTerminateThreads();

}
