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
#include "Kpoint.h"
#include "InputKey.h"
#include "blas.h"
//#include "main.h"
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"




/* Main control structure which is declared extern in main.h so any module */
/* may access it.					                 */
CONTROL ct;

/* PE control structure which is also declared extern in main.h */
PE_CONTROL pct;

KBPSI Kbpsi_str;
unsigned int *perm_ion_index, *perm_state_index, *rev_perm_state_index;
double *projectors, *projectors_x, *projectors_y, *projectors_z;

int *num_nonlocal_ion;
double *kbpsi, *kbpsi_comm, *kbpsi_res, *partial_kbpsi_x, *partial_kbpsi_y, *partial_kbpsi_z;
int kbpsi_num_loop, *kbpsi_comm_send, *kbpsi_comm_recv;
int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;
int *ionidx_allproc;
int max_ion_nonlocal;
STATE *states_tem;
int *state_to_proc;
char *state_overlap_or_not;
STATE *states;
STATE *states1;
ION_ORBIT_OVERLAP *ion_orbit_overlap_region_nl;
double *rho, *rho_old, *rhoc, *vh, *vnuc, *vxc, *rhocore, *eig_rho, *vtot, *vtot_c;
double *rho_oppo, *rho_tot;
int MXLLDA, MXLCOL;
double *sg_twovpsi, *sg_res;
double *statearray, *l_s, *matB, *mat_hb, *mat_X, *Hij, *theta, *work_dis;
double *Hij_00, *Bij_00;
double *work_matrix_row, *coefficient_matrix_row, *nlarray1;
double *work_dis2, *zz_dis, *cc_dis, *gamma_dis, *uu_dis, *mat_Omega;
ORBIT_ORBIT_OVERLAP *orbit_overlap_region;
std::vector<ORBITAL_PAIR> OrbitalPairs;
char *vloc_state_overlap_or_not;
double *orbit_tem;
double *sg_orbit;
double *sg_orbit_res;
int *state_begin;
int *state_end;
double *vtot_global;
double *work_memory;
double *wave_global;
double *rho_global;
double *vxc_old, *vh_old, *vh_corr, *vh_x, *vh_y, *vh_z;


int mpi_nprocs;
int mpi_myrank;

//STATE *states, *states1;

/*Variables from recips.h*/
double b0[3], b1[3], b2[3];
double alat;

std::unordered_map<std::string, InputKey *> ControlMap;

#if ELEMENTAL_LIBS
#include "El.hpp"
using namespace El;
#endif
using namespace std;

MPI_Comm COMM_PEX, COMM_PEY, COMM_PEZ, COMM_3D;
MPI_Comm COMM_EN, COMM_EN1, COMM_EN2;



static void dipole_calculation(double *rhooo, double *dipole);

int main(int argc, char **argv)
{

    // Set branch type
    ct.rmg_branch = RMG_ON;

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

        //  initialize for ELEMENTAl lib
#if ELEMENTAL_LIBS
        Initialize( argc, argv );
#endif

        ReadBranchON(ct.cfile, ct, ControlMap);
        allocate_states();
        get_state_to_proc(states);
        perm_ion_index = new unsigned int[ct.num_ions + 1];
        for(int i = 0; i < ct.num_ions + 1; i++) perm_ion_index[i] = i;

        perm_state_index = new unsigned int[ct.num_states];
        rev_perm_state_index = new unsigned int[ct.num_states];

        switch(ct.runflag)
        {
            case 1:
                if(pct.gridpe == 0) 
                {
                    ReadPermInfo(ct.infile, perm_ion_index);
                }
                MPI_Bcast(perm_ion_index, ct.num_ions, MPI_INT, 0, pct.grid_comm);
                break;
            default:
                if(pct.gridpe == 0) 
                {
                    if(ct.bandwidthreduction)
                        BandwidthReduction(ct.num_ions, ct.ions, perm_ion_index);
                    WritePermInfo(ct.outfile, perm_ion_index);
                }
                MPI_Bcast(perm_ion_index, ct.num_ions, MPI_INT, 0, pct.grid_comm);
                PermAtoms(ct.num_ions, ct.ions, perm_ion_index);
                break;
        }
        ReadOrbitals (ct.cfile, states, ct.ions, pct.img_comm, perm_ion_index);
        GetPermStateIndex(ct.num_ions, ct.ions, perm_ion_index, perm_state_index, rev_perm_state_index);

        init_states();
        my_barrier();




        RmgTimer *RTi = new RmgTimer("1-TOTAL: init");

        init_dimension(&MXLLDA, &MXLCOL);
        init_pe_on();


        if (pct.gridpe == 0)
            printf("\n  MXLLDA: %d ", MXLLDA);

        /* allocate memory for matrixs  */
        allocate_matrix();

        /* Perform some necessary initializations no matter localized or not  
         */
        vxc_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];
        vh_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];

        vh_corr = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)]();
        vh_x = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)]();
        vh_y = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)]();
        vh_z = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)]();


        InitON(vh, rho, rho_oppo, rhocore, rhoc, states, states1, vnuc, vxc, vh_old, vxc_old, ControlMap);


        my_barrier();

        delete(RTi);
        /* Dispatch to the correct driver routine */

        RmgTimer *RTq = new RmgTimer("1-TOTAL: quench");
        switch (ct.forceflag)
        {
            case MD_QUENCH:            /* Quench the electrons */

                quench(states, states1, vxc, vh, vnuc, vh_old, vxc_old, rho, rho_oppo, rhoc, rhocore);
                break;

            default:
                printf("\n undifined MD Method");
                exit(0);
        }

        delete(RTq);

        /* Save data to output file */
        RmgTimer *RTw = new RmgTimer("1-TOTAL: write");

        write_restart(ct.outfile, vh, vxc, vh_old, vxc_old, rho, rho_oppo, &states[0]); 
        double dipole_ion, dipole_ele;
        dipole_calculation(rhoc, &dipole_ion);
        dipole_calculation(rho, &dipole_ele);
        printf("\n duoile    %f %f \n", dipole_ion, dipole_ele);

        /* Save state information to file */
        // write_states_info(ct.outfile, &states[0]);

        my_barrier();
        delete(RTw);

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


    if(pct.imgpe == 0) fclose(ct.logfile);
    RmgPrintTimings(pct.img_comm, ct.logname, ct.scf_steps, pct.num_owned_ions * ct.num_kpts_pe);



    MPI_Finalize();

    return 0;
}


static void dipole_calculation(double *rhooo, double *dipole)
{



    int i, j, k, idx;

    double x, y, z;

    dipole[0] = 0.0;
    dipole[1] = 0.0;
    dipole[2] = 0.0;

    for(i = 0; i < get_FPX0_GRID(); i++)
    {
        x = (get_FPX_OFFSET() + i)*get_hxxgrid() * get_xside();
        for(j = 0; j < get_FPY0_GRID(); j++)
        {
            y = (get_FPY_OFFSET() + j)*get_hyygrid() * get_yside();

            for(k = 0; k < get_FPZ0_GRID(); k++)
            {
                z = (get_FPZ_OFFSET() + k)*get_hzzgrid() * get_zside();

                idx = i * get_FPY0_GRID() * get_FPZ0_GRID() + j*get_FPZ0_GRID() + k;
                dipole[0] += x * rhooo[idx];
                dipole[1] += y * rhooo[idx];
                dipole[2] += z * rhooo[idx];
            }
        }
    }

    dipole[0] *= get_vel_f();
    dipole[1] *= get_vel_f();
    dipole[2] *= get_vel_f();
    
    idx = 3;
    global_sums (dipole, &idx, pct.grid_comm);


}
