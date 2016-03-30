/************************** SVN Revision Information **************************
 **    $Id: main.c 2671 2014-10-13 19:50:57Z luw $    **
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

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "svnrev.h"
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
#include "init_var.h"
#include "transition.h"
#include "prototypes_on.h"
#include "Kbpsi.h"


#include "../Headers/common_prototypes.h"
#include "../Headers/common_prototypes1.h"
#include "prototypes_tddft.h"



#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <string>
#include <boost/algorithm/string.hpp>

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
char *state_overlap_or_not;
int *send_to, *recv_from, num_sendrecv_loop;
int *send_to1, *recv_from1, num_sendrecv_loop1;
int *ionidx_allproc;
int max_ion_nonlocal;
int NPES;
STATE *states_tem;
int *state_to_proc;
STATE *states;
STATE *states1;
STATE *states_distribute;
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
double *vxc_old, *vh_old;


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
int FilenameIncrement1(char *pathname);

int main(int argc, char **argv)
{


    double *Hmatrix, *Smatrix, *Xij_00, *Yij_00, *Zij_00;
    double *Xij_dist, *Yij_dist, *Zij_dist;
    double *rho_matrix_local, *rho_matrix;
    double *Hmatrix_t0, *Hmatrix_old;

    int Ieldyn=1, iprint = 0;
    int n2, numst,i;
    double *Pn0, *Pn1;
    double dipole_ion[3], dipole_ele[3];
    int ione=1;
    FILE *dfi;
    int tot_steps, pre_steps;
    int amode, fhand;
    int FP0_BASIS;
    amode = S_IREAD | S_IWRITE;

    char filename[MAX_PATH+200];
    char char_tem[MAX_PATH+200];
    int ntem, inner_step;
    double t2, tem, half=0.5;




    ct.mpi_threadlevel = MPI_THREAD_SERIALIZED;
    ct.mpi_threadlevel = 0;

    RmgTimer *RT = new RmgTimer("1-TOTAL");
    /* Define a default output stream, gets redefined to log file later */
    ct.logfile = stdout;

    RmgTimer *RT1 = new RmgTimer("1-TOTAL: read_and_alloc");
    ct.images_per_node = 1;
    InitIo(argc, argv, ControlMap);

    //  initialize for ELEMENTAl lib
#if ELEMENTAL_LIBS
    Initialize( argc, argv );
#endif

    ReadBranchON(ct.cfile, ct, ControlMap);
    allocate_states();
    get_state_to_proc(states);
    perm_ion_index = (unsigned int *) malloc(ct.num_ions * sizeof(int));
    for(int i = 0; i < ct.num_ions; i++) perm_ion_index[i] = i;

    perm_state_index = (unsigned int *) malloc(ct.num_states * sizeof(int));
    rev_perm_state_index = (unsigned int *) malloc(ct.num_states * sizeof(int));

    if(pct.gridpe == 0) 
    {
        ReadPermInfo(ct.infile, perm_ion_index);
    }
    MPI_Bcast(perm_ion_index, ct.num_ions, MPI_INT, 0, pct.grid_comm);
    ReadOrbitals (ct.cfile, states, ct.ions, pct.img_comm, perm_ion_index);
    GetPermStateIndex(ct.num_ions, ct.ions, perm_ion_index, perm_state_index, rev_perm_state_index);

    init_states();



    delete(RT1);

    RmgTimer *RTi = new RmgTimer("1-TOTAL: init");

    init_dimension(&MXLLDA, &MXLCOL);
    init_pe_on();


    /* allocate memory for matrixs  */
    allocate_matrix();

    /* Perform some necessary initializations no matter localized or not  
     */
    vxc_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];
    vh_old = new double[Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO)];


    InitON(vh, rho, rho_oppo, rhocore, rhoc, states, states1, vnuc, vxc,
vh_old, vxc_old, ControlMap);


    delete(RTi);


    numst = ct.num_states;
    FP0_BASIS = Rmg_G->get_P0_BASIS(Rmg_G->default_FG_RATIO);

    RmgTimer *RT2 = new RmgTimer("1-TOTAL: orbital_comm");
    orbital_comm(states);
    delete(RT2);

    RmgTimer *RT3 = new RmgTimer("1-TOTAL: kbpsi");
    KbpsiUpdate(states);
    delete(RT3);

    RmgTimer *RT4 = new RmgTimer("1-TOTAL: GetHS");
    GetHS(states, states1, vtot_c, Hij_00, Bij_00);
    delete(RT4);

    RmgTimer *RT5 = new RmgTimer("1-TOTAL: others");
    Cpdgemr2d(numst, numst, Hij_00, ione, ione, pct.descb, Hij, ione, ione,
            pct.desca, pct.desca[1]);
    n2 = MXLLDA * MXLCOL;
    dcopy(&n2, Hij, &ione, matB, &ione);
    PDTRAN(&numst, &numst, &half, matB, &ione, &ione, pct.desca,
                &half, Hij, &ione, &ione, pct.desca);

    Cpdgemr2d(numst, numst, Bij_00, ione, ione, pct.descb, matB, ione, ione,
            pct.desca, pct.desca[1]);


    get_dm_diag_p(states, matB, mat_X, Hij);

    write_eigs(states);

    n2 = (ct.state_end - ct.state_begin) * ct.num_states;
    Xij_00 = new double[n2];
    Yij_00 = new double[n2];
    Zij_00 = new double[n2];
    rho_matrix = new double[n2];


    n2 = MXLLDA * MXLCOL;
    Xij_dist = new double[n2];
    Yij_dist = new double[n2];
    Zij_dist = new double[n2];

    /*  Xij = <phi|x|phi>, Yij = <phi|y|phi>, Zij = <phi|z|phi>  */ 

    get_phi_xyz_phi(states, Xij_00, Yij_00, Zij_00);

    Cpdgemr2d(numst, numst, Xij_00, ione, ione, pct.descb, Xij_dist, ione, ione,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Yij_00, ione, ione, pct.descb, Yij_dist, ione, ione,
            pct.desca, pct.desca[1]);
    Cpdgemr2d(numst, numst, Zij_00, ione, ione, pct.descb, Zij_dist, ione, ione,
            pct.desca, pct.desca[1]);

    n2 = numst * numst;
    Pn0 = new double[2*n2];
    Pn1 = new double[2*n2];
    Hmatrix = new double[n2];
    Hmatrix_t0 = new double[n2];
    Hmatrix_old = new double[n2];
    Smatrix = new double[n2];

    /*  matB: overlap matrix, Hij:  Hamiltonian matrix  distributed in
     *  scalapack way*/

    dipole_calculation(rhoc, dipole_ion);
    dipole_calculation(rho, dipole_ele);



    rmg_printf("\n  x dipolll  %f %f", dipole_ion[0], dipole_ele[0]);
    rmg_printf("\n  y dipolll  %f %f", dipole_ion[1], dipole_ele[1]);
    rmg_printf("\n  z dipolll  %f %f", dipole_ion[2], dipole_ele[2]);

    Cpdgemr2d(numst, numst, mat_X, ione, ione, pct.desca, rho_matrix, ione, ione,
            pct.descb, pct.desca[1]);

    InitStatedistribute ();
    rho_matrix_local = new double[pct.num_local_orbit *pct.num_local_orbit];

    mat_dist_to_global(mat_X, Pn0, pct.desca);
    MatrixToLocal(states_distribute, Pn0, rho_matrix_local);
    for (i = 0; i< FP0_BASIS; i++) rho[i] = 0.0;
    GetNewRhoLocal (states_distribute, rho, rho_matrix_local, rho_matrix);

    dipole_calculation(rho, dipole_ele);
    rmg_printf("\n  x dipoll3  %f %f", dipole_ion[0], dipole_ele[0]);
    rmg_printf("\n  y dipoll3  %f %f", dipole_ion[1], dipole_ele[1]);
    rmg_printf("\n  z dipoll3  %f %f", dipole_ion[2], dipole_ele[2]);

#if 0
    DiagScalapack(states, ct.num_states, Hij_00, Bij_00, rho_matrix, theta);
    GetNewRho_on(states, rho, rho_matrix);
    write_eigs(states);
    dipole_calculation(rho, dipole_ele);
    rmg_printf("\n  x dipoll3  %f %f", dipole_ion[0], dipole_ele[0]);
    rmg_printf("\n  y dipoll3  %f %f", dipole_ion[1], dipole_ele[1]);
    rmg_printf("\n  z dipoll3  %f %f", dipole_ion[2], dipole_ele[2]);

    read_rhomatrix(ct.infile, rho_matrix);
    GetNewRho_on(states, rho, rho_matrix);
    dipole_calculation(rho, dipole_ele);
    rmg_printf("\n  x dipoll4  %f %f", dipole_ion[0], dipole_ele[0]);
    rmg_printf("\n  y dipoll4  %f %f", dipole_ion[1], dipole_ele[1]);
    rmg_printf("\n  z dipoll4  %f %f", dipole_ion[2], dipole_ele[2]);

#endif

    mat_dist_to_global(mat_X, Pn0, pct.desca);
    mat_dist_to_global(matB, Smatrix, pct.desca);

    for(i = 0; i < n2; i++) Pn0[i+n2] = 0.0;



    double time_step = 0.2;
    double efield[3];
    efield[0] = ct.x_field_0 * ct.e_field;
    efield[1] = ct.y_field_0 * ct.e_field;
    efield[2] = ct.z_field_0 * ct.e_field;


    delete(RT5);

    double fs= 0.02418884;  // 1fs = 0.02418884 *10^-15 second 
    if(pct.gridpe == 0)
    {
        sprintf(char_tem, "%s%s", pct.image_path[pct.thisimg], "dipole.dat");
        int name_incr = FilenameIncrement1(char_tem);
        sprintf(filename, "%s.%02d", char_tem, name_incr);

        dfi = fopen(filename, "w");
    }

    if(pct.gridpe == 0)fprintf(dfi, "\n  electric field:  %f  %f  %f ",efield[0], efield[1], efield[2]);

    pre_steps = 0;
    tot_steps = 0;

    if(ct.runflag == Restart_TDDFT) 
    {
        sprintf(filename, "%s%s", ct.infile, ".TDDFT_restart");
        fhand = open(filename, O_RDWR);
        if (fhand < 0)
            rmg_error_handler(__FILE__, __LINE__, " Unable to write file ");
        read(fhand, &pre_steps,  sizeof(int));
        read(fhand, Pn0, 2* numst * numst*sizeof(double));
        close(fhand);
    }

    mat_dist_to_global(Hij, Hmatrix_old, pct.desca);
    if(pre_steps == 0 ) 
    {
        for(i = 0; i < MXLLDA * MXLCOL; i++) Hij[i] = Hij[i] 
            + (efield[0] * Xij_dist[i] + efield[1] * Yij_dist[i] +
                    efield[2] * Zij_dist[i])/time_step *2.0;
    }

    mat_dist_to_global(Hij, Hmatrix, pct.desca);

    for(ct.scf_steps = 0; ct.scf_steps < ct.max_scf_steps; ct.scf_steps++)
    {
        tot_steps = ct.scf_steps + pre_steps;

        dcopy(&n2, Hmatrix, &ione, Hmatrix_t0, &ione);


        for(inner_step = 0; inner_step < 1; inner_step++)
        {


            for(i = 0; i < n2; i++) Hmatrix[i] = 0.5 * time_step*(Hmatrix[i] + Hmatrix_t0[i]);


            RmgTimer *RT2a = new RmgTimer("1-TOTAL: ELDYN");
            eldyn_(&numst, Smatrix, Hmatrix, Pn0, Pn1, &Ieldyn, &iprint);
            delete(RT2a);

            //for(i = 0; i < n2; i++) mat_X[i]= Pn1[i];

            RmgTimer *RT4a = new RmgTimer("1-TOTAL: GetNewRho");
            MatrixToLocal(states_distribute, Pn1, rho_matrix_local);
            GetNewRhoLocal (states_distribute, rho, rho_matrix_local, rho_matrix);
            delete(RT4a);

            int iii = get_FP0_BASIS();

            double tcharge = 0.0;
            for (i = 0; i < get_FP0_BASIS(); i++)
                tcharge += rho[i];
            ct.tcharge = real_sum_all(tcharge, pct.grid_comm);
            ct.tcharge = real_sum_all(ct.tcharge, pct.spin_comm);


            ct.tcharge *= get_vel_f();

            t2 = ct.nel / ct.tcharge;
            dscal(&iii, &t2, &rho[0], &ione);


            if(fabs(t2 -1.0) > 1.0e-11 && pct.gridpe == 0)
                printf("\n Warning: total charge Normalization constant = %e  \n", t2-1.0);

            RmgTimer *RT4b = new RmgTimer("1-TOTAL: Updatepot");
            UpdatePot(vxc, vh, vxc_old, vh_old, vnuc, rho, rho_oppo, rhoc, rhocore);
            delete(RT4b);

            RmgTimer *RT4c = new RmgTimer("1-TOTAL: GetHS");
            for (i = 0; i < get_FP0_BASIS(); i++) vtot[i] = vh[i] +vxc[i] - vh_old[i] - vxc_old[i];
            
            get_vtot_psi(vtot_c, vtot, Rmg_G->default_FG_RATIO);

            HijUpdateNCpp (states_distribute, vtot_c, rho_matrix_local, Hmatrix);

            for(i = 0; i < n2; i++) Hmatrix[i] += Hmatrix_old[i];

            //GetHS(states, states1, vtot_c, Hij_00, Bij_00);

            //Cpdgemr2d(numst, numst, Hij_00, ione, ione, pct.descb, Hij, ione, ione,
             //       pct.desca, pct.desca[1]);
            //mat_dist_to_global(Hij, Hmatrix, pct.desca);
            delete(RT4c);

            tem = 0;
            for(i = 0; i < n2; i++) tem += (Hmatrix[i] - Hmatrix_old[i])*(Hmatrix[i] - Hmatrix_old[i]);
            if(pct.gridpe == 0) printf("\n step = %d, innerstep = %d, Hij - Hij_old = %e\n", tot_steps, inner_step, tem);
            dcopy(&n2, Hmatrix, &ione, Hmatrix_old, &ione);
        }

       // for(i = 0; i < 2*n2; i++) Pn0[i]= Pn1[i]*t2;
        for(i = 0; i < 2*n2; i++) Pn0[i]= Pn1[i];
        dipole_calculation(rho, dipole_ele);

        dipole_ele[0] -= dipole_ion[0];
        dipole_ele[1] -= dipole_ion[1];
        dipole_ele[2] -= dipole_ion[2];


        if(pct.gridpe == 0)fprintf(dfi, "\n  %f  %18.10f  %18.10f  %18.10f ",
                tot_steps*time_step, dipole_ele[0], dipole_ele[1], dipole_ele[2]);
        if((ct.scf_steps +1)%ct.checkpoint == 0)
        {
            write_data(ct.outfile, vh, vxc, vh_old, vxc_old, rho, &states[0]); 
            if(pct.gridpe == 0)
            {
                sprintf(filename, "%s%s", ct.outfile, ".TDDFT_restart");
                fhand = open(filename, O_CREAT |O_TRUNC |O_RDWR, amode);
                if (fhand < 0)
                    rmg_error_handler(__FILE__, __LINE__, " Unable to write file ");

                write(fhand, &tot_steps,  sizeof(int));
                write(fhand, Pn0, 2* numst * numst*sizeof(double));
                close(fhand);
            }

        }


    }



    tot_steps = ct.scf_steps + pre_steps;
    write_data(ct.outfile, vh, vxc, vh_old, vxc_old, rho, &states[0]); 
    if(pct.gridpe == 0)
    {
        sprintf(filename, "%s%s", ct.outfile, ".TDDFT_restart");
        fhand = open(filename, O_CREAT |O_TRUNC |O_RDWR, amode);
        if (fhand < 0)
            rmg_error_handler(__FILE__, __LINE__, " Unable to write file ");

        write(fhand, &tot_steps,  sizeof(int));
        write(fhand, Pn0, 2* numst * numst*sizeof(double));
        close(fhand);
    }

    //    get_dm_diag_p(states, matB, mat_X, Hij);

    //    write_eigs(states);

    delete(RT);
    if(pct.imgpe == 0) fclose(ct.logfile);
    RmgPrintTimings(Rmg_G, ct.logname, ct.scf_steps);

    MPI_Finalize();

    return 0;
}


int FilenameIncrement1(char *pathname)
{

    int lognum = 0;

    boost::filesystem::path current_path(pathname);
    std::string dirname  = current_path.parent_path().string(); 
    std::string basename = boost::filesystem::basename(pathname);
    if(!dirname.length()) dirname = dirname + "./";
    // Does parent path exist?
    if( boost::filesystem::exists(dirname)) {

        // yes so check if it's a file 
        if(!boost::filesystem::is_directory(dirname)) {
            throw RmgFatalException() << "Found " << dirname << "  that is not a directory in " __FILE__ << " at line " << __LINE__ << ".\n";
        }

    }
    else {

        // no so need to make it
        if(!boost::filesystem::create_directory(dirname)) {
            throw RmgFatalException() << "Unable to create logfile directory " << dirname << " in " << __FILE__ << " at line " << __LINE__ << ".\n";
        }

    }


    char lognum_str[4]; 
    snprintf(lognum_str, 3, "%02d", lognum);
    std::string nextname = std::string(pathname) + "." + lognum_str;
    while (boost::filesystem::exists(nextname))
    {
        if (++lognum > 99)
            throw RmgFatalException() << "You have over 100 logfiles, you need to think of a better job naming scenario!" << "\n";
        nextname.erase();
        snprintf(lognum_str, 3, "%02d", lognum);
        nextname = std::string(pathname) + "." + lognum_str;
    }

    // return lognum
    return lognum;

}


