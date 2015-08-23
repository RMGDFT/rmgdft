/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*

   init_nonlocal_comm.c

   Sets up the parameters for communication in calculation of
   <psi_i|kb_ion> <kb_ion|psi_j>

   each processor calculates the part of <psi|kb> in which psi are stored in
   this processor.


   This should be called after get_nlop.c
 */




#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "make_conf.h"
#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "prototypes_on.h"
#include "Kbpsi.h"
#include "init_var.h"



void InitNonlocalComm(void)
{
    int ion, idx, item, nh;
    MPI_Status mstatus;
    int size;
    int send_size, recv_size, size_perstate;
    int proc, proc1, proc2, st1, num_proj;
    int ion1, ion2, ion1_global, ion2_global;
    KBPSI_ION tem_kbpsi_ion;
    std::vector<double> tem_vector;
    MPI_Request request;
    KBPSI_COMM_INFO kbpsi_oneloop;
    Kbpsi_str.kbpsi_comm_loop = kbpsi_num_loop;

    Kbpsi_str.kbpsi_ion = new std::vector<KBPSI_ION>[pct.n_ion_center];
    Kbpsi_str.comm_info = new KBPSI_COMM_INFO[Kbpsi_str.kbpsi_comm_loop];
    

    proc = pct.gridpe;
    for(ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
    {
        ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
        num_proj = pct.prj_per_ion[ion1_global];
        tem_vector.resize(num_proj);

        for (st1 = ct.state_begin; st1 < ct.state_end; st1++)
        {
            idx = (st1 - ct.state_begin) * ct.num_ions + ion1_global;
            if (ion_orbit_overlap_region_nl[idx].flag == 1)
            {
                tem_kbpsi_ion.orbital_index = st1;
                tem_kbpsi_ion.kbpsi = tem_vector;
                Kbpsi_str.kbpsi_ion[ion1].emplace_back(tem_kbpsi_ion);
            }
        }
    }


    for (idx = 0; idx < Kbpsi_str.kbpsi_comm_loop; idx++)
    {

        proc1 = kbpsi_comm_send[idx];
        proc2 = kbpsi_comm_recv[idx];

        kbpsi_oneloop.send_to_pe = proc1;
        kbpsi_oneloop.recv_from_pe = proc2;

        send_size = 0;
        for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
            for (ion2 = 0; ion2 < num_nonlocal_ion[proc1]; ion2++)
            {
                ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                if (ion1_global == ion2_global)
                {
                    kbpsi_oneloop.send_ions.emplace_back(ion1);

                    num_proj = pct.prj_per_ion[ion1_global];
                    size_perstate = sizeof(int) + num_proj * sizeof(double);
                    send_size += Kbpsi_str.kbpsi_ion[ion1].size() * size_perstate;
                }
            }

        for (ion1 = 0; ion1 < num_nonlocal_ion[proc]; ion1++)
            for (ion2 = 0; ion2 < num_nonlocal_ion[proc2]; ion2++)
            {
                ion1_global = ionidx_allproc[proc * max_ion_nonlocal + ion1];
                ion2_global = ionidx_allproc[proc2 * max_ion_nonlocal + ion2];

                if (ion1_global == ion2_global)
                {
                    kbpsi_oneloop.recv_ions.emplace_back(ion1);

                }
            }



        recv_size = 0;
        if(proc1 >=0)
            MPI_Isend(&send_size, 1, MPI_INT, proc1, idx, pct.grid_comm, &request);
        if(proc2 >=0)
            MPI_Recv(&recv_size, 1, MPI_INT, proc2, idx, pct.grid_comm, &mstatus);
        if(proc1 >=0) MPI_Wait(&request, &mstatus);

        kbpsi_oneloop.send_size = send_size;
        kbpsi_oneloop.recv_size = recv_size;
        Kbpsi_str.comm_info[idx] = kbpsi_oneloop;

    }

}

