/************************** SVN Revision Information **************************
 **    $Id$    **
 ******************************************************************************/

/*
get_Hvnlij:

Get the elements of the Hamiltonian matrix due to the non-local
potential, and add them into Aij.


 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "make_conf.h"
#include "params.h"

#include "rmgtypedefs.h"
#include "typedefs.h"
#include "RmgTimer.h"
#include "transition.h"
#include "prototypes_on.h"
#include "init_var.h"

#include "my_scalapack.h"
#include "blas.h"
#include "Kbpsi.h"



void KbpsiComm()
{
    int loop, nh, ion, ip1, ip2, st1, st2, ist;
    MPI_Status mstatus;
    MPI_Request request;
    int num_proj, num_orbital_thision;
    int send_size, recv_size, position;
    char *send_buff, *recv_buff;
    int ion1, ion2, ion1_global, ion2_global;
    int iip1, iip2, iip1a, iip2a;
    int size, proc, proc1, proc2;
    std::vector<double> kbpsi_recv;
    std::vector<int> orbital_index;
    unsigned int idx;


    send_buff = new char[Kbpsi_str.max_send_size];
    recv_buff = new char[Kbpsi_str.max_recv_size];


    for (loop = 0; loop < Kbpsi_str.kbpsi_comm_loop; loop++)
    {

        proc1 = Kbpsi_str.comm_info[loop].send_to_pe; 
        proc2 = Kbpsi_str.comm_info[loop].recv_from_pe; 

        if(proc1 >=0)
        {

            send_size = Kbpsi_str.comm_info[loop].send_size;
            position = 0;
            for(idx=0; idx < Kbpsi_str.comm_info[loop].send_ions.size(); idx++)
            {
                ion = Kbpsi_str.comm_info[loop].send_ions[idx];
                num_proj = pct.prj_per_ion[pct.ionidx[ion]];
                num_orbital_thision = Kbpsi_str.num_orbital_thision[ion];

                MPI_Pack(&num_orbital_thision, 1, MPI_INT, 
                        send_buff, send_size, &position, pct.grid_comm);

                //  vector of int, length of num_orbital_thision
                MPI_Pack(Kbpsi_str.orbital_index[ion].data(),
                        num_orbital_thision, MPI_INT, 
                        send_buff, send_size, &position, pct.grid_comm);
                //  vector of <beta|phi>, 
                //length of num_orbital_thisoon *num_proj
                MPI_Pack(Kbpsi_str.kbpsi_ion[ion].data(), 
                        num_orbital_thision * num_proj, MPI_DOUBLE, 
                        send_buff, send_size, &position, pct.grid_comm);
            }


            MPI_Isend(send_buff, send_size, MPI_BYTE, proc1, loop, pct.grid_comm, &request);
        }
        if(proc2 >=0)
        {
            recv_size = Kbpsi_str.comm_info[loop].recv_size;
            MPI_Recv(recv_buff, recv_size, MPI_BYTE, proc2, loop, pct.grid_comm, &mstatus);

            position = 0;
            for(idx=0; idx < Kbpsi_str.comm_info[loop].recv_ions.size(); idx++)
            {
                ion = Kbpsi_str.comm_info[loop].recv_ions[idx];
                num_proj = pct.prj_per_ion[pct.ionidx[ion]];


                MPI_Unpack(recv_buff, recv_size, &position, 
                        &num_orbital_thision, 1, MPI_INT, pct.grid_comm); 


                orbital_index.resize(num_orbital_thision);
                int num_kbpsi = num_orbital_thision *num_proj;
                kbpsi_recv.resize(num_kbpsi);

                //  vector of int, length of num_orbital_thision
                MPI_Unpack(recv_buff, recv_size, &position, orbital_index.data(),
                        num_orbital_thision, MPI_INT, pct.grid_comm);
                Kbpsi_str.orbital_index[ion].insert(Kbpsi_str.orbital_index[ion].end(),  
                        orbital_index.begin(), orbital_index.end());

                
                //  vector of <beta|phi>, 
                //length of num_orbital_thisoon *num_proj
                MPI_Unpack(recv_buff, recv_size, &position, kbpsi_recv.data(),
                        num_kbpsi, MPI_DOUBLE, pct.grid_comm);

                Kbpsi_str.kbpsi_ion[ion].insert(Kbpsi_str.kbpsi_ion[ion].end(),
                        kbpsi_recv.begin(), kbpsi_recv.end());
            }

        }
        if(proc1 >=0) MPI_Wait(&request, &mstatus);
    }

}

