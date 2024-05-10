#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include "transition.h"
#include "MapElements.h"


template void OutputCubeFile(double *array_3d, int grid, std::string filename);
template void OutputCubeFile(std::complex<double> *array_3d, int grid, std::string filename);
template <typename T> void OutputCubeFile(T *array_dist, int grid, std::string filename)
{

    RmgTimer RT0("CubeFile output");

    int sizes[3], subsizes[3], starts[3];
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET; 

    int ion;
    SPECIES *sp;
    ION *iptr;


    /* this datatype describes the mapping of the local array
     * to the global array (file)
     * */

    sizes[0] = Rmg_G->get_NX_GRID(grid);
    sizes[1] = Rmg_G->get_NY_GRID(grid);
    sizes[2] = Rmg_G->get_NZ_GRID(grid);

    subsizes[0] = Rmg_G->get_PX0_GRID(grid);
    subsizes[1] = Rmg_G->get_PY0_GRID(grid);
    subsizes[2] = Rmg_G->get_PZ0_GRID(grid);

    //Rmg_G->find_node_offsets(Rmg_G->get_rank(), Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO), 
    Rmg_G->find_node_offsets(pct.gridpe, sizes[0], sizes[1], sizes[2],
            &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);

    starts[0] = FPX_OFFSET;
    starts[1] = FPY_OFFSET;
    starts[2] = FPZ_OFFSET;
    int pex, pey, pez;
    pe2xyz (pct.gridpe, &pex, &pey, &pez);

    double *array_d = (double *)array_dist;
    // number of doubles in the last dimension
    // number of lines in the output
    // numbe of char for the last dimension
    int nval = sizeof(T)/sizeof(double);
    int num_d = subsizes[2] * nval;
    // each data has 12 char %12.3e format
    int char_per_data = 12;
    std::vector<int> num_char_pz(pct.pe_z, 0);
    std::vector<int> num_start_pz(pct.pe_z, 0);
    int num_char = num_d * char_per_data;

    for(int iz = 0; iz < num_d; iz++)
    {
        int iz_glob = iz + FPZ_OFFSET * nval;
        if(iz_glob % 6 == 5)  num_char++;
    }

    // new line for each set of z data
    if(pez == pct.pe_z -1) num_char++;

    num_char_pz[pez] = num_char;

    MPI_Allreduce(MPI_IN_PLACE, num_char_pz.data(), pct.pe_z, MPI_INT, MPI_SUM, pct.grid_comm);


    num_start_pz[0] = 0;
    int tot_num_char;
    tot_num_char = num_char_pz[0]/pct.pe_x/pct.pe_y;

    for(int ipz = 1; ipz < pct.pe_z; ipz++)
    {
        num_start_pz[ipz] = num_start_pz[ipz-1] + num_char_pz[ipz-1]/pct.pe_x/pct.pe_y;
        tot_num_char += num_char_pz[ipz]/pct.pe_x/pct.pe_y;
    }


    char *array_print = new char[subsizes[0] * subsizes[1] * num_char];

    if(pct.gridpe == 0 )
    {
        FILE *fhand = fopen(filename.c_str(), "w");


        fprintf(fhand, "#This is a Cube file to be viewed \n");
        fprintf(fhand, "# Atomic coordinates lattice in bohr \n");
        fprintf(fhand, " %d  0.0  0.0  0.0 %d \n", (int)Atoms.size(), nval);

        fprintf(fhand, "%d %12.6f %12.6f %12.6f\n", sizes[0], 
                Rmg_L.get_a0(0) /sizes[0], 
                Rmg_L.get_a0(1) /sizes[0], 
                Rmg_L.get_a0(2) /sizes[0]);
        fprintf(fhand, "%d %12.6f %12.6f %12.6f\n", sizes[1], 
                Rmg_L.get_a1(0) /sizes[1], 
                Rmg_L.get_a1(1) /sizes[1], 
                Rmg_L.get_a1(2) /sizes[1]);
        fprintf(fhand, "%d %12.6f %12.6f %12.6f\n", sizes[2], 
                Rmg_L.get_a2(0) /sizes[2], 
                Rmg_L.get_a2(1) /sizes[2], 
                Rmg_L.get_a2(2) /sizes[2]);

        for(ion = 0; ion < ct.num_ions; ion++)
        {
            iptr = &Atoms[ion];
            sp = &Species[iptr->species];
            int a_num = GetAtomicNumber(sp->atomic_symbol);
            fprintf(fhand, " %d  %f %12.6e   %12.6e   %12.6e\n", a_num,
                    iptr->partial_charge, iptr->crds[0] , 
                    iptr->crds[1] ,   iptr->crds[2] );   
        }
        fclose(fhand);
    }

    char *buffer = array_print;
    for (int i=0; i<subsizes[0]; i++) {
        for (int j=0; j<subsizes[1]; j++) {
            for (int k=0; k<subsizes[2] * nval; k++) {
                int idx = i * subsizes[1] * subsizes[2] *nval + j * subsizes[2] *nval + k;
                int k_glob = k + FPZ_OFFSET * nval; 
                if(k_glob%6 == 5)
                {
                    sprintf(buffer,  "%12.3e\n", array_d[idx]);
                    buffer += 13;
                }
                else
                {
                    sprintf(buffer,  "%12.3e", array_d[idx]);
                    buffer += 12;
                }
            }
            if(pez == pct.pe_z -1 )
            {
                sprintf(buffer,  "\n");
                buffer +=1;
            }


        }
    }

    subsizes[2] = num_char;
    sizes[2] = tot_num_char;
    starts[2] = num_start_pz[pez];
    int order = MPI_ORDER_C;


    MPI_Info fileinfo;
    MPI_Datatype grid_char;
    MPI_Status status;
    MPI_Offset disp;



    MPI_Type_create_subarray(3, sizes, subsizes, starts, order, MPI_CHAR, &grid_char);
    MPI_Type_commit(&grid_char);

    MPI_Info_create(&fileinfo);

    MPI_File mpi_fhand ;

    int pbasis_char = subsizes[0] * subsizes[1] * num_char;
    MPI_Barrier(pct.grid_comm);
    int amode = MPI_MODE_APPEND|MPI_MODE_WRONLY;
    MPI_File_open(pct.grid_comm, filename.c_str(), amode, fileinfo, &mpi_fhand);
    MPI_File_seek(mpi_fhand, 0, MPI_SEEK_END);
    MPI_File_get_position(mpi_fhand, &disp);
    MPI_File_set_view(mpi_fhand, disp, MPI_CHAR, grid_char, "native", MPI_INFO_NULL);
    MPI_File_write_all(mpi_fhand, array_print, pbasis_char, MPI_CHAR, &status);
    MPI_File_close(&mpi_fhand);

}

