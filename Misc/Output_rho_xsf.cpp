#include "portability.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include "transition.h"
#if !(defined(_WIN32) || defined(_WIN64))
#include "cnpy.h"
#endif
#include "main.h"



void Output_rho_xsf(double *array_3d, MPI_Comm comm)
{

#if !(defined(_WIN32) || defined(_WIN64))

    char newname[MAX_PATH + 20];



    int sizes[3], subsizes[3], starts[3];
    MPI_Info fileinfo;
    MPI_Status status;
    MPI_Offset disp = 0;
    MPI_Datatype num_as_string;
    MPI_Datatype localarray;
    int FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET, idx; 
    int idx1;
    char *const fmt="%12.6e ";
    char *const endfmt="%12.6e\n";
    const int charspernum=13;

    int ion;
    SPECIES *sp;
    ION *iptr;


    /* this datatype describes the mapping of the local array
     * to the global array (file)
     * */

    sizes[0] = Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO);
    sizes[1] = Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO);
    sizes[2] = Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO);

    subsizes[0] = Rmg_G->get_PX0_GRID(Rmg_G->default_FG_RATIO);
    subsizes[1] = Rmg_G->get_PY0_GRID(Rmg_G->default_FG_RATIO);
    subsizes[2] = Rmg_G->get_PZ0_GRID(Rmg_G->default_FG_RATIO);

    //Rmg_G->find_node_offsets(Rmg_G->get_rank(), Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO), 
    Rmg_G->find_node_offsets(pct.gridpe, Rmg_G->get_NX_GRID(Rmg_G->default_FG_RATIO), 
            Rmg_G->get_NY_GRID(Rmg_G->default_FG_RATIO), Rmg_G->get_NZ_GRID(Rmg_G->default_FG_RATIO),
            &FPX_OFFSET, &FPY_OFFSET, &FPZ_OFFSET);

    starts[0] = FPX_OFFSET;
    starts[1] = FPY_OFFSET;
    starts[2] = FPZ_OFFSET;


    /* each number is represented by charspernum chars */
    MPI_Type_contiguous(charspernum, MPI_CHAR, &num_as_string); 
    MPI_Type_commit(&num_as_string); 

    char *array_in_char;

    array_in_char = new char[sizes[0] * sizes[1] * sizes[2] *charspernum];

    int count = 0;
    for (int k=0; k<subsizes[2]; k++) {
        for (int j=0; j<subsizes[1]; j++) {
            for (int i=0; i<subsizes[0]; i++) {
                idx = i * subsizes[1] * subsizes[2] + j * subsizes[2] + k;
                if( (((i+starts[0] + 1)%5) == 0) | ((i +starts[0]) == sizes[0] -1) )
                    sprintf(&array_in_char[count*charspernum], endfmt, array_3d[idx]);
                else
                    sprintf(&array_in_char[count*charspernum], fmt, array_3d[idx]);
                count++;
            }
        }
    }


    MPI_Barrier(comm);



    int order = MPI_ORDER_FORTRAN;
    //int order = MPI_ORDER_C;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, order, num_as_string, &localarray);
    MPI_Type_commit(&localarray);


    snprintf(newname, MAX_PATH, "%s%s", pct.image_path[pct.thisimg], "rho.xsf");

    if(pct.gridpe == 0)
    {
        FILE *fhand = fopen(newname, "w");


        fprintf(fhand, "\n#This is a file to be viewed by xcrysden \n");

        fprintf(fhand, "\n CRYSTAL \n");
        fprintf(fhand, "\n# lattice parameters in Angstroms \n");
        fprintf(fhand, "\n PRIMVEC \n");
        fprintf(fhand, "%12.6f %12.6f %12.6f\n",Rmg_L.get_a0(0) * a0_A, Rmg_L.get_a0(1) * a0_A, Rmg_L.get_a0(2) * a0_A);
        fprintf(fhand, "%12.6f %12.6f %12.6f\n",Rmg_L.get_a1(0) * a0_A, Rmg_L.get_a1(1) * a0_A, Rmg_L.get_a1(2) * a0_A);
        fprintf(fhand, "%12.6f %12.6f %12.6f\n",Rmg_L.get_a2(0) * a0_A, Rmg_L.get_a2(1) * a0_A, Rmg_L.get_a2(2) * a0_A);

        fprintf(fhand, "\n# Atomic coordinates in Angstroms \n");
        fprintf(fhand, "\n PRIMCOORD");
        fprintf(fhand, "\n  %d %d\n", ct.num_ions, 1);
        

        for(ion = 0; ion < ct.num_ions; ion++)
        {
            iptr = &ct.ions[ion];
            sp = &ct.sp[iptr->species];
            fprintf(fhand, " %2s  %18.6e   %18.6e   %18.6e\n", sp->atomic_symbol, iptr->crds[0] * a0_A, 
                  iptr->crds[1] * a0_A,   iptr->crds[2] * a0_A);   
        }

        fprintf(fhand, "BEGIN_BLOCK_DATAGRID_3D\n");                        
        fprintf(fhand, "   charge_density_view_with_Xcrysdena\n");
        fprintf(fhand, "   BEGIN_DATAGRID_3D_this_is_3Dgrid#1\n");           
        fprintf(fhand, "%6d %6d %6d\n", sizes[0], sizes[1], sizes[2]);
        fprintf(fhand, "%8.3f %8.3f %8.3f\n", 0.0, 0.0, 0.0);
        fprintf(fhand, "%12.6f %12.6f %12.6f\n",Rmg_L.get_a0(0) * a0_A, Rmg_L.get_a0(1) * a0_A, Rmg_L.get_a0(2) * a0_A);
        fprintf(fhand, "%12.6f %12.6f %12.6f\n",Rmg_L.get_a1(0) * a0_A, Rmg_L.get_a1(1) * a0_A, Rmg_L.get_a1(2) * a0_A);
        fprintf(fhand, "%12.6f %12.6f %12.6f\n",Rmg_L.get_a2(0) * a0_A, Rmg_L.get_a2(1) * a0_A, Rmg_L.get_a2(2) * a0_A);
        fclose(fhand);
    }
    MPI_Barrier(comm);

    MPI_File mpi_fhand ;
    MPI_Info_create(&fileinfo);

    //int amode = MPI_MODE_CREATE|MPI_MODE_WRONLY;
    int amode = MPI_MODE_APPEND|MPI_MODE_WRONLY;
    MPI_File_open(comm, newname, amode, fileinfo, &mpi_fhand);

    
    MPI_File_seek(mpi_fhand, 0, MPI_SEEK_END);
    MPI_File_get_position(mpi_fhand, &disp);

    MPI_File_set_view(mpi_fhand, disp,  MPI_CHAR, localarray, "native", MPI_INFO_NULL);
     
    MPI_File_write_all(mpi_fhand, array_in_char, subsizes[0] * subsizes[1] * subsizes[2], num_as_string, &status);

    MPI_Type_free(&localarray);
    MPI_Type_free(&num_as_string);

    MPI_File_close(&mpi_fhand);

    if(pct.gridpe == 0)
    {
        FILE *fhand = fopen(newname, "a");
        fprintf(fhand, "   END_DATAGRID_3D\n");                    
        fprintf(fhand, " END_BLOCK_DATAGRID_3D\n");
        fclose(fhand);
    }

    snprintf(newname, MAX_PATH, "%s%s", pct.image_path[pct.thisimg], "electron_density.npz");
    delete [] array_in_char;
    double cell[9];
    for(int i = 0; i < 3; i++) 
    {  
        cell[i] = Rmg_L.a0[i];
        cell[i+3] = Rmg_L.a1[i];
        cell[i+6] = Rmg_L.a2[i];
    }
    const unsigned int shape[] = {3,3};
    cnpy::npz_save(newname, "cell", cell, shape, 2, "w");

    int length = sizes[0] * sizes[1] * sizes[2];
    double *g_array = new double[sizes[0] * sizes[1] * sizes[2]];
    for(idx = 0; idx < sizes[0] * sizes[1] * sizes[2]; idx++) g_array[idx] = 0.0;
    for(int i = 0; i < subsizes[0]; i++)
    for(int j = 0; j < subsizes[1]; j++)
    for(int k = 0; k < subsizes[2]; k++)
    { 

      idx = i * subsizes[1] * subsizes[2] + j * subsizes[2] + k;
      idx1 = (i+starts[0]) * sizes[1] * sizes[2] + (j+ starts[1]) * sizes[2] + k + starts[2];
      g_array[idx1] = array_3d[idx];
    }

    MPI_Allreduce(MPI_IN_PLACE, g_array, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);
 
    const unsigned int shape1[] = {(unsigned int)sizes[0], (unsigned int)sizes[1], (unsigned int)sizes[2]};
    cnpy::npz_save(newname, "density", g_array, shape1, 3, "a");
#endif

}

