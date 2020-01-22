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

    T *array_glob = new T[sizes[0] * sizes[1] * sizes[2]]();

    for (int i=0; i<subsizes[0]; i++) {
        for (int j=0; j<subsizes[1]; j++) {
            for (int k=0; k<subsizes[2]; k++) {
                int idx = i * subsizes[1] * subsizes[2] + j * subsizes[2] + k;
                int idx_g = (i+starts[0]) * sizes[1] * sizes[2] + 
                    (j+starts[1]) * sizes[2] + k + starts[2];
                array_glob[idx_g] = array_dist[idx];
            }
        }
    }

    int length = sizes[0] * sizes[1] * sizes[2] * sizeof(T)/sizeof(double);
    MPI_Allreduce(MPI_IN_PLACE, array_glob, length, MPI_DOUBLE, MPI_SUM, pct.grid_comm);

    if(pct.gridpe == 0)
    {
        FILE *fhand = fopen(filename.c_str(), "w");


        fprintf(fhand, "#This is a Cube file to be viewed \n");
        fprintf(fhand, "# Atomic coordinates in Angstroms \n");
        int nval = sizeof(T)/sizeof(double);
        fprintf(fhand, " %d  0.0  0.0  0.0 %d \n", (int)Atoms.size(), nval);

        fprintf(fhand, "%d %12.6f %12.6f %12.6f\n",sizes[0], 
                Rmg_L.get_a0(0) * a0_A/sizes[0], 
                Rmg_L.get_a0(1) * a0_A/sizes[0], 
                Rmg_L.get_a0(2) * a0_A/sizes[0]);
        fprintf(fhand, "%d %12.6f %12.6f %12.6f\n",sizes[1], 
                Rmg_L.get_a1(0) * a0_A/sizes[1], 
                Rmg_L.get_a1(1) * a0_A/sizes[1], 
                Rmg_L.get_a1(2) * a0_A/sizes[1]);
        fprintf(fhand, "%d %12.6f %12.6f %12.6f\n",sizes[2], 
                Rmg_L.get_a2(0) * a0_A/sizes[2], 
                Rmg_L.get_a2(1) * a0_A/sizes[2], 
                Rmg_L.get_a2(2) * a0_A/sizes[2]);

        for(ion = 0; ion < ct.num_ions; ion++)
        {
            iptr = &Atoms[ion];
            sp = &Species[iptr->species];
            int a_num = GetAtomicNumber(sp->atomic_symbol);
            fprintf(fhand, " %d  %f %12.6e   %12.6e   %12.6e\n", a_num,
                    iptr->partial_charge, iptr->crds[0] * a0_A, 
                    iptr->crds[1] * a0_A,   iptr->crds[2] * a0_A);   
        }

        for (int i=0; i<sizes[0]; i++) {
            for (int j=0; j<sizes[1]; j++) {
                for (int k=0; k<sizes[2]; k++) {
                    int idx = i * sizes[1] * sizes[2] + j * sizes[2] + k;
                    fprintf(fhand, " %g ", std::real(array_glob[idx]));
                    if(nval == 2) 
                        fprintf(fhand, " %g ", std::imag(array_glob[idx]));
                    if(k * nval % 6 == 5) fprintf(fhand, "\n");
                }

                fprintf(fhand, "\n");
            }
        }

        fclose(fhand);

    }


}

