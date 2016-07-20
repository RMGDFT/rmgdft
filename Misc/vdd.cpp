//kmlively 7/18/2016 1:00pm 

#include <vector>
#include <cstdlib>
#include <cmath>
#include "transition.h"
#include "common_prototypes.h"
#include "mpi.h"

/*Relative ion index or absolute*/
#define RELATIVE 1

using namespace std;

void vdd(double * rho)
{

    double hxgrid, hygrid, hzgrid, vel, *atomic_rho, distmin, del_vdd;
    double xc, yc, zc;
    int i, j, FPX_OFFSET, FPY_OFFSET, FPZ_OFFSET, FP0_BASIS, dimx, dimy, dimz, min_ion, num_ions;
    int *rel_index, rel_index_size;

    hxgrid = get_hxxgrid(); //The fine grid spacing in x
    hygrid = get_hyygrid(); //The fine grid spacing in y
    hzgrid = get_hzzgrid(); //The fine grid spacing in z

    vel = get_vel_f(); //the volume element of the grid

    dimx = get_FPX0_GRID(); //The number of fine grid points in the x direction for this processor
    dimy = get_FPY0_GRID(); //The number of fine grid points in the y direction for this processor
    dimz = get_FPZ0_GRID(); //The number of fine grid points in the z direction for this processor

    FPX_OFFSET = get_FPX_OFFSET(); //The processor grid offset from pe 0 of this processor in the x direction
    FPY_OFFSET = get_FPY_OFFSET(); //The processor grid offset from pe 0 of this processor in the y direction
    FPZ_OFFSET = get_FPZ_OFFSET(); //The processor grid offset from pe 0 of this processor in the z direction

    FP0_BASIS = get_FP0_BASIS(); //The number of grid points per processor

    num_ions = ct.num_ions;

    atomic_rho = new double[FP0_BASIS]; //dynamically make a new array the size of which is the number of gridpoints in this processor
    lcao_get_rho (atomic_rho); //populates the atomic_rho array with the charge density at each gridpoint due to the atoms in this processor 

    rel_index = new int[num_ions];
    rel_index_size = 0;


    double * loc_array;
    loc_array = new double[num_ions]; //this will be contiguous

    //Initalize loc_array to zeros
    for (int i = 0; i < num_ions; i++)
	loc_array[i] = 0.0;

    double Iondistances[num_ions];//array to store the distance of every Ion from the gridpoint under question in the below loop
    double Gridcoor[3]; //the current gridpoint coordinates




    /* Get the offset of the x coordinate of this processor*/

    xc = FPX_OFFSET * hxgrid;
    //Go through each each grid point by looping through each cartesian direction
    j = 0; //index for use in all three dimensions

    for (int ix = 0; ix < dimx; ix++)
    {

	yc = FPY_OFFSET * hygrid; //Get the offset of the y coordinate of this processor
	for (int iy = 0; iy < dimy; iy++)
	{

	    zc = FPZ_OFFSET * hzgrid; //Get the offset of the z coordinate of this processor
	    for (int iz = 0; iz < dimz; iz++)
	    {
		//populate gridcoor array with this gridpoint's coordinates
		Gridcoor[0] = xc;
		Gridcoor[1] = yc;
		Gridcoor[2] = zc;

		del_vdd = (rho[j] - atomic_rho[j]); //the portion of the voronoi density deformation integral for this gridpoint 

		//populate the Iondistances vector with the distance of every ION / ION image from this gridpoint
		for (int ionct = 0; ionct < ct.num_ions; ionct++) //for each ion number
		{
		    //populate Iondistances vector with the distance between this gridpoint and every ION and ION image
		    //minvoroimage returns the smallest distance between these gridpoint coordinates and every image of the ion in the first input.
		    //this ion has images due to periodic conditions, minvoroimage finds the smallest distance between this gridpoint and these images
		    Iondistances[ionct] = minvoroimage(&ct.ions[ionct], Gridcoor); 

		}

		//find the smallest distance in Iondistances to find the Ion nearest this gridpoint

		distmin = Iondistances[0]; //set distmin initially to the first distance

		min_ion = 0; //minimum ion number also set initially to the first ion number


		for (int i=1; i <  ct.num_ions ; i++) //go through the Iondistances vector to just before the last entry
		{
		    if (Iondistances[i] < distmin) //if the next entry is less than the current, set the next entry to distmin
		    {
			distmin = Iondistances[i]; 
			min_ion = i; //record the position in Iondistances of the minimum distance
			//this should correspond to the ion_number in ct.ions[ion_number]
		    }
		}


#if RELATIVE
		/*Check if index is already recorded*/
		for (i = 0; i < rel_index_size; i++) 
		{
		    if (rel_index[i] == min_ion)
		    {
			loc_array[i] += del_vdd;
			break;
		    }
		}

		/*If not*/
		if (i == rel_index_size)
		{
		    rel_index_size += 1;
		    rel_index[i] = min_ion; 
		    loc_array[i] += del_vdd;
		}
#else

		loc_array[min_ion] += del_vdd;
#endif


		zc += hzgrid; //advance to the next gridpoint in the z direction
		j++; //incrememnt 3-d index


	    }
	    yc += hygrid; //advance to the next gridpoint in the y direction
	}
	xc += hxgrid; //advance to the next gridpoint in the x direction
    }

#if RELATIVE
    printf("\n PE %d: Found %d voronoi domains, vel %e", pct.imgpe, rel_index_size, vel);
		
    for (i = 0; i < rel_index_size; i++) 
	printf("\n PE %d: Domain %d Ion %d Partial integral %e", pct.imgpe, i, rel_index[i], vel*loc_array[i]);
#else

    printf("\n PE %d: Partial integrals are %e %e %e, vel is %e", pct.imgpe, vel*loc_array[0], vel*loc_array[1], vel*loc_array[2], vel);

#endif

#if 0
    MPI_Status Stat;
    MPI_Request reqs;

    //send loc_array to processor 0 for summation
    if (pct.imgpe == 0)
    {	

	double *vdd_final;
	vdd_final = static_cast<double*>(calloc(ct.num_ions, sizeof(double)));
	for (int i = 0; i < ct.num_ions; i++)
	    vdd_final[i] = 0.0;

	MPI_Recv(&loc_array, ct.num_ions, MPI_DOUBLE, MPI_ANY_SOURCE, pct.imgpe, MPI_COMM_WORLD, &Stat);

	for (int i = 0; i < ct.num_ions; i++)
	    vdd_final[i] += loc_array[i];

	std::cout << "vdd_final is " << std::endl;
	for (int i = 0; i < ct.num_ions; i++)
	    std::cout << "ION " << i << ": " << vdd_final[i] << std::endl;
    }
    else
	MPI_Isend(&loc_array, ct.num_ions, MPI_DOUBLE, 0, pct.imgpe, MPI_COMM_WORLD, &reqs); //all other processors send loc_array

    if (pct.imgpe == 0)
    {

    }

#endif

    delete [] atomic_rho;
    delete [] loc_array;
    delete [] rel_index;

}
